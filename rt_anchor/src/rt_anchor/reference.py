"""The bundled reference dataset — the column every result is expressed on.

The v2 method calibrates one column *onto another*. Something has to be the
"other", and for a user who has only their own two runs that something must
ship with the package. ``reference_data/col35/`` holds the Column 25 pair
(internally column ``25``): the serum run that defines the ``Cal_RT_min`` time
axis, and the standards run that defines the iRT ruler.

Two consequences worth stating out loud:

* Numbers from different reference pairs are **not comparable**. Substituting
  your own reference is supported and sometimes right (a different lab, a
  different gradient you want everything expressed on), but it moves the axis,
  so :func:`resolve_reference` records exactly which pair was used and the model
  JSON carries that record.
* Both halves of a pair must come from the same column and gradient. Supplying
  one and not the other is always a mistake, and is refused.

Path resolution has to work three ways: from a source checkout, from an
installed wheel (the data is package data, possibly inside a zip), and from a
PyInstaller bundle (which unpacks data to ``sys._MEIPASS`` and does not
participate in ``importlib.resources`` at all). :func:`reference_dir` tries all
three, in that order of trustworthiness.
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .errors import ConfigError, PanelError
from .io.loader import load_feature_table
from .io.schema import FeatureTable

#: Key of the reference pair used when the caller does not choose one.
DEFAULT_REFERENCE_KEY = "col35"

#: Registered bundled reference datasets: key -> static description.
#: ``label`` is what the UI shows; the files are resolved lazily so importing
#: this module never touches the filesystem.
REFERENCE_SETS: Dict[str, Dict] = {
    "col35": {
        "key": "col35",
        "label": "Column 25",
        "sub": "Human serum + Mix 4.4 standards, QTOF positive",
        "directory": "col35",
        "sample_file": "reference_sample.csv",
        "standards_file": "reference_standards.csv",
        "provenance_file": "provenance.json",
    },
}


@dataclass
class ReferencePair:
    """A resolved reference: where the two runs are and where they came from."""

    key: str                       # "col35", or "custom" for user-supplied files
    label: str
    sample_path: str
    standards_path: str
    is_default: bool
    provenance: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        """The ``reference`` block of ``model.json`` (spec §5.1)."""
        return {
            "key": self.key,
            "label": self.label,
            "sample": self.sample_path,
            "standards": self.standards_path,
            "is_default": bool(self.is_default),
            "provenance": dict(self.provenance),
        }


@dataclass
class LoadedReference:
    """A resolved reference with both runs parsed into feature tables."""

    pair: ReferencePair
    sample: FeatureTable
    standards: FeatureTable

    @property
    def provenance(self) -> Dict:
        return self.pair.to_dict()


def reference_root() -> str:
    """Directory holding the bundled ``reference_data`` tree.

    Order of attempts: PyInstaller's unpack directory, then package data via
    ``importlib.resources``, then the source tree next to this file. The frozen
    case comes first because in a bundle the other two can *appear* to succeed
    while pointing at a stale copy inside the zipped library.
    """
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        for cand in (os.path.join(meipass, "rt_anchor", "reference_data"),
                     os.path.join(meipass, "reference_data")):
            if os.path.isdir(cand):
                return cand

    try:
        from importlib import resources

        with resources.as_file(resources.files("rt_anchor") / "reference_data") as p:
            if os.path.isdir(str(p)):
                return str(p)
    except Exception:
        # zipped distribution without as_file support, or rt_anchor not
        # importable as a package resource — fall through to the path below
        pass

    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference_data")


def reference_dir(key: str = DEFAULT_REFERENCE_KEY) -> str:
    """Directory of one bundled reference set."""
    spec = _spec(key)
    return os.path.join(reference_root(), spec["directory"])


def available_references() -> List[Dict]:
    """Every bundled reference set, with its provenance, for ``rt-anchor references``."""
    out = []
    for key in REFERENCE_SETS:
        spec = dict(REFERENCE_SETS[key])
        spec["directory_path"] = reference_dir(key)
        spec["provenance"] = _provenance(key)
        spec["available"] = os.path.isfile(_path(key, "sample_file")) and \
            os.path.isfile(_path(key, "standards_file"))
        out.append(spec)
    return out


def resolve_reference(key: Optional[str] = None,
                      sample: Optional[str] = None,
                      standards: Optional[str] = None) -> ReferencePair:
    """Decide which reference pair to use.

    ``sample`` / ``standards`` override the bundle, and must be given together
    — a run calibrated against one lab's serum and another lab's standards is
    meaningless, and the resulting numbers would look perfectly plausible, so
    the half-specified case is refused rather than filled in.
    """
    if (sample is None) != (standards is None):
        raise ConfigError(
            "reference_sample and reference_standards must be given together "
            "(they must come from the same column and gradient). Got "
            f"reference_sample={sample!r}, reference_standards={standards!r}. "
            "Omit both to use the bundled Column 25 reference."
        )

    if sample is not None and standards is not None:
        for label, path in (("reference_sample", sample), ("reference_standards", standards)):
            if not os.path.isfile(path):
                raise ConfigError(f"{label} '{path}' does not exist.")
        return ReferencePair(
            key="custom",
            label="user-supplied reference pair",
            sample_path=os.path.abspath(sample),
            standards_path=os.path.abspath(standards),
            is_default=False,
            provenance={"note": "user-supplied; results are on this column's "
                                "time axis and are not comparable with "
                                "default-reference runs."},
        )

    k = key or DEFAULT_REFERENCE_KEY
    spec = _spec(k)
    sample_path = _path(k, "sample_file")
    standards_path = _path(k, "standards_file")
    missing = [p for p in (sample_path, standards_path) if not os.path.isfile(p)]
    if missing:
        raise PanelError(
            f"Bundled reference '{k}' is incomplete: {missing} not found under "
            f"{reference_dir(k)}. In a frozen build, check that "
            f"rt_anchor/reference_data/** is listed as data files in the "
            f".spec; from source, check the package was not installed with the "
            f"data excluded."
        )
    return ReferencePair(
        key=k,
        label=spec["label"],
        sample_path=sample_path,
        standards_path=standards_path,
        is_default=(k == DEFAULT_REFERENCE_KEY),
        provenance=_provenance(k),
    )


def load_reference(key: Optional[str] = None,
                   sample: Optional[str] = None,
                   standards: Optional[str] = None,
                   polarity: Optional[str] = None,
                   source_format: Optional[str] = None,
                   rt_unit: Optional[str] = None) -> LoadedReference:
    """Resolve a reference pair and load both runs as feature tables.

    ``source_format`` / ``rt_unit`` are forwarded to the loader; they apply to a
    user-supplied pair (the bundled pair is a plain MS-DIAL export in minutes
    and auto-detects correctly).
    """
    pair = resolve_reference(key=key, sample=sample, standards=standards)
    ref_sample = load_feature_table(pair.sample_path, source_format=source_format,
                                    polarity=polarity, rt_unit=rt_unit)
    ref_standards = load_feature_table(pair.standards_path, source_format=source_format,
                                       polarity=polarity, rt_unit=rt_unit)
    return LoadedReference(pair=pair, sample=ref_sample, standards=ref_standards)


# ---------------------------------------------------------------- internals ---

def _spec(key: str) -> Dict:
    if key not in REFERENCE_SETS:
        raise ConfigError(
            f"Unknown reference key '{key}'. Known: {sorted(REFERENCE_SETS)}. "
            f"To use your own runs pass reference_sample= and reference_standards=."
        )
    return REFERENCE_SETS[key]


def _path(key: str, which: str) -> str:
    return os.path.join(reference_dir(key), _spec(key)[which])


def _provenance(key: str) -> Dict:
    """Provenance JSON shipped alongside the data; ``{}`` when absent."""
    path = _path(key, "provenance_file")
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}
