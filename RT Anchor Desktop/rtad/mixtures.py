"""Bundled standard mixtures — a thin re-export of :mod:`rt_anchor.mixtures`.

The panel definitions used to live here, which meant the desktop app and the
calibration engine each held their own idea of what "Mix 21" is. Under the v2
cross-column method the panel choice moves the iRT ruler (spec §6, §10.2), so
two definitions that quietly disagree would produce two different scales with no
error anywhere. They now live in the engine — ``rt_anchor/mixtures.py`` — and
this module exists only so that ``from rtad.mixtures import ...`` keeps working
for the app's own code and for ``v2/5 chroms/iRT/run_iRT.py``, which imports
:data:`MIX21LPC` from here.

Nothing is redefined below. If you are looking for the manifests, they are in
the engine; edit them there.
"""

from __future__ import annotations

from rt_anchor.mixtures import (  # noqa: F401  (re-export)
    DEFAULT_MIXTURE,
    MIX15,
    MIX21,
    MIX21LPC,
    NONE,
    OFFERED_MIXTURES,
    get_manifest,
    get_meta,
    mixture_keys,
    n_standards,
    offered_keys,
    preview_payload,
)

__all__ = [
    "DEFAULT_MIXTURE", "MIX15", "MIX21", "MIX21LPC", "NONE", "OFFERED_MIXTURES",
    "get_manifest", "get_meta", "mixture_keys", "n_standards", "offered_keys",
    "preview_payload",
]
