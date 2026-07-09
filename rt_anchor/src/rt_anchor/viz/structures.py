"""Molecular-structure previews for the default standard panel (optional).

Renders the 15 built-in standards to small transparent SVGs via rdkit, for the
structure-on-hover feature in the interactive HTML report. rdkit is imported
lazily: if it is not installed, this returns ``{}`` and the report simply omits
the hover previews (everything else still works).
"""

from __future__ import annotations

from typing import Dict


def render_default_structures(size=(230, 165)) -> Dict[str, str]:
    """Return {standard_name: inline_svg} for the default panel, or {} if rdkit absent."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
    except Exception:
        return {}
    from ..panel import DEFAULT_SMILES

    out: Dict[str, str] = {}
    for name, smi in DEFAULT_SMILES.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        d = rdMolDraw2D.MolDraw2DSVG(*size)
        opt = d.drawOptions()
        opt.clearBackground = False          # transparent
        try:
            opt.bondLineWidth = 1
            opt.addStereoAnnotation = False
        except Exception:
            pass
        rdMolDraw2D.PrepareAndDrawMolecule(d, m)
        d.FinishDrawing()
        svg = d.GetDrawingText()
        i = svg.find("<svg")                 # strip the xml declaration for clean innerHTML
        out[name] = svg[i:] if i >= 0 else svg
    return out
