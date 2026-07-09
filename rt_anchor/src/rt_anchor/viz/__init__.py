"""Visualisation: blue-themed report in interactive HTML + static PDF.

Three modules — metrics, TIC (before/after), performance — assembled into one
report. matplotlib renders the PDF, Plotly the interactive HTML. Import is lazy
so the core package does not require plotting libraries unless a report is made.
"""

from .report import build_html, build_pdf, write_report

__all__ = ["write_report", "build_html", "build_pdf"]
