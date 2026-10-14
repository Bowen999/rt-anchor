# Negative ESI, 20-min method

Serum and standards-mixture runs of the 20-min LC method on the same QTOF as
[`../6_columns/`](../6_columns/), but acquired in **negative mode**, exported
from MS-DIAL (area tables). Use `polarity="negative"`.

| file | content | features |
|---|---|---|
| `samples.csv` | serum | 399 |
| `standards.csv` | standards mixture, default threshold | 55 |
| `standards_low_threshold.csv` | same standards run re-exported with a lower detection threshold | 113 |

Lipid names in the tables are m/z-only MS-DIAL guesses; none is MS/MS-confirmed.
The manually vetted lipid list for these data is
[`../lipid_vetted_list.xlsx`](../lipid_vetted_list.xlsx).
