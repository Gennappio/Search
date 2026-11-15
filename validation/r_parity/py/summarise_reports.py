"""One table per report: for every case and R variant, how many compared quantities are equal,
how many differ, whether R stopped with an error, and the largest difference among the equal
ones. Written to ``summary.md`` next to the reports. Usage::

    python validation/r_parity/py/summarise_reports.py <results directory>
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


def summarise(report: pd.DataFrame) -> pd.DataFrame:
    rows = []
    checked = report[report["status"] != "info"]
    for (case, variant), group in checked.groupby(["case", "variant"], sort=False):
        equal = group[group["status"] == "equal"]
        rows.append(
            {
                "case": case,
                "R run": variant,
                "quantities": len(group),
                "equal": len(equal),
                "differ": int((group["status"] == "differs").sum()),
                "R error": int((group["status"] == "r_error").sum()),
                "values compared": int(equal["n"].sum()),
                "largest difference among the equal": (
                    f"{equal['max_abs_diff'].max():.1e}" if len(equal) else "-"
                ),
            }
        )
    return pd.DataFrame(rows)


def to_markdown(table: pd.DataFrame) -> str:
    header = "| " + " | ".join(table.columns) + " |"
    rule = "|" + "|".join("---" for _ in table.columns) + "|"
    body = ["| " + " | ".join(str(v) for v in row) + " |" for row in table.itertuples(index=False)]
    return "\n".join([header, rule, *body])


def main() -> None:
    results = Path(sys.argv[1])
    parts = ["# R parity reports, summary", ""]
    for path in sorted(results.glob("report_*.csv")):
        report = pd.read_csv(path)
        table = summarise(report)
        total = table[["quantities", "equal", "differ", "R error"]].sum()
        parts += [
            f"## {path.name}",
            "",
            f"{int(total['quantities'])} compared quantities: {int(total['equal'])} equal, "
            f"{int(total['differ'])} different, {int(total['R error'])} where R stopped with an "
            "error.",
            "",
            to_markdown(table),
            "",
        ]
        info = report[report["status"] == "info"]
        if len(info):
            parts += ["Information rows:", ""]
            parts += [
                f"- {r.case}, {r.variant}, {r.quantity}: {r.detail}" for r in info.itertuples()
            ]
            parts.append("")
    text = "\n".join(parts)
    (results / "summary.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
