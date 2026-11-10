# Turn the measured errors, exact degrees, and interval certificates into LaTeX rows.
from csv import DictReader
from pathlib import Path
import math

HERE = Path(__file__).resolve().parent
L_CASES = {
    "ellipse_1_025": r"Ellipse \((1,1/4)\)",
    "ellipse_1_05": r"Ellipse \((1,1/2)\)",
    "upper_ellipse_1_05": "Upper ellipse",
    "five_lobed": "Five lobes",
    "superellipse_6": "Superellipse",
    "analytic_inlet": "Inlet",
}
AXES = {"real": "Real axis", "imaginary": "Imaginary axis"}
ODD = {"smooth_s": "Smooth S", "two_semicircle_s": "Two-semicircle S", "odd_sine": r"\(\sin x\)"}
LABELS = {
    "QTCF": "QTCF (midpoint)",
    "forced_Classic": "Paired TCF",
    "native_squared": "Native square TCF",
}
L_METHODS = ("QTCF", "forced_Classic")
P_METHODS = (*L_METHODS, "native_squared")
REL_TOL = 4e-11


def rows(path):
    with (HERE / path).open(newline="") as handle:
        return list(DictReader(handle))


def texnum(value, digits=3):
    value = float(value)
    if not math.isfinite(value):
        return r"\infty" if math.isinf(value) else r"\text{--}"
    label = f"{value:.{digits}g}"
    if "e" in label:
        base, exp = label.split("e")
        return rf"{base}\cdot10^{{{int(exp)}}}"
    return label


def tablerow(cells):
    return " & ".join(cells) + r"\\"


def main():
    summary = {(r["case"], r["method"]): r for r in rows("summary.csv")}
    degrees = {(r["case"], r["method"]): r["reduced_degree"] for r in
               rows("reduced_history_degrees.csv") if r["returned"] == "true"}
    certificates = {(r["case"], r["method"]): r for r in
                    rows("interval_certification/certificates.csv")}
    assert set(degrees) == {key for key, row in summary.items() if row["returned_block"] != "0"}

    def relative(case, method):
        r = summary[case, method]
        return float(r["error"]) / (float(r["validation_target"]) / REL_TOL)

    def fit_cells(case, method, include_poles):
        r = summary[case, method]
        if r["returned_block"] == "0":
            ncols = 5 if include_poles else 4
            label = "No complete pair" if method == "forced_Classic" and r["stop_reason"] == "nonfinite_weight" else "--"
            return [rf"\multicolumn{{{ncols}}}{{c}}{{{label}}}"]
        error = texnum(relative(case, method), 3)
        time = texnum(1000 * float(r["fit_seconds"]), 3)
        cells = [r["nodes"], degrees[case, method], "$" + error + "$", "$" + time + "$"]
        if include_poles:
            cells.append(poles(case, method))
        return cells

    def poles(case, method):
        value = certificates.get((case, method), {}).get("axis_poles", "-1")
        return "--" if value == "-1" else value

    def group_rows(cases, methods, include_poles=False):
        blocks = []
        for case, label in cases.items():
            blocks.append("\n".join(
                tablerow([label if i == 0 else "", LABELS[method], *fit_cells(case, method, include_poles)])
                for i, method in enumerate(methods)
            ))
        return ("\n" + r"\addlinespace[3pt]" + "\n").join(blocks)

    for name, table in (("group_l_rows", group_rows(L_CASES, L_METHODS)),
                        ("group_p_interval_rows", group_rows(AXES, P_METHODS, include_poles=True)),
                        ("group_p_odd_rows", group_rows(ODD, P_METHODS))):
        (HERE / f"{name}.tex").write_text(table + "\n")


if __name__ == "__main__":
    main()
