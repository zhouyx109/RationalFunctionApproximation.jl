# Replay the original absolute-value reductions exactly on the returned nodes.
from csv import DictReader, DictWriter
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_rows(name):
    with (HERE / name).open(newline="") as stream:
        return list(DictReader(stream))


def main():
    # supports.csv contains the returned fit. The trace can include later blocks.
    supports = read_rows("supports.csv")
    reductions = read_rows("reductions.csv")
    rows = []
    for axis, coordinate in (("real", "real"), ("imaginary", "imag")):
        points = [
            Fraction(float(row[coordinate]))
            for row in supports
            if row["case"] == axis and row["method"] == "QTCF"
        ]
        assert points[0] == 0 and len(points) % 2 == 1
        peaks = {}
        for row in reductions:
            if row["case"] == axis and row["method"] == "QTCF":
                block = int(row["block"])
                peaks[block] = max(peaks.get(block, 0), float(row["reduced_max"]))

        # In coordinates along either axis, the positive-half-axis
        # magnitude is d*x+c. The sign on the negative half-axis reverses.
        d, c = Fraction(0), Fraction(1)
        for block, i in enumerate(range(1, len(points), 2), start=2):
            assert points[i] == -points[i+1] and points[i] != 0
            radius = abs(points[i])
            rows.append(dict(
                axis=axis, block=block, radius=float(radius),
                exact_d_before_rounded=float(d), exact_c_before_rounded=float(c),
                exact_reduced_abs_rounded=float(abs(d*radius+c)),
                classic_float64_reduced_abs=peaks[block],
            ))
            d, c = -radius/c, -radius*radius/c

    with (HERE / "interval_reduction_replay.csv").open("w", newline="") as stream:
        writer = DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
