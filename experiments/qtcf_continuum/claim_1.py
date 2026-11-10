"""Check Claim 1 and same-order interpolation on recorded Float64 nodes.

Read supports.csv from the continuum run and write the audit CSVs beside it.
"""

import argparse
import csv
from hashlib import sha256
from pathlib import Path

from sympy import QQ, QQ_I, ring
from sympy.polys.matrices import DomainMatrix


HERE = Path(__file__).resolve().parent
ODD_CASES = ("smooth_s", "two_semicircle_s", "odd_sine")
MIDPOINT = "QTCF"
POLYNOMIALS, T = ring("t", QQ_I)


# Preserve the exact binary value of each recorded Float64 component.
def gaussian(real, imag):
    return QQ_I.new(*(QQ(*float(value).as_integer_ratio()) for value in (real, imag)))


def read_supports(path):
    fits = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["method"] != MIDPOINT or int(row["index"]) > 11:
                continue
            node = gaussian(row["real"], row["imag"])
            value = gaussian(row["value_real"], row["value_imag"])
            fits.setdefault(row["case"], []).append((node, value))
    assert set(ODD_CASES) <= set(fits)
    return fits


def reduce_weight(data, weights, node, value):
    p, q = value, QQ_I.one
    for (selected, _), weight in zip(data, weights):
        p, q = (node - selected) * q, p - weight * q
    return p, q


def scalar_weights(data):
    weights = []
    for j, (node, value) in enumerate(data):
        p, q = reduce_weight(data[:j], weights, node, value)
        if not q:
            assert p
            return weights
        weights.append(p / q)
    return weights


def scalar_form(data, weights):
    p, q = POLYNOMIALS(weights[-1]), POLYNOMIALS.one
    for (node, _), weight in reversed(list(zip(data, weights))[:-1]):
        p, q = p * weight + (T - node) * q, p
    return p, q


def block_forms(data):
    previous = POLYNOMIALS.one, POLYNOMIALS.zero
    node, value = data[0]
    factor = T - node
    current = POLYNOMIALS(value), POLYNOMIALS.one
    yield current
    for j in range(1, len(data) - 1, 2):
        a, b = data[j:j + 2]
        reduced = []
        for node, value in (a, b):
            p, q = (polynomial.evaluate(T, node) for polynomial in current)
            pp, qq = (polynomial.evaluate(T, node) for polynomial in previous)
            error = p - value * q
            assert error, "An exact QTCF block was unavailable on the Float64 supports"
            reduced.append(-factor.evaluate(T, node) * (pp - value * qq) / error)
        line = (T - a[0]) * (reduced[1] - reduced[0]) / (b[0] - a[0]) + reduced[0]
        next_form = tuple(line * x + factor * y for x, y in zip(current, previous))
        previous, current = current, next_form
        factor = (T - a[0]) * (T - b[0])
        yield current


def coprime(form):
    p, q = form
    common = p.gcd(q)
    return p.exquo(common), q.exquo(common)


def homogeneous(data, p_degree, q_degree):
    rows = [[node**i for i in range(p_degree + 1)] +
            [-value * node**i for i in range(q_degree + 1)] for node, value in data]
    width = p_degree + q_degree + 2
    basis = DomainMatrix(rows, (len(rows), width), QQ_I).nullspace().to_list()
    assert len(basis) == 1
    coefficients = basis[0]
    return tuple(POLYNOMIALS.from_dict({(i,): c for i, c in enumerate(part)})
                 for part in (coefficients[:p_degree + 1], coefficients[p_degree + 1:]))


def same_quotient(left, right):
    p, q = left
    pp, qq = right
    return p * qq == pp * q


def attained(form, datum):
    p, q = form
    node, value = datum
    denominator = q.evaluate(T, node)
    return bool(denominator and p.evaluate(T, node) == value * denominator)


def is_odd(form):
    p, q = form
    return p * q.compose(T, -T) + p.compose(T, -T) * q == 0


def seeded_rows(cases):
    for case in ODD_CASES:
        data = cases[case]
        assert len(data) == 11 and data[0] == (QQ_I.zero, QQ_I.zero)
        assert all((-z, -y) == data[j + 1] for j, (z, y) in enumerate(data) if j % 2)
        blocks = list(block_forms(data))
        for k in range(5):
            old = data[:2 * k + 1]
            new, mate = data[2 * k + 1:2 * k + 3]
            old_form = coprime(homogeneous(old, k, k))
            assert is_odd(old_form) and same_quotient(old_form, blocks[k])
            next_form = coprime(homogeneous(old + [new], k + 1, k))
            assert is_odd(next_form)
            row = {
                "case": case,
                "complete_pairs": k,
                "prefix_numerator_degree": old_form[0].degree(),
                "prefix_denominator_degree": old_form[1].degree(),
                "homogeneous_new_node_attained": attained(next_form, new),
                "homogeneous_partner_attained": attained(next_form, mate),
                "homogeneous_next_quotient_odd": True,
            }
            if k % 2:
                # On these supports no finite scalar prefix reaches this degree pair.
                assert old_form[0].degree() == k and old_form[1].degree() == k - 1
                assert not attained(next_form, new) and not attained(next_form, mate)
                yield dict(row, seed_order="", first_weight="", partner_weight="unreachable")
                continue
            # Seed each prefix separately. The greedy paired run stops at k=0.
            order = [0]
            for j in range(1, 2 * k + 1, 4):
                order.extend((j, j + 2, j + 1, j + 3))
            seed = [old[j] for j in order]
            weights = scalar_weights(seed)
            assert len(weights) == len(seed)
            seeded = scalar_form(seed, weights)
            assert same_quotient(seeded, old_form) and all(attained(seeded, row) for row in seed)
            p, q = reduce_weight(seed, weights, *new)
            assert p and q
            first_weight = p / q
            first_form = scalar_form(seed + [new], weights + [first_weight])
            assert same_quotient(first_form, next_form)
            assert attained(first_form, new) and attained(first_form, mate)
            p, q = reduce_weight(seed + [new], weights + [first_weight], *mate)
            assert p and not q
            yield dict(row, seed_order=" ".join(str(i + 1) for i in order),
                       first_weight="finite", partner_weight="infinite")


def same_order_rows(cases):
    for case, data in cases.items():
        weights = scalar_weights(data)
        agreements = 0
        for k, quadratic in enumerate(block_forms(data)):
            count = 2 * k + 1
            if count > len(weights):
                break
            assert same_quotient(scalar_form(data[:count], weights[:count]), quadratic)
            agreements += 1
        yield {"case": case, "method": MIDPOINT, "startup": "midpoint",
               "proposed_supports": len(data),
               "finite_scalar_supports": len(weights),
               "first_infinite_scalar_insertion": len(weights) + 1 if len(weights) < len(data) else "",
               "equal_completed_prefixes": agreements}


def save_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("supports", nargs="?", type=Path, default=HERE / "supports.csv")
    args = parser.parse_args()
    cases = read_supports(args.supports)
    seeded = list(seeded_rows(cases))
    replay = list(same_order_rows(cases))
    save_csv(args.supports.parent / "seeded_claim.csv", seeded)
    save_csv(args.supports.parent / "same_order_replay.csv", replay)
    for row in seeded:
        print(row["case"], row["complete_pairs"], row["partner_weight"])
    count = sum(row["equal_completed_prefixes"] for row in replay)
    print(f"Verified {count} equal completed same-order prefixes for {MIDPOINT}.")
    print(f"Input supports.csv SHA-256 {sha256(args.supports.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
