#!/usr/bin/env python3
"""Regression checks for the LRA volume pipeline.

The script runs the ``ttc`` binary on a given SMT-LIB input and validates that
the reported polytope statistics match the expectations supplied on the command
line. It is intentionally lightweight so that it can be reused for multiple
scenarios via CTest.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple


Row = dict[str, float | int]


def parse_output(output: str) -> Tuple[Optional[int], List[Row], Optional[int]]:
    polytopes = None
    rows: List[Row] = []
    final_volume = None
    pattern = re.compile(r"^c\s+(\d+)\s+([0-9.]+)\s+(\d+)\s+(\d+)")
    for line in output.splitlines():
        line = line.strip()
        if line.startswith("c polytopes:"):
            try:
                polytopes = int(line.split(":", 1)[1])
            except ValueError as exc:  # pragma: no cover - defensive guard
                raise AssertionError(
                    f"Unable to parse polytopes count from '{line}'"
                ) from exc
        else:
            match = pattern.match(line)
            if match:
                idx = int(match.group(1))
                volume = float(match.group(2))
                deleted = int(match.group(3))
                total = int(match.group(4))
                rows.append({
                    "index": idx,
                    "volume": volume,
                    "deleted": deleted,
                    "total": total,
                })
        if line.startswith("s volume"):
            try:
                final_volume = int(line.split()[2])
            except (IndexError, ValueError) as exc:
                raise AssertionError(
                    f"Unable to parse final volume from '{line}'"
                ) from exc
    return polytopes, rows, final_volume


def _normalize_expected_sequence(
    expected: Optional[Sequence[float | int]],
    count: int,
) -> Optional[List[float | int]]:
    if expected is None or len(expected) == 0:
        return None
    if len(expected) == 1:
        return [expected[0]] * count
    if len(expected) != count:
        raise AssertionError(
            f"Expected {count} entries but received {len(expected)}"
        )
    return list(expected)


def normalize_expected_volumes(
    expected: Optional[Sequence[float]],
    count: int,
) -> Optional[List[float]]:
    normalized = _normalize_expected_sequence(expected, count)
    if normalized is None:
        return None
    return [float(value) for value in normalized]


def normalize_expected_integers(
    expected: Optional[Sequence[int]],
    count: int,
) -> Optional[List[int]]:
    normalized = _normalize_expected_sequence(expected, count)
    if normalized is None:
        return None
    return [int(value) for value in normalized]


def check_volumes(
    rows: Iterable[Row],
    expected_volumes: Sequence[float],
    tolerance: float,
) -> None:
    for row, expected in zip(rows, expected_volumes):
        lower = expected - tolerance
        upper = expected + tolerance
        actual = float(row["volume"])
        if not (lower <= actual <= upper):
            raise AssertionError(
                f"Polytope {row['index']} volume {actual:.3f} outside expected range"
            )


def check_integer_field(
    rows: Iterable[Row],
    field: str,
    expected_values: Sequence[int],
) -> None:
    for row, expected in zip(rows, expected_values):
        actual = int(row[field])
        if actual != expected:
            raise AssertionError(
                f"Polytope {row['index']} {field} {actual} does not match expected {expected}"
            )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path, help="Path to the ttc executable")
    parser.add_argument("input_file", type=Path, help="SMT-LIB input file")
    parser.add_argument(
        "--expected-polytopes",
        type=int,
        required=True,
        help="Expected number of polytopes in the output",
    )
    parser.add_argument(
        "--expected-rows",
        type=int,
        required=True,
        help="Expected number of rows in the volume table",
    )
    parser.add_argument(
        "--expected-volumes",
        type=float,
        nargs="*",
        help="Expected volume for each row (single value applies to all rows)",
    )
    parser.add_argument(
        "--expected-deleted",
        type=int,
        nargs="*",
        help="Expected number of deleted samples per row",
    )
    parser.add_argument(
        "--volume-tolerance",
        type=float,
        default=80.0,
        help="Tolerance applied when comparing expected and actual volumes",
    )
    parser.add_argument(
        "--final-range",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        help="Inclusive range for the final volume",
    )
    parser.add_argument(
        "--require-final-equals-total",
        action="store_true",
        help="Assert that the final total samples matches the reported volume",
    )
    parser.add_argument(
        "--skip-final-total-check",
        dest="require_final_equals_total",
        action="store_false",
        help="Disable the final total/volume consistency check",
    )
    parser.set_defaults(require_final_equals_total=True)

    args = parser.parse_args(argv)

    if not args.binary.exists():
        raise AssertionError(f"Binary '{args.binary}' does not exist")
    if not args.input_file.exists():
        raise AssertionError(f"Input file '{args.input_file}' does not exist")

    completed = subprocess.run(
        [str(args.binary), str(args.input_file)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    polytopes, rows, final_volume = parse_output(completed.stdout)

    if polytopes != args.expected_polytopes:
        raise AssertionError(
            f"Expected {args.expected_polytopes} polytopes, found {polytopes}"
        )
    if len(rows) != args.expected_rows:
        raise AssertionError(
            f"Expected {args.expected_rows} table rows, found {len(rows)}"
        )

    normalized = normalize_expected_volumes(args.expected_volumes, len(rows))
    if normalized is not None:
        check_volumes(rows, normalized, args.volume_tolerance)

    expected_deleted = normalize_expected_integers(
        args.expected_deleted, len(rows)
    )
    if expected_deleted is not None:
        check_integer_field(rows, "deleted", expected_deleted)

    if final_volume is None:
        raise AssertionError("Missing final volume line in output")

    if args.final_range is not None:
        lower, upper = args.final_range
        if not (lower <= final_volume <= upper):
            raise AssertionError(
                f"Final volume {final_volume} outside expected range [{lower}, {upper}]"
            )

    if args.require_final_equals_total:
        final_total = int(rows[-1]["total"]) if rows else 0
        if final_total != final_volume:
            raise AssertionError(
                f"Final total samples {final_total} does not match reported volume {final_volume}"
            )

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as exc:  # pragma: no cover - used for CTest diagnostics
        print(str(exc), file=sys.stderr)
        sys.exit(1)
