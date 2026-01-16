#!/usr/bin/env python3
"""
Extract Crosslink Residue Pairs

This script extracts residue pairs and spectral counts from raw XL-MS output
files that contain formatted strings like "XP_037295533.1 (1500)-XP_037295533.1 (1509)".

Converts these to a simple 3-column format: Residue1, Residue2, Spectral Count

Usage:
    python extract_residues.py --input raw_crosslinks.csv --output crosslinks.csv
"""

import argparse
import csv
import re
from pathlib import Path


def extract_residues(s: str):
    """
    Extract residue numbers from a formatted crosslink string.

    Args:
        s: String like 'XP_037295533.1 (1500)-XP_037295533.1 (1509)'

    Returns:
        Tuple of (residue1, residue2) as strings, or (None, None) if not found.
    """
    nums = re.findall(r'\((\d+)\)', s)
    if len(nums) >= 2:
        return nums[0], nums[1]
    return None, None


def process_file(input_path: str, output_path: str, delimiter: str = ","):
    """
    Process a raw crosslink file and extract residue pairs.

    Args:
        input_path: Path to input CSV file.
        output_path: Path for output CSV file.
        delimiter: CSV delimiter character.
    """
    rows_to_write = []

    with open(input_path, newline='', encoding="utf-8") as f_in:
        reader = csv.reader(f_in, delimiter=delimiter)

        for row in reader:
            if len(row) < 2:
                continue

            col1 = row[0].strip()
            col2 = row[1].strip()

            res1, res2 = extract_residues(col1)
            if res1 is None or res2 is None:
                continue

            if not col2:
                continue

            rows_to_write.append((res1, res2, col2))

    with open(output_path, "w", newline='', encoding="utf-8") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["Residue1", "Residue2", "Spectral Count"])
        writer.writerows(rows_to_write)

    print(f"Extracted {len(rows_to_write)} crosslinks")
    print(f"Output saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract residue pairs from raw XL-MS output"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input CSV file with raw crosslink data"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output CSV file for extracted residue pairs"
    )
    parser.add_argument(
        "--delimiter", "-d",
        default=",",
        help="CSV delimiter (default: comma)"
    )

    args = parser.parse_args()
    process_file(args.input, args.output, args.delimiter)


if __name__ == "__main__":
    main()
