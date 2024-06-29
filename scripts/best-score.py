#!/usr/bin/env python3
"""
Author : Ken Youens-Clark <kyclark@arizona.edu>
Date   : 2024-06-28
Purpose: Find best score
"""

import argparse
import csv
import re
from pprint import pprint
from collections import defaultdict
from typing import NamedTuple, TextIO


class Args(NamedTuple):
    """Command-line arguments"""

    file: TextIO


# --------------------------------------------------
def get_args() -> Args:
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Find best score",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "file", help="Input file", metavar="FILE", type=argparse.FileType("rt")
    )

    args = parser.parse_args()

    return Args(args.file)


# --------------------------------------------------
def main() -> None:
    """Make a jazz noise here"""

    args = get_args()
    scores = {}

    reader = csv.DictReader(
        args.file, delimiter="\t", fieldnames=["score", "target", "query"]
    )

    for rec in reader:
        target = re.sub('__.*', '', rec['target'])
        query = rec['query']
        if rec['score'].isdigit():
            if target not in scores:
                scores[target] = defaultdict(int)
            scores[target][query] = max([int(rec['score']), scores[target][query]])

    pprint(scores)

    for target, tscores in scores.items():
        (score, query) = sorted((v, k) for (k, v) in tscores.items())[-1]
        print(f'target {target:10} query {query:10} score {score:>}')
    print("Done")


# --------------------------------------------------
if __name__ == "__main__":
    main()
