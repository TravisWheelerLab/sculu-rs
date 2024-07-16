#!/usr/bin/env python3
"""
Author : Ken Youens-Clark <kyclark@arizona.edu>
Date   : 2024-06-28
Purpose: Find best score
"""

import argparse
import csv
from pprint import pprint
from collections import defaultdict
from itertools import repeat, permutations
from typing import List, NamedTuple, TextIO


class Args(NamedTuple):
    """Command-line arguments"""

    file: TextIO
    lambda_value: float


# --------------------------------------------------
def get_args() -> Args:
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Find best score",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "file",
        help="Input file",
        metavar="FILE",
        type=argparse.FileType("rt"),
    )

    # TODO: Should be loaded from matrix?
    parser.add_argument(
        "-l",
        "--lamb",
        help="Lambda valuw",
        metavar="LAMBDA",
        type=float,
        default=0.1227
    )

    args = parser.parse_args()

    return Args(args.file, args.lamb)


# --------------------------------------------------
def main() -> None:
    """Make a jazz noise here"""

    args = get_args()

    reader = csv.DictReader(
        args.file, delimiter="\t", fieldnames=["score", "target", "query"]
    )

    scores = {}
    for rec in reader:
        target = rec["target"]
        query = rec["query"]
        if rec["score"].isdigit():
            score = int(rec["score"])

            if target not in scores:
                scores[target] = defaultdict(int)
            scores[target][query] = max([score, scores[target][query]])

    print("TARGET/QUERY")
    pprint(scores)
    clear_winners = defaultdict(int)
    winning_sets = defaultdict(int)

    for target, bit_scores in scores.items():
        print(">>>", target)
        families = list(bit_scores.keys())
        print(families)
        print(list(bit_scores.values()))
        conf = mk_conf(bit_scores.values(), args.lambda_value)
        print(list(map(lambda v: f"{v:0.04f}", conf)))
        pos = list(range(len(conf)))
        all_comps = []
        for i in pos:
            val = conf[i]
            others = [conf[j] for j in pos if j != i]
            pairs = zip(repeat(val, len(others)), others)
            comparisons = [(x * 0.3) > y for x, y in pairs]
            print(
                "{:8} {:0.04f} {:5} => {}".format(
                    families[i], val, str(all(comparisons)), comparisons
                )
            )
            all_comps.append(all(comparisons))

        if winner := [fam for fam, win in zip(families, all_comps) if win]:
            winner = winner[0]
            clear_winners[winner] += 1
            print("Clear Winner:", winner)
        else:
            print("NO WINNER")
            fam_comps = sorted(zip(conf, families))
            top_conf = fam_comps[-1][0]
            thresh = top_conf / 5
            print(f"Top Conf: {top_conf}")
            print(f"Thresh: {thresh}")
            winning_set = [f for c, f in fam_comps if c > thresh]
            print('Winning Set:', winning_set)
            for pair in permutations(winning_set, 2):
                winning_sets[tuple(sorted(pair))] += 1

        print()

    print('Clean Winners')
    pprint(clear_winners)
    print()

    print('Winning Sets')
    pprint(winning_sets)

    independence = []
    for f1, f2 in permutations(families, 2):
        ind = 1
        num_shared = winning_sets[tuple(sorted([f1, f2]))]

        if num_shared > 0:
            num_wins = clear_winners[f1]
            ind = num_wins / (num_wins + num_shared)

        print(f"ind {ind} f1 {f1} f2 {f2} "
              f"num_shared {num_shared} num_wins {num_wins}")

        # num_shared is used to break ties
        independence.append((ind, num_shared, f1, f2))

    for ind, _, f1, f2 in sorted(independence):
        print(f'Independence {f1:6}/{f2:6}: {ind:0.02f}')

    print("Done")


# --------------------------------------------------
def mk_conf(
    vals: List[int],
    lamb: float,
) -> List[float]:
    """
    same function from PolyA to calculate confidence
    """
    converted = list(map(lambda v: 2 ** (v * lamb), vals))
    total = sum(converted)
    return list(map(lambda v: v / total, converted))


# --------------------------------------------------
if __name__ == "__main__":
    main()
