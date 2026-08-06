#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

"""Plot the flip-probability matrix produced by `avalanche --matrix-output` as a heatmap.

Requires: pandas, numpy, matplotlib (pip install pandas numpy matplotlib)

Usage:
    python3 plot_avalanche.py matrix.hash0.csv
    python3 plot_avalanche.py matrix.hash0.csv -o out.png --annotate
"""

import argparse
import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import TwoSlopeNorm


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", help="CSV file written by `avalanche --matrix-output`")
    parser.add_argument("-o", "--output", help="Output image path (default: <csv-stem>.png)")
    parser.add_argument(
        "--all-columns",
        action="store_true",
        help="Show all 64 output-bit columns, including ones that are structurally always 0 "
        "(by default only 'live' columns, i.e. those with at least one nonzero entry, are shown)",
    )
    parser.add_argument(
        "--clip",
        type=float,
        default=0.5,
        help="Clamp the color scale to 0.5 +/- CLIP (default: 0.5). "
        "Values beyond this are shown fully saturated.",
    )
    parser.add_argument("--annotate", action="store_true", help="Write the probability value into each cell")
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args()

    df = pd.read_csv(args.csv, index_col=0)
    df = df.iloc[::-1]
    df.columns = df.columns.astype(int)
    df.index = df.index.astype(int)

    if not args.all_columns:
        live_cols = df.columns[(df != 0).any(axis=0)]
        if len(live_cols) == 0:
            sys.exit("No non-zero columns found in the matrix; pass --all-columns to see the raw grid anyway.")
        dropped = df.shape[1] - len(live_cols)
        df = df[live_cols]
        if dropped:
            print(f"Hiding {dropped} column(s) that are structurally always 0 (pass --all-columns to show them).")
    df = df.transpose()
    data = df.to_numpy()
    vmin, vmax = 0.5 - args.clip, 0.5 + args.clip

    fig_h = max(4.0, 0.2 * len(df.index))
    fig_w = max(6.0, 0.2 * len(df.columns))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.5, vmax=vmax)
    im = ax.imshow(data, cmap="RdBu_r", norm=norm, aspect="auto")

    x_labs = list(df.columns)
    x_labs[1::2] = ("" for _ in x_labs[1::2])
    ax.set_xticks(range(len(df.columns)), x_labs)
    # ax.set_xticks(range(len(df.columns)))
    # ax.set_xticklabels(df.columns[::2])
    y_labs = list(str(x) for x in df.index)
    y_labs[1::2] = ("" for _ in y_labs[1::2])
    ax.set_yticks(range(len(df.index)), y_labs)
    ax.set_ylabel("Output bit")
    ax.set_xlabel("Input bit")
    ax.set_title("P(output bit flips | input bit flipped)")

    if args.annotate:
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                color = "white" if abs(data[i, j] - 0.5) > args.clip * 0.6 else "black"
                ax.text(j, i, f"{data[i, j]:.2f}", ha="center", va="center", fontsize=5, color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    if args.clip == 0.5:
        cbar.set_label(f"P(flip)")
    else:
        cbar.set_label(f"P(flip)  (clamped to [{vmin:.2f}, {vmax:.2f}])")

    fig.tight_layout()
    out = args.output or (args.csv.rsplit(".", 1)[0] + ".png")
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
