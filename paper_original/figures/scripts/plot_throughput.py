import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 15

CATS = ["TC", "WG", "WT", "CP", "LJ", "TW", "FS", "Avg"]

INSERTION = {
    "Terrace": [1.7, 1.5, 4.0, 3.9, 2.0, 2.8, 2.7, 2.7],
    "Aspen":   [1.3, 1.2, 6.5, 2.8, 6.3, 3.8, 2.0, 3.4],
    "DZiG":    [9.4, 9.6, 9.6, 17.2, 11.2, 6.7, 5.6, 9.9],
    "Spruce":  [8.7, 8.7, 9.0, 16.2, 10.3, 6.0, 5.1, 9.1],
}

DELETION = {
    "Terrace": [2.9, 2.1, 4.6, 3.3, 4.7, 4.0, 2.5, 3.4],
    "Aspen":   [2.5, 1.7, 7.8, 3.8, 9.1, 3.5, 2.2, 4.4],
    "DZiG":    [8.0, 9.1, 10.2, 14.7, 11.9, 6.5, 4.6, 9.3],
    "Spruce":  [7.3, 8.4, 9.5, 13.5, 11.0, 5.8, 4.0, 8.5],
}

SERIES_ORDER = ["Terrace", "Aspen", "DZiG", "Spruce"]

COLORS = {
    "Terrace": "#B5D5A4",
    "Aspen":   "#A8C8E8",
    "DZiG":    "#C9B3D5",
    "Spruce":  "#F4B8C1",
}

def plot(data, out_path, show_legend):
    fig, ax = plt.subplots(figsize=(10.72, 2.72))

    n_groups = len(CATS)
    n_series = len(SERIES_ORDER)
    bar_w = 0.23
    x = np.arange(n_groups, dtype=float)  # uniform spacing, including before "Avg"

    mid = (n_series - 1) / 2.0
    for i, series in enumerate(SERIES_ORDER):
        vals = data[series]
        offset = (i - mid) * bar_w
        bars = ax.bar(x + offset, vals, width=bar_w, color=COLORS[series],
                       edgecolor="black", linewidth=0.9, label=series, zorder=3)
        for xv, v in zip(x + offset, vals):
            ax.text(xv, v + 0.3, f"{v:.1f}", ha="center", va="bottom", fontsize=10.3)

    ax.axhline(1, color="black", linestyle="--", linewidth=0.9, zorder=5)

    ax.set_xticks(x)
    ax.set_xticklabels(CATS, fontsize=15)
    ax.set_ylabel("Speedup", fontsize=15)
    ax.set_ylim(0, 20.5)
    ax.set_yticks(range(0, 19, 2))
    ax.tick_params(axis="y", labelsize=13)
    ax.set_xlim(x[0] - 0.55, x[-1] + 0.55)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_linewidth(1.0)

    if show_legend:
        handles, labels = ax.get_legend_handles_labels()
        order = [3, 2, 1, 0]  # Spruce, DZiG, Aspen, Terrace
        ax.legend([handles[i] for i in order], [labels[i] for i in order],
                  loc="upper right", frameon=True, fontsize=12,
                  handlelength=1.3, handleheight=1.0, labelspacing=0.2,
                  borderpad=0.35, bbox_to_anchor=(1.01, 1.1))

    fig.tight_layout(pad=0.3)
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)

plot(INSERTION, os.path.join(OUT_DIR, "insertion_throughput_regen.pdf"), True)
plot(DELETION, os.path.join(OUT_DIR, "deletion_throughput_regen.pdf"), False)
print("done")
