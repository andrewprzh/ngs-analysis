#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

# --- read TSV: label\tvalue per line ---
tsv_file = sys.argv[1] if len(sys.argv) > 1 else "data.tsv"
labels, values = [], []
with open(tsv_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        labels.append(parts[0])
        values.append(int(parts[1]))

values_m = [v / 1_000_000 for v in values]

def_colors = ["#b0b0b0", "#c4a951", "#c9a0c9", "#a8b4be", "#7baed4", "#e8d76a", "#a8c98a"]

colors = list(reversed(list(reversed(def_colors))[:len(values_m)]))

fig, ax = plt.subplots(figsize=(6, 5))
bars = ax.bar(range(len(labels)), values_m, color=colors[:len(labels)], width=0.7)

ax.set_ylabel("Read Number (Million)", fontsize=17)
ax.set_xticks([])
ax.set_ylim(0, max(values_m))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# legend matching bar colors and multi-word labels
legend_labels = [l.replace(" ", "\n") for l in labels]
ax.legend(bars, legend_labels, loc="center left", bbox_to_anchor=(1.02, 0.64),
          frameon=False, fontsize=17, handlelength=1.5)

plt.tight_layout()
out = tsv_file.rsplit(".", 1)[0] + ".pdf"
plt.savefig(out, bbox_inches="tight")
print(f"Saved to {out}")
