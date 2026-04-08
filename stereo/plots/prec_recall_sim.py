#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams['pdf.fonttype'] = 42

# --- Example data matching the plot ---
# Adjust these or load from a file
groups = ["Standard\nsimulation", "Concatenated\nreads"]
categories = ["Precision", "Recall"]
values = {
    groups[0]:     [99.95, 76.62],
    groups[1]: [99.94, 78.07],
}

# Pale versions inspired by image 2 but softer
pale_colors = ["#aad0f0", "#dde7a0"]  # pale blue, pale gold

x = np.arange(len(groups))
width = 0.3

fig, ax = plt.subplots(figsize=(4, 4))

for i, cat in enumerate(categories):
    vals = [values[g][i] for g in groups]
    ax.bar(x + i * width, vals, width, label=cat, color=pale_colors[i])

ax.set_ylabel("", fontsize=14)
ax.set_ylim(0, 110)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])

ax.set_xticks(x + width / 2)
ax.set_xticklabels(groups, fontsize=13)

#ax.set_title("Simulated reads", fontsize=15)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(frameon=False, fontsize=15, loc="upper center", bbox_to_anchor=(0.5, 1.2))

plt.tight_layout()
plt.savefig("simulated_reads.pdf", bbox_inches="tight")
print("Saved to simulated_reads.pdf")
