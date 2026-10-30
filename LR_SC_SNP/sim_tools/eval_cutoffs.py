"""Precision / recall of loose phase 1 calls at several error-test cut-offs, with and without unreliable flags,
and after the cell-type filter (run per cut-off, since BH depends on the tested set).

Usage: eval_cutoffs.py DATASET_DIR CELL_TYPES
Reads DATASET_DIR/lr_sc_snp/loose.{sites,eval.sites,eval.truth}.tsv, writes loose.p*.{sites,celltype_*}.tsv.
"""
import os
import subprocess
import sys

import pandas as pd

LR = "/home/andreyp/ngs-analysis/LR_SC_SNP"
CUTOFFS = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
d, cell_types = sys.argv[1], sys.argv[2]
o = os.path.join(d, "lr_sc_snp", "loose")

sites = pd.read_csv(o + ".sites.tsv", sep="\t", low_memory=False)
ev = pd.read_csv(o + ".eval.sites.tsv", sep="\t", low_memory=False).rename(columns={"#variant_id": "variant_id"})
truth = pd.read_csv(o + ".eval.truth.tsv", sep="\t", low_memory=False)
n_truth = truth.groupby("type").size()
ev["flagged"] = ev.unreliable != "."

rows = []
for cut in CUTOFFS:
    sel = ev[ev.error_pvalue <= cut]
    name = "p%g" % cut
    sub = sites[sites.error_pvalue <= cut]
    prefix = o + "." + name
    sub.to_csv(prefix + ".sites.tsv", sep="\t", index=False)
    subprocess.run(["/usr/bin/python3", LR + "/celltype_filter.py", "--sites", prefix + ".sites.tsv",
                    "--barcode_counts", o + ".barcode_counts.tsv.gz", "--cell_types", cell_types, "-o", prefix],
                   check=True, stdout=subprocess.DEVNULL)
    ct = pd.read_csv(prefix + ".celltype_sites.tsv", sep="\t").rename(columns={"#variant_id": "variant_id"})
    passed = set(ct.variant_id[ct["pass"] == 1])
    for vtype in ["SNV", "INS", "DEL"]:
        x = sel[sel.type == vtype]
        tp, fp = (x.label == "TP").sum(), (x.label == "FP").sum()
        u = x[~x.flagged]
        utp, ufp = (u.label == "TP").sum(), (u.label == "FP").sum()
        c = x[x.variant_id.isin(passed)]
        ctp, cfp = (c.label == "TP").sum(), (c.label == "FP").sum()
        nt = n_truth.get(vtype, 0)
        rows.append({"cutoff": name, "type": vtype, "TP": tp, "FP": fp,
                     "precision": tp / max(tp + fp, 1), "recall": x[x.label == "TP"].truth_id.nunique() / max(nt, 1),
                     "unflagged_TP": utp, "unflagged_FP": ufp, "unflagged_precision": utp / max(utp + ufp, 1),
                     "TP_flagged": 1 - utp / max(tp, 1), "FP_flagged": 1 - ufp / max(fp, 1),
                     "ct_pass_TP": ctp, "ct_pass_FP": cfp, "ct_precision": ctp / max(ctp + cfp, 1),
                     "ct_recall": c[c.label == "TP"].truth_id.nunique() / max(nt, 1)})
r = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(r.round(3).to_string(index=False))

print("\nflags by reason, p <= 1e-4 (share of calls)")
for reason in ["read_end", "editing", "homopolymer"]:
    ev[reason] = ev.unreliable.str.contains(reason)
print(ev.groupby(["label", "type"])[["read_end", "editing", "homopolymer", "flagged"]].mean().round(3).to_string())
