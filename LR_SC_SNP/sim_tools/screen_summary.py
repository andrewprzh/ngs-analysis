"""Summaries of screen_sites.py output. Usage: screen_summary.py PREFIX"""
import sys
import numpy
import pandas as pd

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
prefix = sys.argv[1]
d = pd.read_csv(prefix + ".sites.tsv", sep="\t", low_memory=False)
s = d[d.screened == 1].copy()


def section(title):
    print("\n=== " + title)


# ---------------------------------------------------------------- recall funnel
section("Alt template funnel (screened TP + FN): fractions of alt templates")
t = s[s.group.isin(["TP", "FN"])].copy()
t["bin"] = pd.cut(t.alt_templates, [0, 9, 19, 49, 99, 1e9], labels=["5-9", "10-19", "20-49", "50-99", ">=100"])
g = t.groupby(["type", "bin"], observed=True)
f = pd.DataFrame({"sites": g.size(), "called": g.apply(lambda x: (x.group == "TP").mean()),
                  "alt_tpl": g.alt_templates.sum(),
                  "covered": g.tpl_covered.sum() / g.alt_templates.sum(),
                  "alt_given_cov": g.tpl_alt.sum() / g.tpl_covered.sum(),
                  "kept_given_cov": g.tpl_kept.sum() / g.tpl_covered.sum(),
                  "in_alt_mol": g.tpl_in_alt_mol.sum() / g.alt_templates.sum(),
                  "alt_umis/alt_tpl": g.est_alt_umis.sum() / g.alt_templates.sum()})
print(f.round(3).to_string())

section("Coverage of all templates: depth UMIs / covering templates (truth), TP + FN")
t["cov_ratio"] = t.est_depth_umis / t.covering_templates
t["tx3_bin"] = pd.cut(t.dist_tx3, [-1, 200, 500, 1000, 2000, 4000, 1e10])
print(t.groupby("tx3_bin", observed=True).agg(sites=("site", "size"), depth_per_cov_tpl=("cov_ratio", "median"),
                                               alt_cov=("tpl_covered", "sum"), alt_tpl=("alt_templates", "sum"))
      .assign(alt_covered=lambda x: x.alt_cov / x.alt_tpl).round(3).to_string())

section("Missed variants (FN screened): first failing cut-off, caller's view")
fn = s[s.group == "FN"].copy()


def reason(r):
    if r.est_alt_umis < 10:
        return "alt_umis<10"
    if r.est_alt_barcodes < 10:
        return "alt_barcodes<10"
    if r.est_depth_umis and r.est_alt_umis / r.est_depth_umis < 0.005:
        return "vaf<0.005"
    if r.est_pvalue > 1e-9:
        return "error_test"
    return "passes?"


fn["reason"] = fn.apply(reason, axis=1)
fn["bin"] = pd.cut(fn.alt_templates, [0, 9, 19, 49, 99, 1e9], labels=["5-9", "10-19", "20-49", "50-99", ">=100"])
print(pd.crosstab([fn.type, fn.bin], fn.reason, margins=True).to_string())
print("\nFN failing only the error test: est alt UMIs / depth / VAF / class")
e = fn[fn.reason == "error_test"]
if len(e):
    print(e.assign(vaf=e.est_alt_umis / e.est_depth_umis).groupby("type")
          .agg(n=("site", "size"), alt=("est_alt_umis", "median"), depth=("est_depth_umis", "median"),
               vaf=("vaf", "median")).to_string())
    print(e.error_class.value_counts().head(10).to_string())
print("\nFN passes? (should be few):")
print(fn[fn.reason == "passes?"][["site", "type", "alt_templates", "est_alt_umis", "est_alt_barcodes",
                                  "est_depth_umis", "est_pvalue", "dist_truth"]].head(15).to_string())

section("Lowering cut-offs: FN (screened, all alt template bins) that would pass")
for mu, mb, pv in [(10, 10, 1e-9), (5, 5, 1e-9), (5, 5, 1e-6), (3, 3, 1e-6), (10, 10, 1e-6), (8, 8, 1e-7)]:
    ok = (fn.est_alt_umis >= mu) & (fn.est_alt_barcodes >= mb) & (fn.est_pvalue <= pv) & \
         (fn.est_alt_umis >= 0.005 * fn.est_depth_umis)
    print("alt_umis>=%d barcodes>=%d p<=%g: %d of %d screened FN (%s)" %
          (mu, mb, pv, ok.sum(), len(fn), ", ".join("%s %d" % (k, v) for k, v in fn[ok].type.value_counts().items())))

# ---------------------------------------------------------------- FP vs TP features
section("FP vs TP: site features (all sites, not only screened)")
c = d[d.group.isin(["TP", "FP"])].copy()
feats = {
    "vaf>=0.5": c.vaf >= 0.5,
    "vaf>=0.9": c.vaf >= 0.9,
    "vaf<0.02": c.vaf < 0.02,
    "not exonic (expressed tx)": c.expr_tx_exonic == 0,
    "splice site <=3bp": c.dist_splice <= 3,
    "splice site <=10bp": c.dist_splice <= 10,
    "tx 3' end <=50bp": c.dist_tx3 <= 50,
    "tx 3' end <=200bp": c.dist_tx3 <= 200,
    "tx 5' end <=50bp": c.dist_tx5 <= 50,
    "several genes span": c.genes_spanning > 1,
    "both strands span": c.strands_spanning.astype(str).str.len() > 1,
    "true variant <=10bp": c.dist_truth <= 10,
    "homopolymer>=4": c.hp >= 4,
    "editing-like": c.editing_like == 1,
    "depth<30": c.depth_umis < 30,
}
rows = []
for name, m in feats.items():
    row = {"feature": name}
    for (grp, typ), x in c.groupby(["group", "type"]):
        row["%s_%s" % (grp, typ)] = round(float(m[x.index].mean()), 3)
    rows.append(row)
print(pd.DataFrame(rows).set_index("feature").to_string())
print("\ncounts:", c.groupby(["group", "type"]).size().to_dict())

section("FP vs TP: read features at screened sites (medians of per-site values)")
sc = s[s.group.isin(["TP", "FP"])]
cols = ["alt_mapq_lt60", "ref_mapq_lt60", "alt_mapq_lt20", "alt_near_end20", "ref_near_end20", "alt_qdist_med",
        "ref_qdist_med", "alt_near_intron10", "ref_near_intron10", "alt_src_exonic", "ref_src_exonic", "alt_src_gene",
        "ref_src_gene", "alt_fwd", "ref_fwd", "alt_read_err", "ref_read_err", "secondary", "supplementary"]
print(sc.groupby(["group", "type"])[cols].median().T.round(3).to_string())
print("\nmeans:")
print(sc.groupby(["group", "type"])[cols].mean().T.round(3).to_string())

section("FP SNV by VAF: read features (means)")
fs = s[(s.group == "FP")].copy()
fs["vbin"] = pd.cut(fs.vaf, [0, 0.02, 0.1, 0.5, 0.9, 1.01])
print(fs.groupby(["type", "vbin"], observed=True)[["alt_src_exonic", "alt_src_gene", "alt_mapq_lt60", "alt_near_end20",
                                                    "alt_near_intron10", "alt_fwd", "dist_splice", "dist_tx3"]]
      .agg(["mean"]).round(3).to_string())
print(fs.groupby(["type", "vbin"], observed=True).size().to_string())

section("tx_change of FP and TP SNVs")
sn = d[d.group.isin(["TP", "FP"]) & (d.type == "SNV")]
print(pd.crosstab(sn.tx_change, sn.group).to_string())
