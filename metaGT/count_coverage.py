import sys
import pysam

cov90 = set()
match90 = set()
unaligned = set()
all_contigs = set()

samf = pysam.AlignmentFile(sys.argv[1], "r")
for a in samf:
    all_contigs.add(a.query_name)
    if a.is_unmapped:
        unaligned.add(a.query_name)
        continue

    ctg_len = a.query_length
    ctg_aligned = a.query_alignment_length
    ref_len = samf.get_reference_length(a.reference_name)
    ref_aligned = a.reference_length

    if ref_aligned >=  ref_len * 0.975:
        cov90.add(a.reference_name)
        if ctg_aligned >= ctg_len * 0.975:
            match90.add(a.reference_name)

print("Total: ", len(all_contigs))
print("Unmapped: ", len(unaligned))
print("Genes covered: ", len(cov90))
print("Genes covered no extra", len(match90))
