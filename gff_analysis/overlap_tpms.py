import sys


def load_expressed_known_transcripts(infile):
    transcripts = set()
    for l in open(infile):
        if l.startswith("#"): continue
        if l.find("nic") != -1: continue

        v = l.strip().split("\t")
        tid = v[0]
        if float(v[1]) > 0:
            transcripts.add(tid)

    return transcripts


transcripts1 = load_expressed_known_transcripts(sys.argv[1])
transcripts2 = load_expressed_known_transcripts(sys.argv[2])

print("First file: %d" % len(transcripts1))
print("Second file: %d" % len(transcripts2))
print("Overlap: %d" % len(transcripts1.intersection(transcripts2)))

diff12 = transcripts1.difference(transcripts2)
print("First only: %d" % len(diff12))
for t in sorted(diff12):
    print("  %s" % t)

diff21 = transcripts2.difference(transcripts1)
print("Second only: %d" % len(diff21))
for t in sorted(diff21):
    print("  %s" % t)