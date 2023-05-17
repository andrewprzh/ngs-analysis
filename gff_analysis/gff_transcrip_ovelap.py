import sys


def load_transcript_ids(inf):
    transcript_ids = set()
    for l in open(inf):
        if l.startswith("#"): continue
        v = l.strip().split("\t")
        if len(v) < 9: continue
        if v[2] != "transcript": continue

        atts = v[8].split(" ")
        tid = None
        for i in range(len(atts) - 1):
            if atts[i] == "transcript_id":
                tid = atts[i+1][1:-2]
                break
        if tid: transcript_ids.add(tid)

    return transcript_ids


transcripts1 = load_transcript_ids(sys.argv[1])
transcripts2 = load_transcript_ids(sys.argv[2])

print("GTF1: %d" % len(transcripts1))
print("GTF2: %d" % len(transcripts2))
print("Overlap: %d" % len(transcripts1.intersection(transcripts2)))
print(transcripts2.difference(transcripts1))
