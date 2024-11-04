import sys


def fix_ref(reference_name):
    if reference_name == "MT":
        return 'chrM'
    elif '.' not in reference_name:
        return 'chr' + reference_name
    else:
        return reference_name


for l in sys.stdin:
    if l.startswith("@"):
        if l.startswith("@SQ"):
            v = l.strip().split('\t')
            new_ref = fix_ref(v[1].split(':')[1])
            sys.stdout.write("@SQ\tSN:%s\t%s\n" % (new_ref, '\t'.join(v[2:])))
        else:
            sys.stdout.write(l)
        continue

    v = l.strip().split('\t')
    new_ref = fix_ref(v[2])
    sys.stdout.write("%s\t%s\t%s\n" % ('\t'.join(v[:2]), new_ref, '\t'.join(v[3:])))
