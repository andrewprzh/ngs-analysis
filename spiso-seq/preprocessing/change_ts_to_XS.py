import os
import sys

f = open(sys.argv[1])

for l in f:
    if (l.find("ts:A:") != -1):
        tokens = l.strip().split('\t')
        ts_flag_pos = 0
        for i in range(0,len(tokens)):
            if tokens[i].startswith("ts:A:"):
                ts_flag_pos = i
                break

        if len(tokens) < 6 or ts_flag_pos == 0:
            print(l.strip())
            continue

	flag = int(tokens[1])
        xs_flag = "XS:A:"
        if (flag & 16 == 0 and tokens[ts_flag_pos][5] == '+') or (flag & 16 == 16 and tokens[ts_flag_pos][5] == '-'):
            xs_flag += '+'
        elif (flag & 16 == 0 and tokens[ts_flag_pos][5] == '-') or (flag & 16 == 16 and tokens[ts_flag_pos][5] == '+'):
            xs_flag += '-'
        else:
            #print("Strange flag, will set to ? " + str(flag) + " " + token[ts_flag_pos])
            xs_flag += '?'
        tokens[ts_flag_pos] = xs_flag
        print("\t".join(tokens))
    else:
        print(l.strip())
