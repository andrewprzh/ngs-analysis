import matplotlib.pyplot as plt
import sys
import typing
from collections import OrderedDict


COLORS = ["#000000", "#990000", "#76a5af", "#674ea7", "#f65001", "#1155cc"]
SHAPES = ['o', '^', 'v', 'p', 's', 'D']
SIZES = [20 + i * 5 for i in range(7)]


def read_table(inf:typing.TextIO):
    l = inf.readline()
    tools = list(filter(None, l.strip().split('\t')))
    l = inf.readline()
    table = OrderedDict()
    row_index = 0
    while l:
        v = l.split('\t')
        label = v[0]
        table[label] = OrderedDict()
        for i in range(len(tools)):
            table[label][tools[i]] = (float(v[2*i+2]), float(v[2*i+1]), SIZES[row_index], COLORS[i], SHAPES[row_index])
        l = inf.readline()
        row_index += 1

    # print(table)
    return table


def plot_table(talbe:OrderedDict):
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    for label in talbe.keys():
        for tool in talbe[label].keys():
            point = talbe[label][tool]
            plt.scatter(x=point[0], y=point[1], s=point[2], c=point[3], marker=point[4])

    plt.show()


inf = open(sys.argv[1])
table = read_table(inf)
plot_table(table)