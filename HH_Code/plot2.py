import os
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

xlimit = 8

filename = sys.argv[1]

basename = os.path.basename(filename)
dt = float('0.' + basename[2:basename.index('.')])

aaa = []
sss = []
# parse data from file
with open(filename) as f:
    for ln in f.readlines():
        d = ln.split('\t')
        aaa.append(d[2])
        sss.append(d[3])

plt.plot(aaa, 'k-', label='A', linewidth=1)
plt.plot(sss, 'r-', label='S', linewidth=1)
plt.legend()

plt.axis([0, xlimit * 1000 / dt, 0, 1])
plt.xticks([x * 1000 / dt for x in range(xlimit + 1)], range(xlimit + 1))

plt.title('80% Excitatory Using Double Precision')
plt.xlabel('time (s)')

plt.axes().set_aspect(400000)

plt.savefig(sys.argv[2], bbox_inches='tight', pad_inches=0.1)
