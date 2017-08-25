import glob
import matplotlib.pyplot as plt

data = []
filename = glob.glob('./results/*.m')[0]

# parse data from .m file
with open(filename) as f:
    for line in f.readlines():
        s = line.index('[') + 2
        e = line.index(']')

        data.append([float(x) / 1000 for x in line[s:e].split(' ')])

xs = [x for neuron in data for x in neuron]
ys = [i + 1 for i in range(len(data)) for _ in data[i]]

plt.plot(xs, ys, 'r,')
plt.show()
