import pathlib
import matplotlib.pyplot as plt

data = []
files = list(pathlib.Path('.').glob('**/*.m'))

for i, f in enumerate(files):
    print('%s) %s' % (i, f.as_posix()))

filename = files[int(input('Choose a number. Which file? '))]

# parse data from .m file
with open(filename.as_posix()) as f:
    for line in f.readlines():
        s = line.index('[') + 2
        e = line.index(']')

        data.append([float(x) / 1000 for x in line[s:e].split(' ')])

xs = [x for neuron in data for x in neuron]
ys = [i + 1 for i in range(len(data)) for _ in data[i]]

plt.plot(xs, ys, 'r.')
if input('Save the figure? (y/n) ') == 'y':
    plt.savefig(filename.name[:-2] + '.pdf')

plt.show()
