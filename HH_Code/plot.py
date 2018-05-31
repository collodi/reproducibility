import sys
import pathlib
import matplotlib.pyplot as plt

files = None

if len(sys.argv) > 1:
    filenames = list(pathlib.Path('.').glob('**/{}/*.m'.format(sys.argv[1].lower())))
else:
    files = list(pathlib.Path('.').glob('**/*.m'))
    for i, f in enumerate(files):
        print('%s) %s' % (i, f.as_posix()))

    filenames = [files[int(input('Choose a number. Which file? '))]]

for filename in filenames:
    data = []

    # parse data from .m file
    with open(filename.as_posix()) as f:
        for line in f.readlines():
            s = line.index('[') + 2
            e = line.index(']')

            data.append([float(x) / 1000 for x in line[s:e].split(' ')])

    xs = [x for neuron in data for x in neuron]
    ys = [i + 1 for i in range(len(data)) for _ in data[i]]

    plt.plot(xs, ys, 'r.')
    plt.xlabel('time (s)')
    plt.ylabel('Neuron #')

    if len(sys.argv) > 1:
        plt.title('{} dt=0.{}'.format(sys.argv[2], filename.name[2:-2]))
        plt.savefig(filename.name[:-2] + '.pdf')
    elif input('Save the figure? (y/n) ') == 'y':
        plt.title(input('Graph title? '))
        plt.savefig(filename.name[:-2] + '.pdf')
        plt.show()

    print('finished {}'.format(filename.name[:-2]))
    plt.close()
