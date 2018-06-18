import pandas
import argparse
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

# https://stackoverflow.com/questions/8130823/set-matplotlib-3d-plot-aspect-ratio
def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xy'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

    getattr(ax, 'set_zlim')(extents[2,0], extents[2,0] + maxsize)

def main(args):
    df = pandas.read_csv(args.file, sep=',')
    data = df[['rx','ry', 'rz']].values;
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    mpl.rcParams['legend.fontsize'] = 10

    fig = plt.figure(figsize=(10,6))
    ax = fig.gca(projection='3d')

    ax.plot(x, y, z, label='trajectory')
    ax.plot(x, y, z[0]*np.ones(len(z)), label='projection-xy', linestyle=':', color='k')
    ax.legend()

    ax.set_xlabel('$r_x$ [m]')
    ax.set_ylabel('$r_y$ [m]')
    ax.set_zlabel('$r_z$ [m]')
    plt.title('Simulation trajectory')

    axisEqual3D(ax)

    ax.view_init(elev=30., azim=225)

    plt.savefig('trajectory_3d.pdf')
    if (args.plot):
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    main(args)
