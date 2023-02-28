from numpy.random import RandomState
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
from mpl_toolkits.mplot3d import *
import random

import mpl_toolkits.mplot3d
import HashPrimerIndexMatrix.HashUtils as hashUtils
from HashPrimerIndexMatrix.BaiscObjects import HashPrimerItem, HashPrimerMatrix
from matplotlib.animation import FuncAnimation


def plot_3d(points, points_color, title)-> tuple:
    x, y, z = points.T

    fig, ax = plt.subplots(
        figsize=(12, 12),
        facecolor="white",
        tight_layout=True,
        subplot_kw={"projection": "3d"},
    )
    fig.suptitle(title, size=16)
    ax.scatter(x, y, z, c=points_color, s=500)
    ax.view_init(azim=-60, elev=9)

    # x = [0.5, 1.5]
    # y = [0.5, 1.5]
    # z = [0.5, 1.5]
    # u = [1, 1]
    # v = [1, -1]
    # w = [1, -1]
    # ax.quiver(x, y, z, u, v, w, length=1.732/2, arrow_length_ratio=0.3)
    # fig.colorbar(col, ax=ax, orientation="horizontal", shrink=0.6, aspect=60, pad=0.01)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.zaxis.set_major_locator(ticker.MultipleLocator(1))

    return fig, ax

def anim(_hashPrimerMatrix: HashPrimerMatrix):
    ax.clear()
    points = create_3dim_blocks_array(PATTERN)
    x, y, z = points.T
    points_color = create_colmap_array_by_hashPrimerMatrix(_hashPrimerMatrix)
    ax.scatter(x, y, z, c=points_color, s=250)
    set_3d_axies(ax, PATTERN)

def set_3d_axies(ax: Axes3D, pattern):

    ax.set_xlim3d(0, pattern[0])
    ax.set_ylim3d(0, pattern[1])
    ax.set_zlim3d(0, pattern[2])
    ax.set_xticks([0.5+i for i in range(0, pattern[0], 5)])
    ax.set_xticklabels([i for i in range(0, pattern[0], 5)], rotation=40)

    ax.set_yticks([0.5+i for i in range(pattern[1])])
    ax.set_yticklabels([i+1 for i in range(pattern[1])])

    ax.set_zticks([0.5 + i for i in range(pattern[2])])
    ax.set_zticklabels([str((i + 1)*100)+"×"+str((i + 1)*100) for i in range(pattern[2])])
    ax.set_xlabel("FileName Dimension", labelpad=15)
    ax.set_ylabel("Series Dimension", labelpad=15)
    ax.set_zlabel("Resolution Dimension", labelpad=20)

def plot_2d(points, points_color, title):
    fig, ax = plt.subplots(figsize=(3, 3), facecolor="white", constrained_layout=True)
    fig.suptitle(title, size=16)
    add_2d_scatter(ax, points, points_color)
    plt.show()


def add_2d_scatter(ax, points, points_color, title=None):
    x, y = points.T
    ax.scatter(x, y, c=points_color, s=50, alpha=0.8)
    ax.set_title(title)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())


# 创建一个三维矩阵，用于绘制3D图
def create_3dim_blocks_array(pattern):
    x, y, z = pattern[0], pattern[1], pattern[2]
    array = np.ones((x*y*z, 3), dtype=np.float)
    m = 0
    for i in range(x):
        for j in range(y):
            for k in range(z):
                array[m] = [i+0.5, j+0.5, k+0.5]
                m += 1
    return array

def create_colmap_array(pattern: tuple, red_list, blue_list):
    x, y, z = pattern[0], pattern[1], pattern[2]
    cmap = ["#00aa0000"] * (x*y*z)

    for position in blue_list:
        cmap[position[0] * (y * z) + position[1] * (z) + position[2]] = "#0000aaa0"
    for position in red_list:
        cmap[position[0] * (y*z) + position[1] * (z) + position[2]] = "#ff0000a0"


    return cmap

def create_colmap_array_by_hashPrimerMatrix(hashPrimerMatrix: HashPrimerMatrix):
    x, y, z = hashPrimerMatrix.pattern[0], hashPrimerMatrix.pattern[1], hashPrimerMatrix.pattern[2]
    cmap = ["#00aa0020"] * (x * y * z)
    for item in hashPrimerMatrix.items:
        if item.getFile() is not None:
            index = item.getCoordinateToColorMapIndex()
            cmap[index] = "#0000aa50"
        if item.checkIfInOriginalCoordinate() is False:
            index = item.getCoordinateToColorMapIndex()
            cmap[index] = "#ff000080"

    return cmap

PATTERN = (8, 8, 8)
HashPrimerItem.pattern = PATTERN
hashPrimerMatrix = HashPrimerMatrix(PATTERN, True)

with open("20merPrimerLibraryFinal.txt", 'r') as files:
    files = eval(files.read())

random.seed(195)
#files = files[:500]
files = [str(i) for i in range(40)]
ress = [random.randint(0, PATTERN[2]-1) for _ in range(len(files))]
series = [random.randint(0, PATTERN[1]-1) for _ in range(len(files))]
anims = []
from copy import deepcopy
for f, r, s in zip(files, ress, series):
    hash_primer_matrix = hashPrimerMatrix.trySaveFile(f, r, s)
    anims.append(deepcopy(hash_primer_matrix))

S_color = create_colmap_array_by_hashPrimerMatrix(hashPrimerMatrix)
fig, ax = plot_3d(create_3dim_blocks_array(PATTERN), S_color, "Primer Storage Array")
set_3d_axies(ax, PATTERN)


animation = FuncAnimation(fig, anim, frames=anims, repeat=False, blit=False, interval=100)
#animation.save(str(PATTERN)+".gif", writer='pillow')

plt.tick_params(pad=10)
plt.show()
