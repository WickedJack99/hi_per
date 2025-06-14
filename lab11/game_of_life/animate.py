import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def load_all_grids(folder):
    i = 0
    grids = []
    while os.path.exists(f"{folder}/step_{i}.txt"):
        grid = np.loadtxt(f"{folder}/step_{i}.txt", dtype=int)
        grids.append(grid)
        i += 1
    print(f"Loaded {len(grids)} grids from {folder}")
    return grids

filename = "output"
grids = load_all_grids(filename)

fig, ax = plt.subplots()
ax.set_aspect('equal')
img = ax.imshow(grids[0], cmap='binary')

def update(frame):
    img.set_data(grids[frame])
    ax.set_title(f"Generation {frame}")
    return img

ani = animation.FuncAnimation(fig, update, frames=len(grids), interval=150)
plt.show()