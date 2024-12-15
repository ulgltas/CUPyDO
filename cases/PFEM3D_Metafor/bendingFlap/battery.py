from matplotlib import pyplot as plt
from natsort import natsorted
import numpy as np
import os, gmsh

# Results from Bano

plt.plot(*np.transpose(
[[0.00000, 0.000000], [0.024242, 0.006306], [0.044444, 0.015855],
[0.062626, 0.025765], [0.076767, 0.037297], [0.092929, 0.047387],
[0.117171, 0.058378], [0.139393, 0.070090], [0.169696, 0.080180],
[0.208080, 0.088828], [0.258585, 0.095495], [0.327272, 0.098558],
[0.389898, 0.098198], [0.462626, 0.097657], [0.535353, 0.096216],
[0.632323, 0.094774], [0.741414, 0.093693], [0.868686, 0.092792],
[1.000000, 0.092252], [1.185858, 0.091891], [1.381818, 0.091351],
[1.503030, 0.091351], [1.717171, 0.090990], [1.854545, 0.090990],
[2.002020, 0.090990]]))

# List the output files in workspace

gmsh.initialize()
gmsh.option.setNumber('General.Terminal', 0)
directory = natsorted(os.listdir('metafor'))

# Extract the results from Gmsh

time = list()
result = list()

for file in directory:

    gmsh.open(f'metafor/{file}')

    time.append(float(file.replace('.msh', '').replace('output_', '')))
    distance = gmsh.model.mesh.getNode(41)[0]-[0.15, 0.125, 0.1]
    result.append(np.linalg.norm(distance))

# Plot the results on the graph

plt.grid(color='lightgray')
plt.plot(time, result, '--k')
plt.savefig('plot.pdf', format='pdf')
