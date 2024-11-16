import ase.build
from ase.build import molecule
from ase import Atoms
from ase.visualize import view
from ase.io import read, write
import numpy as np

### build slab model
element = 'Au'
size = (1, 1, 4)
facet = '111'
slab = ase.build.fcc111(element, size=size, vacuum=14.0, a=4.098) # insert the lattice constant after bulk structural relaxation Au 4.098; Ag 4.066;
view(slab)

### read in slab model
# element = 'Au'
# slab = read(r'E:\1Research\2024年夏\Program\Structure\Slab\Result\CONTCAR_444_Hex_Cart.vasp')
# # slab = slab.repeat((1, 2, 1))  #
# view(slab)



########## Get the geometric info ##########

# get the z-direction coordinates
z_coords = slab.get_positions()[:, 2]
# calculate the thickness of slab
slab_thickness = z_coords.max() - z_coords.min()

# get the sidelength of slab
cell = slab.get_cell()

# print
print("Cell dimensions (in Angstroms):")
print(f"a: {cell[0, 0]}")
print(f"b: {cell[1, 1]}")
print(f"c (cell height including vacuum): {cell[2, 2]}")
print(f"Slab thickness (without vacuum): {slab_thickness}")


############## save the data ##############
# Write the atoms data to the file
# write(rf'E:\1Research\2024年夏\Program\Structure\{element}({facet})_{size[0]}by{size[1]}by{size[2]}.xyz', slab)


# import ase.io
file_path = rf'E:\1Research\2024年夏\Program\Structure\Slab\{element}\{element}-({facet})_{size[0]}{size[1]}{size[2]}.vasp'
# ase.io.write(file_path, slab, format='vasp')


