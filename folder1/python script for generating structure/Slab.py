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


''''''
# 创建只含有一个原子的unit cell
# primUnit = ase.build.bulk('Au','fcc', a=4.078)
# ase.visualize.view(primUnit)


######### 在slab上加CO吸附物的代码在下面 #########

# create adsorbate and put it on slab
# adsorbate = molecule('CO')
# # adsorbate = Atoms("H")
# ase.build.add_adsorbate(slab, adsorbate,2.5, position='ontop', offset=(0.5, 0.5),mol_index=0) #offset设置了H的位置
# ase.visualize.view(slab)

#### rotate CO molecure ####
# # print the position of CO before rotation
# print("\ninitial positions:")
# print(adsorbate.positions)
#
# # 假设需要将CO分子绕Z轴顺时针旋转45度
# theta = np.radians(45)  # 将角度转换为弧度
# rotation_matrix = np.array([[1, 0, 0],
#                             [0, np.cos(theta), -np.sin(theta)],
#                             [0, np.sin(theta), np.cos(theta)]])
#
# # 找到CO分子的中心位置
# center = np.mean(adsorbate.positions, axis=0)
#
# # 将CO分子的每个原子绕中心旋转
# for atom in adsorbate:
#     atom.position = np.dot(rotation_matrix, atom.position - center) + center
#
# ase.build.add_adsorbate(slab, adsorbate,2.5, position='ontop', offset=(0.5, 0.5),mol_index=0) #offset设置了H的位置
# ase.visualize.view(slab)
# # 打印旋转后CO分子的原子位置
# print("\nRotated positions:")
# print(adsorbate.positions)