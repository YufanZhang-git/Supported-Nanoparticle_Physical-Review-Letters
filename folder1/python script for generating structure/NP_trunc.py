import numpy as np
# import ase
import ase.build
from ase.visualize import view
import ase.io
from ase.cluster.wulff import wulff_construction
from ase.io import write
from ase.cluster.cubic import FaceCenteredCubic
from ase.cluster.octahedron import Octahedron
from ase import Atoms

# Starting from here, we create our NP
# truncated 方法1, 可创建SupNP，通过在layers里设置-1
element = 'Ag'
NP = FaceCenteredCubic(element,
                       surfaces=[(1, 0, 0), (1, 1, 1), (1, 1, -1)],  #, (0, 0, -1) 定义不同的晶面 (1, 1, -1)
                       layers=[4, 3, -1],  #, 2
                       latticeconstant=4.0857)  # Au实验值CRC：4.078, Ag实验值CRC: 4.0857

# NP = Octahedron(element,
#                            length=7,
#                            cutoff=2, # cuboctahedral: length=2*cutoff+1; regular octahedron: length=3*cutoff+1
#                            latticeconstant=4.085)

# NP.rotate(10, 'z')
NP.rotate([1, 1, -1], [0, 0, -1])
# view(NP)


#################### get NP info ##################
# Get the number of atoms
print("Number of atoms:", len(NP))

# Get the positions of all atoms
positions = NP.get_positions()

# Calculate the distance between all pairs of atoms
distances = np.linalg.norm(positions[:, np.newaxis] - positions, axis=2)

# Find the maximum distance (diameter)
print("Diameter of the atoms model:", np.max(distances))

# get the number of layers of the NP
z_positions = NP.positions[:, 2]  # get all the z coordinates
unique_z_positions = np.unique(np.around(z_positions, 5))  # 获取唯一的 z 坐标值
num_layers_NP = len(unique_z_positions)  # 纳米颗粒的层数
print("unique_z_positions:", unique_z_positions)
print("Number of layers in NP:", num_layers_NP)

# calculate the 层间距
layer_spacings = np.diff(unique_z_positions)  # 计算相邻层之间的距离差
average_layer_spacing = np.mean(layer_spacings)  # 计算平均层间距
print("Average layer spacing in NP:", average_layer_spacing)

## the following codes manually add one more layer of atoms to the created NP
# The x, y coordinates are read from VESTA by putting the cursor on the corresponding atom
new_atom = Atoms('Ag', positions=[(-3.970, -1.927, -1.180+3*average_layer_spacing)])
NP += new_atom
new_atom = Atoms('Ag', positions=[(-3.223, 0.863, -1.180+3*average_layer_spacing)])
NP += new_atom
new_atom = Atoms('Ag', positions=[(-2.475, 3.654, -1.180+3*average_layer_spacing)])
NP += new_atom
new_atom = Atoms('Ag', positions=[(0.316, 4.402, -1.180+3*average_layer_spacing)])
NP += new_atom
new_atom = Atoms('Ag', positions=[(2.359, 2.359, -1.180+3*average_layer_spacing)])
NP += new_atom
new_atom = Atoms('Ag', positions=[(4.402, 0.316, -1.180+3*average_layer_spacing)])
NP += new_atom
new_atom = Atoms('Ag', positions=[(3.654, -2.475, -1.180+3*average_layer_spacing)])
NP += new_atom
new_atom = Atoms('Ag', positions=[(0.863, -3.223, -1.180+3*average_layer_spacing)])
NP += new_atom
new_atom = Atoms('Ag', positions=[(-1.927, -3.970, -1.180+3*average_layer_spacing)])
NP += new_atom
new_atom = Atoms('Ag', positions=[(-0.432, 1.611, -1.180+3*average_layer_spacing)])
NP += new_atom
new_atom = Atoms('Ag', positions=[(-1.180, -1.180, -1.180+3*average_layer_spacing)])
NP += new_atom
new_atom = Atoms('Ag', positions=[(1.611, -0.432, -1.180+3*average_layer_spacing)])
NP += new_atom

view(NP)
#################### Write the atoms data to the file #####################
# write(rf'E:\1Research\2024年夏\Program\Structure\NP\(111) terminated\{element}NP_octa_trunc_half_{len(NP)}.xyz', NP) # , format='vasp'


#####
# This code block is an example for Ru NP. We do not really use it
# atoms = wulff_construction('Ru',
#                            surfaces=[(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 0), (2, 1, 1)],#定义不同的晶面
#                            energies=[0.180, 0.190, 0.160, 0.200, 0.190],#定义不同晶面的表面能，需和晶面一一对应
#                            size=100, #
#                            structure='fcc', #
#                            latticeconstant=3.7476, #
#                            rounding='closest') #or, 'above' 'below'
#