# 2024-7-4
# First, read in (1) the NP model and (2) the slab model. Then, use slab.extend to combine them into SupNP

import ase.build
from ase.visualize import view
import numpy as np
from ase.cluster.wulff import wulff_construction
from ase import Atoms
from ase.build import molecule, rotate
from ase.cluster.cubic import FaceCenteredCubic
from ase.cluster.octahedron import Octahedron
from ase.io import write
from ase.io import read

############### read NP model ###############
NP_element = 'Ag'
NP = read(rf'E:\1Research\2024年夏\Program\Structure\NP\(111) terminated\{NP_element}NP_octa_trunc_half_94.xyz')#

# NP = FaceCenteredCubic(NP_element,
#                            surfaces=[(1, 0, 0), (1, 1, 1), (0, 0, -1)],#定义不同的晶面 (1, 1, -1)
#                            layers=[2, 2, 0],
#                            latticeconstant=4.086) # 4.086源自CRC Handbook

NP.rotate(45, 'z')
# view(NP)

#### get the NP info ####
Nr_atom_NP = len(NP) # NP里原子的个数
print("Number of atoms in NP:", Nr_atom_NP)

# 获取纳米颗粒的层数,为后面计算层间距
z_positions_NP = NP.positions[:, 2]  # 提取所有原子的 z 坐标
unique_z_positions_NP = np.unique(z_positions_NP)  # 获取唯一的 z 坐标值
num_layers_NP = len(unique_z_positions_NP)  # 纳米颗粒的层数

# 计算层间距
layer_spacings_NP = np.diff(unique_z_positions_NP)  # 计算相邻层之间的距离差
average_layer_spacing_NP = np.mean(layer_spacings_NP)  # 计算平均层间距
print(rf'average_layer_spacing_NP={average_layer_spacing_NP}')

############### read in the slab model ###############
slab_element = 'Au'
slab = read(rf'E:\1Research\2024年夏\Program\Structure\Slab\{slab_element}\Result\CONTCAR_954_Orth_Cart.vasp')
# view(slab)

#### get slab info ####
Nr_atom_slab = len(slab.positions) # slab.positions这个列表里的元素的数目，即slab里原子的个数
a = (slab.positions[1,1] - slab.positions[0,1])/np.sqrt(3) # Au 111面上原子间距离 2.89772; Ag
print("Number of atoms in slab:", len(slab))

# 获取slab的层数,为后面计算层间距
z_positions_slab = slab.positions[:, 2]  # 提取所有原子的 z 坐标
unique_z_positions_slab = sorted(set(z_positions_slab))  # 获取唯一的 z 坐标值
num_layers_slab = len(unique_z_positions_slab)  # 纳米颗粒的层数
# 计算层间距
layer_spacings_slab = np.diff(np.array(list(unique_z_positions_slab)))  # 计算相邻层之间的距离差
average_layer_spacing_slab = np.mean(layer_spacings_slab)  # 计算平均层间距


############### build the SupNP model ################
d_gap = (average_layer_spacing_NP + average_layer_spacing_slab)/2

# 计算d_gap，即NP最低点与slab最高点的纵向距离
slab_top_z = np.max(slab.positions[:, 2])
NP_bottom_z = np.min(NP.positions[:, 2])
d_z_gap = slab_top_z - NP_bottom_z

# 将NP移动到slab顶部，使得NP最低点与slab最高点的纵向距离为d_gap
NP.positions[:, 2] += d_z_gap + d_gap
# 将NP移动到slab中央，使得NP最低点与slab最高点的纵向距离为d_gap
NP.positions[:, 0] += 4*a # NP所有原子的x坐标移动5个a(1164) 3(744) 4(844) 4(954) 7(1484)
NP.positions[:, 1] += 4*a/2*np.sqrt(3) # NP所有原子的y坐标移动6个二分之根号三a(1164) 4(744) 4(844) 4(954) 8(1484)

slab.extend(NP) # 把移动后的NP放到slab上
view(slab) # 检查合成的supNP的形状

# 为NP建立一个cell，此cell与supNP的cell一致
NP.set_cell(slab.cell.lengths())
# view(NP)

'''
############### 建立SupNP+adsorbate模型 ################

adsorbate = molecule('CO')

#### 旋转CO分子 ####
# a1 is the vector from C to O (in the initial configuration)
a1 = tuple(adsorbate.positions[0] - adsorbate.positions[1])
# a2 is the desired final direction of the CO bond
a2 = (0, -1, np.sqrt(3))  #

# Rotate CO molecule around its center of mass (COM), rotating Z-axis to X-axis
rotate(adsorbate, a1, a2, (1, 0, 0), (1, 0, 0), center='COM')

# 把CO分子放到SupNP上
height = d_gap+num_layers_NP*d_gap+0.5
ase.build.add_adsorbate(slab, adsorbate, height, position='ontop', offset=(3.5, 1.85), mol_index=0) #offset设置了H的位置
'''

# Get the number of atoms in Supported NP
print("Number of all atoms:", len(slab))

# 保存移动后的NP文件，cell尺寸与supNP的尺寸一样
# write(rf'E:\1Research\2024年夏\Program\Structure\NP\(111) terminated\AgNP 输入-输出文件\POSCAR.vasp', NP, format='vasp') #AuNP SR的POSCAR和CONTCAR
# 保存supNP文件
# write(rf'E:\1Research\2024年夏\Program\Structure\SupNP\SupNP_{NP_element}NP{len(NP)}_Sup{len(slab)-len(NP)}.vasp', slab, format='vasp')

