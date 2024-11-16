# 20240703

import ase.build
from ase.build import surface
from ase.visualize import view
from ase.io import write
import numpy as np

# 创建slab模型
element = 'Au'
layers = 4
size = (1, 1, layers)
facet = '111'
slab = ase.build.fcc111(element, size=size, vacuum=14.0, a=4.098) # 4.098是经过bulk struc relax后得到的
# view(slab)
slab = slab.repeat((1, 2, 1))  # 重复以形成正交晶胞
# view(slab)

# 调整晶胞使其成为正交晶胞
cell = slab.get_cell()
cell[0, 1] = 0.0  # 将非对角元设为0以形成正交晶胞
cell[1, 0] = 0.0
slab.set_cell(cell, scale_atoms=False)

# # 可视化和保存模型
# view(slab)
# write(rf'E:\1Research\2024年夏\Program\Structure\Slab\{element}-({facet})_法3.vasp', slab, format='vasp') # {size[0]}{size[1]}{size[2]}

# 扩大
n1 = 11
n2 = 6
slab = slab.repeat((n1, n2, 1))  # 重复以形成正交晶胞
view(slab)
print(f"Total number of atoms in the slab: {len(slab)}")
print(f"Length of the slab in the a direction: {slab.cell.cellpar()[0]:.2f} Å")
print(f"Length of the slab in the b direction: {slab.cell.cellpar()[1]:.2f} Å")

# # 可视化和保存模型
# write(rf'E:\1Research\2024年夏\Program\Structure\Slab\{element}-({facet})_{n1}-{2*n2}-{layers}.vasp', slab, format='vasp') # {size[0]}{size[1]}{size[2]}
