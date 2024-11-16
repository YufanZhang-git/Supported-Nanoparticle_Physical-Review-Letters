# 2024-8-26
# read CHGCAR, plot 2D slice

import numpy as np
from ase.io import read
from ase.io import write
from ase.calculators.vasp import VaspChargeDensity
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.integrate import trapz

# Load the CHGCAR file
#charge = VaspChargeDensity(rf'E:\1Research\2024年夏\Program\Structure\SupNP\SupAgNP 输入-输出文件\charge density difference\CHGDIFF_cov15%_sol.vasp')

# convert the charge density from .vasp to .npy format
#charge_density = charge.chg[0]  # 3D numpy array, For non-spin-polarized calculations
#np.save(rf'E:\1Research\2024年夏\Program\Structure\SupNP\SupAgNP 输入-输出文件\charge density difference\CHGDIFF_cov15%_sol.npy', charge_density)
# Load the charge density .npy file
charge_density_vac = np.load(rf'E:\1Research\2024年夏\Program\Structure\SupNP\SupAgNP 输入-输出文件\charge density difference\CHGDIFF_cov15%_vac.npy') # (e/Å³)
charge_density_sol = np.load(rf'E:\1Research\2024年夏\Program\Structure\SupNP\SupAgNP 输入-输出文件\charge density difference\CHGDIFF_cov15%_sol.npy') # (e/Å³)

# choose vac or sol
charge_density = -1000 * charge_density_vac # (e/nm³)
# charge_density = -1000 * charge_density_sol # (e/nm³)

# Get grid dimensions
nx, ny, nz = charge_density_vac.shape
lx = 4.057 # [nm]
ly = 4.015 # [nm]
lz = 3.600 # [nm]

# Create coordinate axes
x = np.linspace(0, lx, nx) # [nm]
y = np.linspace(0, ly, ny) # [nm]
z = np.linspace(0, lz, nz) # [nm]
# X, Z = np.meshgrid(x, z) # [nm]
Y, Z = np.meshgrid(y, z) # [nm]

S_cell = lx*ly # [nm^2]
dx = x[1] - x[0] # [nm]
dy = y[1] - y[0] # [nm]
dz = z[1] - z[0] # [nm]
dV = dx * dy * dz # [nm^3]

## define the 2D slice
#y_index = round(0.45926*ny)  # You can choose any index between 0 and nz-1
x_index = round(0.5*nx)  # You can choose any index between 0 and nz-1

# Extract the 2D plane from the 3D charge density data
plane_density = charge_density[x_index, :, :]  # (e/nm³) 2D numpy array of shape (nx, ny)

# Mask definition
mask = (Z >= 1.515) & (Z <= 2) & (Y >= 1.598) & (Y <= 2.336) # [nm]
plane_density_masked = np.where(mask, plane_density.T, 0)

# Colormap settings
# Blue-red colormap
color1 = np.array([1, 5, 185]) / 255   # Blue
color2 = np.array([100, 136, 240]) / 255  # Light Blue
color3 = np.array([255, 255, 255]) / 255  # Light gray
color4 = np.array([225, 120, 90]) / 255  # Light red
color5 = np.array([165, 0, 0]) / 255    # Red
colors = [color1, color2, color3, color4, color5]

# Create a colormap with equal spacing
custom_colormap = LinearSegmentedColormap.from_list('custom_colormap', colors, N=256)

# Frame the data
upper = 6 # (e/nm³)
lower = -6 #
plane_density_clipped = np.clip(plane_density, lower, upper) # (e/nm³)
plane_density_masked_clipped = np.clip(plane_density_masked, lower, upper) # (e/nm³)

# Plotting the 2D charge density
import importlib
plt.rcParams.update({'font.size': 14})
plt.figure(1, figsize=(8, 4.75))
plt.margins(y=0)
plt.contourf(Y, Z, plane_density_clipped.T, levels=100, cmap=custom_colormap, vmin=lower, vmax=upper) # (e/nm³)
# plt.contourf(Y, Z, plane_density_masked_clipped, levels=100, cmap=custom_colormap, vmin=lower, vmax=upper) # (e/nm³)
cbar = plt.colorbar(label='Charge Density (e/nm³)')
cbar.ax.set_ylim(lower, upper)
cbar.set_ticks([lower, lower/2, 0, upper/2, upper])
cbar.set_ticklabels([lower, lower/2, 0, upper/2, upper])
plt.axis('equal')
plt.xlabel('r (nm)')
plt.ylabel('z (nm)')
plt.ylim(0, 3)
# plt.title('Charge Density on XZ-plane at Y = {:.2f}'.format(0.45926))
## # plt.savefig(rf'E:\1Research\2024年夏\Program\Structure\SupNP\SupAgNP 输入-输出文件\charge density difference\delta_rho_sol-vac_x.png', dpi=300)
plt.show()

############# 1D plot #############
# Choose the x-index (e.g., middle)
ny_range = np.array([0.331, 0.457, 0.583, 0.709])
y_range = ly * ny_range
print(y_range)
# Plotting the 2D charge density with several vertical lines
import importlib
# plt.rcParams.update({'font.size': 14})
# plt.figure(1, figsize=(8, 3.96))
# plt.margins(y=0)
# plt.contourf(X, Z, plane_density_clipped.T, levels=100, cmap=custom_colormap, vmin=lower, vmax=upper)
# plt.plot([y_range[0], y_range[0]], [0, 2.5], color='k', linestyle='-', linewidth=1) # visualize the path by which we integrate
# plt.plot([y_range[1], y_range[1]], [0, 2.5], color='m', linestyle='-', linewidth=1)
# plt.plot([y_range[2], y_range[2]], [0, 2.5], color='b', linestyle='-', linewidth=1)
# plt.plot([y_range[3], y_range[3]], [0, 2.5], color='g', linestyle='-', linewidth=1)
# cbar = plt.colorbar(label='Charge Density (e/nm³)')
# cbar.ax.set_ylim(lower, upper)
# cbar.set_ticks([lower, lower/2, 0, upper/2, upper])
# cbar.set_ticklabels([lower, lower/2, 0, upper/2, upper])
# plt.axis('equal')
# plt.xlabel('r (nm)')
# plt.ylabel('z (nm)')
# plt.ylim(0, 2.5)

y_index = round(ny*ny_range[3])  # 调整ny_range的序数

# Extract the 2D plane
line_density = plane_density[y_index, :]  # 1D numpy array of shape (1, ny)

# Plotting the 1D charge density
# upper = 60  # 6 e0/A^3
# lower = -60  #
# plt.rcParams.update({'font.size': 14})
# plt.figure(2, figsize=(6, 4))
# plt.margins(x=0)
# plt.plot([0, 2.5], [0, 0], color='k', linestyle='-', linewidth=1)
# plt.plot(z, line_density, color='g')
# plt.ylim(lower, upper)
# plt.xlim(0, 2.5)
# plt.xlabel('z (nm)')
# plt.ylabel('Charge Density (e/nm³)')

### calculate the charge density per area by integration ###
# rho_2D = np.trapz(line_density, z)  # (e/nm^2)
# print("charge density:", rho_2D)

## 一维折线图，plot 2D charge density as a function of position
# rho_2D_array = np.array([5.101, 4.169, 5.720, 6.279]) #[6.30, 4.31, 5.43, 6.64]
# plt.figure(3, figsize=(6, 3.9))
# # plt.margins(x=0)
# plt.plot(y_range, rho_2D_array, color='k', marker='o', linewidth=2, markersize=6)
# plt.ylim([0, 7])
# plt.xlim([0, 4])
# plt.xlabel('r (nm)')
# plt.ylabel('Charge Density (e/nm^2)')
# plt.show()

## #积分得到NP最上层的二维电荷密度
# integrated_y = trapz(plane_density_masked.T, y, axis=0) # [e/nm2]
# integrated_y_z = trapz(integrated_y, z) # [e/nm]
# rho_2D_masked = integrated_y_z/0.738 # [e/nm2] 0.991
# print(f"Total charge in the masked region trapz: {rho_2D_masked} e/nm2")

## 求整个cell的电荷密度和偶极矩密度

## 求整个cell平均的面电荷密度 [e/nm2]: 对整个cell的电荷密度积分，除以cell表面积，
# Perform trapezoidal integration along each axis
charge_z = np.trapz(charge_density, z, axis=2)  # [e/nm^2] Integrate over z-axis
charge_yz = np.trapz(charge_z, y, axis=1)     # [e/nm] Integrate over y-axis
charge_total_trapz = np.trapz(charge_yz, x, axis=0)   # [e] Integrate over x-axis
charge_2D_density_trapz = charge_total_trapz / S_cell  # [e/nm^2]
print(f'The average charge density of a cell is {charge_2D_density_trapz:.8f} [e/nm^2]')

## 求整个cell的dipole moment，再除以cell表面积 [e/nm]
dipole_z = np.trapz(z*charge_density, z, axis=2)  # [e/nm] Integrate over z-axis
dipole_yz = np.trapz(dipole_z, y, axis=1)     # [e] Integrate over y-axis
dipole_total_trapz = np.trapz(dipole_yz, x, axis=0)   # [e*nm] Integrate over x-axis
dipole_2D_density_trapz = dipole_total_trapz / S_cell  # [e/nm]
print(f'The average dipole moment density of a cell is {dipole_2D_density_trapz:.4f} [e/nm^2]')

plt.show()
