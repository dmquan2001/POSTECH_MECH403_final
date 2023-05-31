#!/usr/bin/env python

# %%
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

# %%
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
os.getcwd()

# %%
# E_file = 'E_rp_rn.dat'
# r_p_file = 'r_p.dat'
# r_n_file = 'r_n.dat'

# E = np.loadtxt(E_file)
# r_p = np.loadtxt(r_p_file)
# r_n = np.loadtxt(r_n_file)

# r_p = np.linspace(0, 1, 6)
# r_n = np.linspace(0, 1.4, 7)

# rp_rn_E = np.empty((0, 3), float)

# for i in range(len(r_p)):
#     for j in range(len(r_n)):
#         # print(r_p[i], r_n[j], E[i, j])
#         rp_rn_E = np.vstack([rp_rn_E, [r_p[i], r_n[j], E[i, j]]])

# np.savetxt("rp_rn_E.dat", rp_rn_E)

# %%
x_y_E_files = [
    'rp_rn_E_0.dat',
    'rp_rn_E_1.dat',
    'rp_rn_E_2.dat',
    'rp_rn_E_3.dat',
    'rp_rn_E_4.dat',
    'rp_rn_E_5.dat',
]

x_y_E = np.vstack([np.loadtxt(f) for f in x_y_E_files])
x_y_E = np.unique(x_y_E, axis=0)

x = x_y_E[:, 0]
y = x_y_E[:, 1]
E = x_y_E[:, 2] / 599.2415764262199

x_max = x[E == E.max()][0]
y_max = y[E == E.max()][0]
E_max = E.max()
print(x_max, y_max, E_max)

np.savetxt('rp_rn_E_all.dat', x_y_E)

# %%

xi = np.linspace(x.min(), x.max(), 201)
yi = np.linspace(y.min(), y.max(), 201)

triang = mpl.tri.Triangulation(x, y)
interpolator = mpl.tri.LinearTriInterpolator(triang, E)
Xi, Yi = np.meshgrid(xi, yi)
Ei = interpolator(Xi, Yi)

# vmin = E.min()
# vmax = E.max()

vmin = 0
vmax = 0.85
n_levels = 50
levels = np.linspace(vmin, vmax, n_levels)

fig, ax = plt.subplots(dpi=300)
# plt.axis('off')
plt.gca().set_aspect('equal')
# plt.gca().invert_yaxis()
plt.xlabel(r"$r_A / r$")
plt.ylabel(r"$r_B / r$")
plt.title("Percent of power transmitted")

cntr = ax.contourf(xi, yi, Ei, n_levels, cmap='viridis', vmin=vmin, vmax=vmax)
# cntr = ax.contourf(xi, yi, p_angle_i, 100, cmap='hsv', vmin=vmin, vmax=vmax, levels=levels)
# cntr = ax.contourf(xi, yi, p_angle_i, levels=100, cmap='hsv')

cbar = fig.colorbar(
    mpl.cm.ScalarMappable(norm=cntr.norm, cmap=cntr.cmap),
    ticks=np.array([0, 0.2, 0.4, 0.6, 0.85]),
    boundaries=levels,
    values=(levels[:-1] + levels[1:]) / 2,
)

# plt.scatter(x, y)
plt.show()

# %%
