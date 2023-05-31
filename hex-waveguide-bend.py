#!/usr/bin/env python3

# %%
import numpy as np
import matplotlib.pyplot as plt

import meep as mp
from meep.materials import Si, aSi, cSi
from tqdm import tqdm

# %%


def gen_Simulation(r_p, r_n, parameters):

    r_over_a, omega_a_over_2pic, resolution, a = parameters

    r = a * r_over_a
    f = 1 / a * omega_a_over_2pic
    wl = 1 / f

    sx = 20 * a
    sy = sx

    material = Si
    # material = mp.Medium(epsilon=Si.epsilon(freq=f)[0, 0])

    dpml = sx * 0.1

    cell_size = mp.Vector3(sx + 2 * dpml, sy + 2 * dpml)

    geometry = []

    geometry.append(
        mp.Block(
            size=mp.Vector3(sx + 2 * dpml, sy + 2 * dpml, mp.inf),
            center=mp.Vector3(),
            material=material,
        )
    )

    square2hex_mat = np.array([[1, 0.5], [0, np.sqrt(3) / 2]])

    x_range = np.arange(-sx / 2, sx / 2 + a / 2, a)

    for i, x1 in enumerate(x_range):
        for j, y1 in enumerate(x_range):
            # x1, y1 : square; x2, y2 : hex

            x2, y2 = np.dot(square2hex_mat, np.array([x1, y1]))

            if x2 >= sx / 2:
                x2 -= sx
            if x2 <= -sx / 2:
                x2 += sx

            p_defects = False
            n_defects = False

            n = len(x_range)
            n_mid = int(n / 2)

            p_defect1 = i == n_mid and j == n_mid
            p_defects = p_defect1

            n_defect2 = i == n_mid - 1 and j == n_mid + 1
            n_defects = n_defect2

            if p_defects:
                geometry.append(mp.Cylinder(
                    radius=r * r_p,
                    center=mp.Vector3(x2, y2),
                ))
                continue

            if n_defects:
                if np.isclose(r_n, 0):
                    continue
                else:
                    geometry.append(mp.Cylinder(
                        radius=r * r_n,
                        center=mp.Vector3(x2, y2),
                    ))
                continue

            if np.isclose(y2, 0) and x2 < 0:
                continue

            if np.isclose(x1, 0) and y2 > 0:
                continue

            geometry.append(mp.Cylinder(
                r,
                center=mp.Vector3(x2, y2),
            ))

    sources = [
        mp.Source(
            mp.ContinuousSource(wavelength=wl),
            component=mp.Ez,
            center=mp.Vector3(-sx / 2 + 1 * a, 0),
        ),
    ]

    pml_layers = [mp.PML(dpml)]

    sim = mp.Simulation(
        cell_size=cell_size,
        boundary_layers=pml_layers,
        geometry=geometry,
        sources=sources,
        resolution=resolution,
    )
    return sim


def cal_E(sim, parameters, plot=False):
    r_over_a, omega_a_over_2pic, resolution, a = parameters

    r = a * r_over_a
    f = 1 / a * omega_a_over_2pic
    wl = 1 / f

    sx = 20 * a
    sy = sx

    E_t = []

    field_energy = lambda sim: E_t.append(
        sim.field_energy_in_box(
            center=mp.Vector3(5 * a, 10 * np.sqrt(3) / 2 * a),
            size=mp.Vector3(np.sqrt(3) * a + 2 * r, 0),
        )
    )

    sim.reset_meep()
    sim.run(
        mp.after_time(40, mp.at_every(0.01, field_energy)),
        until=50,
    )

    if plot:
        plt.figure(dpi=100)
        plt.plot(E_t)
        plt.show()

    return np.sum(E_t)


# %%


def main_single(parameters, r_p, r_n):
    sim = gen_Simulation(r_p, r_n, parameters)
    E = cal_E(sim, parameters)
    print(E)


def main_loop(parameters, r_p_, r_n_):

    ni = len(r_p_)
    nj = len(r_n_)
    E_rp_rn = np.zeros((ni, nj))

    mp.verbosity(0)
    with tqdm(total=ni * nj) as pbar:
        for i in range(ni):
            for j in range(nj):
                # print(f'Progress: {i}/{ni}, {j}/{nj}')
                r_p = r_p_[i]
                r_n = r_n_[j]
                # print(r_p, r_n)

                sim = gen_Simulation(r_p, r_n, parameters)
                E_rp_rn[i, j] = cal_E(sim, parameters)
                # print(f'Progress: {i}/{ni}, {j}/{nj}')

                pbar.update(1)

    # np.savetxt('out/E_rp_rn.dat', E_rp_rn)
    # np.savetxt('out/r_p.dat', r_p_)
    # np.savetxt('out/r_n.dat', r_n_)

    rp_rn_E = np.empty((0, 3), float)

    for i in range(len(r_p_)):
        for j in range(len(r_n_)):
            r_p = r_p_[i]
            r_n = r_n_[j]
            # print(r_p[i], r_n[j], E[i, j])
            rp_rn_E = np.vstack([rp_rn_E, [r_p, r_n, E_rp_rn[i, j]]])

    np.savetxt("out/rp_rn_E.dat", rp_rn_E)


# %%
# r / a = 0.35
# ωa / (2πc) = fa / c = 0.225 -> 0.325

r_over_a = 0.35
omega_a_over_2pic = 0.29

resolution = 100
a = 0.4

r = a * r_over_a
f = 1 / a * omega_a_over_2pic
wl = 1 / f

sx = 20 * a
sy = sx
dpml = sx * 0.1
r = a * r_over_a

parameters = [r_over_a, omega_a_over_2pic, resolution, a]

# %%
# r_p = 0.66
# r_n = 0.86
# main_single(parameters, r_p, r_n)

# %%
# r_p_ = np.linspace(0, 1.4, 14 * 2 + 1)
# r_n_ = np.linspace(0, 1, 41)

# r_p_ = np.array([0, 0.2, 0.4, 0.5, 0.55, 0.58, 0.6, 0.62, 0.65, 0.8, 1, 1.2])
# r_n_ = np.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 1])

# r_p_ = np.array([0.54, 0.56, 0.58, 0.6, 0.62, 0.64])
# r_n_ = np.array([0.9, 0.91, 0.92, 0.93, 0.94, 0.95])

# r_p_ = np.array([0.62, 0.64, 0.66, 0.68, 0.7])
# r_n_ = np.array([0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91])

# r_p_ = np.array([0.64, 0.66, 0.68, 0.7])
# r_n_ = np.array([0.80, 0.81, 0.82, 0.83, 0.84])

# r_p_ = np.array([0.54, 0.56, 0.58, 0.6, 0.62])
# r_n_ = np.array([0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89])

# r_p_ = np.array([0.66, 0.68, 0.70])
# r_n_ = np.array([0.91, 0.92, 0.93, 0.94, 0.95])

# r_p_ = np.linspace(0, 1, 6)
# r_n_ = np.linspace(0, 1.2, 7)

# main_loop(parameters, r_p_, r_n_)

# %%

# sim.plot2D()
# plt.show()

# %%
# fluxr = mp.FluxRegion(
#     center=mp.Vector3(5 * a, 10 * np.sqrt(3) / 2 * a),
#     size=mp.Vector3(np.sqrt(3) * a + 2 * r, 0),
# )
# sim.add_flux(1, 1, 1, fluxr)
# sim.plot2D()
# plt.show()

# %%
# r_p = 0.66
# r_n = 0.85
# sim = gen_Simulation(r_p, r_n, parameters)
# sim.run(until=50)

# %%
# Ez = sim.get_array(
#     center=mp.Vector3(),
#     size=mp.Vector3(sx, sy),
#     component=mp.Ez,
# )
# eps = sim.get_array(
#     center=mp.Vector3(),
#     size=mp.Vector3(sx, sy),
#     component=mp.Dielectric,
# )

# %%
# fig, ax = plt.subplots(dpi=100)

# ax.set_title(r"Scalar map of $E_z$")
# ax.set_xlabel(r"$X\, (\mu \mathrm{m})$ ")
# ax.set_ylabel(r"$Y\, (\mu \mathrm{m})$ ")

# ax.imshow(eps.transpose(), interpolation='spline36', cmap='binary')
# ax.imshow(Ez.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9, vmin=-1, vmax=1)

# # plt.savefig("Ez_max.svg")
# plt.show()

# %%
r_p = 0
r_n = 1
# r_p = 0.66
# r_n = 0.85
sim = gen_Simulation(r_p, r_n, parameters)

# %%
mp.verbosity(0)
sim.reset_meep()

sim.use_output_directory("out")
sim.run(
    # mp.at_beginning(mp.in_volume(mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)), mp.output_epsilon)),
    # until=1,
    mp.at_every(
        0.1,
        mp.in_volume(
            mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)),
            mp.to_appended("ez", mp.output_efield_z),
        )
    ),
    until=50,
)

print(sim.get_array(mp.Ez).max())
print(sim.get_array(mp.Ez).min())

# %%
