#!/usr/bin/env python3

# %%
import numpy as np
import matplotlib.pyplot as plt

import meep as mp
from meep.materials import Si, aSi, cSi

# %%
# r / a = 0.35
# ωa / (2πc) = fa / c = 0.225 -> 0.325

r_over_a = 0.35
omega_a_over_2pic = 0.29

a = 0.4
r = a * r_over_a

f = 1 / a * omega_a_over_2pic

wl = 1 / f

# material = Si
material = mp.Medium(epsilon=Si.epsilon(f)[0, 0])

# %%
sx = 20 * a
# sy = 12 * np.sqrt(3) * a
sy = sx

dpml = sx * 0.1

cell_size = mp.Vector3(sx + 2 * dpml, sy + 2 * dpml)

resolution = 100

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

for x1 in x_range:
    for y1 in x_range:
        # x1, y1 : square; x2, y2 : hex

        x2, y2 = np.dot(square2hex_mat, np.array([x1, y1]))

        if x2 >= sx / 2:
            x2 -= sx
        if x2 <= -sx / 2:
            x2 += sx

        if np.isclose(y2, 0):
            continue

        geometry.append(mp.Cylinder(
            r,
            center=mp.Vector3(x2, y2),
        ))

sources = [
    mp.Source(
        mp.ContinuousSource(frequency=f),
        component=mp.Ez,
        center=mp.Vector3(-sx / 2 + 1 * a, 0),
    ),
]

pml_layers = [mp.PML(dpml)]

# %%
sim = mp.Simulation(
    cell_size=cell_size,
    boundary_layers=pml_layers,
    geometry=geometry,
    sources=sources,
    resolution=resolution,
)

# %%

# %%
# sim.plot2D()
# plt.show()

# %%
# sim.run(until=40)

# print(sim.get_array(mp.Ez).max())
# print(sim.get_array(mp.Ez).min())

# %%
# fig, ax = plt.subplots(dpi=100)
# sim.plot2D(
#     ax,
#     fields=mp.Ez,
#     output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)),
# )
# plt.show()

# %%

E_t = []

field_energy = lambda sim: E_t.append(
    sim.field_energy_in_box(
        center=mp.Vector3(sx / 2 - 0.5 * a, 0),
        size=mp.Vector3(0,
                        np.sqrt(3) * a + 0 * r),
    )
)

# %%
mp.verbosity(1)
sim.reset_meep()

sim.use_output_directory("out")
sim.run(
    mp.at_beginning(mp.in_volume(mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)), mp.output_epsilon)),
    until=5,
    # mp.at_every(
    #     0.1,
    #     mp.in_volume(
    #         mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)),
    #         mp.to_appended("ez", mp.output_efield_z),
    #     ),
    # ),
    # mp.after_time(40, mp.at_every(0.01, field_energy)),
    # until=50,
)

print(sim.get_array(mp.Ez).max())
print(sim.get_array(mp.Ez).min())

print(np.sum(E_t))

# %%
