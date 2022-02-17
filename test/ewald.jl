using JuMD
using JuMD: MMSystem, ForceGroups, ForceGroup, CoulombForce, Ewald
using Test

natoms = 2
box = (10.0, 10.0, 10.0)
cutoff = 3.0
alpha = 0.8
kmax = 11

positions = [(0.5, 0.0, 0.0), (-0.5, 0.0, 0.0)]
masses = [1.0, 1.0]
atomic_numbers = [1, 1]

ewald = Ewald(alpha, natoms, kmax, box)
elec_force = CoulombForce([1.0, -1.0], ewald)
force_group_elec = ForceGroup(elec_force)
force_groups = ForceGroups((elec=force_group_elec,))
system = MMSystem(box, positions, masses, atomic_numbers, force_groups, cutoff=cutoff)
JuMD.force!(system)

e = (system.force_groups.energies / JuMD.KE)[1]
e_ref = -1.0021255 # DOI: 10.1063/1.481216

@test round(e, sigdigits=8) == e_ref
