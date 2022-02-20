using JuMD
using JuMD: MMSystem, ForceGroups, CoulombForce, EwaldRecip
using Test

natoms = 2
box = (10.0, 10.0, 10.0)
cutoff = 3.0
alpha = 0.8
kmax = 11

positions = [(0.5, 0.0, 0.0), (-0.5, 0.0, 0.0)]
charges = [1.0, -1.0]
masses = [1.0, 1.0]
atomic_numbers = [1, 1]

recip = EwaldRecip(alpha, natoms, kmax, box)
elec = CoulombForce(charges=charges, cutoff=cutoff, recip=recip)
force_groups = ForceGroups(groups=(elec=elec,))
system = MMSystem(box, positions, masses, atomic_numbers, force_groups)
JuMD.force!(system)

e = (system.force_groups.energies / JuMD.KE)[1]
f = JuMD.force(system)[1][1] / JuMD.KE

# DOI: 10.1063/1.481216
e_ref = -1.0021255 
f_ref = -0.9956865

@test round(e, sigdigits=8) == e_ref
@test round(f, sigdigits=7) == f_ref
