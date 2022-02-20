using JuMD
using JuMD: MMSystem
using Test

system = MMSystem("test/ala_di.parm7", "test/ala_di.rst7")
JuMD.force!(system)

e = system.force_groups.energies

# account for the different Coulomb constant from AmberTools
KE_AMBER = 332.0522173
S = (JuMD.KE / 4.184 * 10) / KE_AMBER

# calcuated from sander's Python API
# energies were converted from kcal/mol to kJ/mol
e_ref = [0.020598314980793887,
         0.3619937665763452,
         9.643998432066596,
         2.8119866647276215,
         5.015691768207229,
         -80.12379878458636 * S,
         48.93546385709723 * S] .* 4.184

@test e ≈ e_ref
