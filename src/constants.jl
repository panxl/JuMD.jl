import PhysicalConstants.CODATA2018: k_B, N_A, ε_0

const KB = ustrip(uconvert(u"kJ/mol/K", k_B * N_A))
const KE = ustrip(uconvert(u"kJ*nm", 1 / (4 * pi * uconvert(u"C^2*N^-1*m^−2", ε_0) / u"q^2")) * N_A)
