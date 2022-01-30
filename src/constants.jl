const J_PER_CAL = 4.184                  # This is defined as the thermochemical calorie
const JPKC = J_PER_CAL * 1000.0          # kilocalories per joule
const BOLTZMANN = 1.380658e-23           # Boltzmann's constant in J/K
const AVOGADRO = 6.0221367e+23           # Avogadro's number
const KB = (BOLTZMANN * AVOGADRO) / JPKC # Boltzmann's constant in internal units