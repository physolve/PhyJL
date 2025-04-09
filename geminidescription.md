Starting from: Q (SCFH) = 63.3 * Cv * sqrt((P1(psia)² - P2(psia)²) / (SG * T(Rankine)))
Substituting the conversion factors:
Q (SLPM) = 0.4475 * 63.3 * Cv * sqrt(((14.5038 * P1(bar))² - (14.5038 * P2(bar))²) / (SG * 1.8 * T(Kelvin)))
Q (SLPM) = 28.32 * Cv * sqrt((14.5038² * (P1(bar)² - P2(bar)²)) / (SG * 1.8 * T(Kelvin)))
Q (SLPM) = 28.32 * 14.5038 / sqrt(1.8) * Cv * sqrt((P1(bar)² - P2(bar)²) / (SG * T(Kelvin)))
Q (SLPM) ≈ 306.91 * Cv * sqrt((P1² - P2²) / (SG * T))
Where P1 and P2 are in bar, T is in Kelvin, and SG is the specific gravity.


At the initial state, the pressure in Reservoir 1 (P1) is 7 bar, and the pressure in Reservoir 2 (P2) is 0.5 bar. The temperature (T) is assumed to be 293 K, and the flow coefficient (Cv) is 0.001. The specific gravity (SG) of hydrogen is 0.0695.

First, we determine the flow regime by calculating the ratio P2 / P1:
P2 / P1 = 0.5 bar / 7 bar ≈ 0.0714

Since 0.0714 is less than 0.528, the flow is initially choked. Therefore, we use the choked flow equation:

Q_initial (SLPM) = 152.98 * Cv * P1 / sqrt(SG * T)
Q_initial (SLPM) = 152.98 * 0.001 * 7 / sqrt(0.0695 * 293)
Q_initial (SLPM) = 0.15298 * 7 / sqrt(20.3335)
Q_initial (SLPM) = 1.07086 / 4.509
Q_initial ≈ 0.2375 SLPM

The initial flow rate of hydrogen from Reservoir 1 to Reservoir 2 is approximately 0.2375 standard liters per minute.