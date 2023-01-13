# flux_rope_analysis

This repository collects a set of scripts employed for the analysis of kink unstable flux ropes. A new version of the file `pipreadmods.py` that includes the reading of $\eta$ for resistive cases and of the neutral pressure $P_n$ is also present within this repository.

#### Energy balance and heating processes

1.  The script `energy_components.py` integrates the energy terms across the volume and allow to plot the global contribution of each component.
2.  The script `growth_rate.py` estimates the growth rate of the kink instability by fitting the logarithm of the total kinetic energy (WORK IN PROGRESS).
3.  The script `heating.py` calculates the global contributions of each heating term (Ohmic heating $\eta J^2$ for all simulations, collisional frictional heating $\frac{1}{2} \alpha_c \rho_n \rho_p v_{D}^{2}$ and ionisation-recombination frictional heating for the PIP simulations) as a function of time.
4.  The script `average_tempeature.py` calculates the average plasma and neutral temperatures of both MHD and PIP simulations, and produces a plot comparing the temperatures of different cases.

#### Current density profiles

1.  The script `current_density_frames.py` saves frames of the current density magnitude captured at the centre of the flux rope. This particular script was generated in order to produce animations.
2.  The script `current_density_mosaic.py` produces three figures:
  - The first figure is a 1D plot of the initial conditions ($t = 0$) for the magnetic field components $B_y$ and $B_z$ and the current density components $J_y$ and $J_z$ as a function of the $x$ direction.
  - The second figure is a 2D plot of the initial conditions ($t = 0$) for the magnetic field components $B_y$ and $B_z$ and the current density components $J_y$ and $J_z$ in the $xy$-plane at the centre of the flux rope.
  - The thrid figure is a mosaic of 2D frames of the current density magnitude captured at the centre of the flux rope of both MHD and PIP simulations.

#### Analysis of ionisation-recombination processes

1.  The script `ionrec.py` produces a 1D plot of the global ionization/recombination rates in PIP (partially ionized plasmas) simulations as a function of time.

