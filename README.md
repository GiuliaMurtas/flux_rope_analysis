# flux_rope_analysis

This repository collects a set of scripts employed for the analysis of kink unstable flux ropes.

The script "energy_components.py" integrates the energy terms across the volume and allow to plot the global contribution of each component.

The script "current_density_frames.py" saves frames of the current density magnitude captured at the centre of the flux rope.

The script "ionrec.py" produces a 1D plot of the global ionization/recombination rates in PIP (partially ionized plasmas) simulations as a function of time.

The script "heating.py" calculates the global contributions of each heating term (Ohmic heating $\eta J^2$ for all simulations, frictional heating $\frac{1}{2} \alpha_c \rho_n \rho_p v_{D}^{2}$ for the PIP simulations) as a function of time.
