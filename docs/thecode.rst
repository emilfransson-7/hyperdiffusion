

The functions of the script
===========================

In this section we will describe the function in the script. 

* initililization
	Creates the initial density profile

* implicit_d_and_v_solver 
	The main bulk of the script is in this function. It evolves the density profile. 
	As it is an implicit solver, the function set up an equations system and solves it. Returns the density profile, the iterations per timestep, 
	the diffusion and the pinch. 

* plotter
	As the name suggest, this plot all simulations.

* check_convergence
	A small function that check if we have convergence in the solver.

* norm_dens_der 
	Calculate the normalized density gradient.

* delta_v_func 
	Calculates the volume between two grid points.

* vprime 
	Calculates the surface area at a grid point. 

* d_hyper_source 
	This function calculate the value of the hyperdiffusion "source" if it is chosen to be used. 
	This is described in detail in the previous section.

* source 
	Calculates the internal source, in practise the Neutral Beam Injection, for the simulation. This is described by an algebraic expression in
	"source_function" and its shape mimics realistic source distribution. This function can be played around with to assure the user
	that the simulations work properly.

* pinch
	Returns the pinch to the solver. If hyperdiffusion is added as diffusion and pinch this function adds the "hyper pinch".

* diff 
	Calculates the diffusion, with or without hyperdiffusion depending on the settings. Includes the different diffusion models.

* mtanh 
	A modified tangens function used when creating the initial density profile.

* f_ped	
	Creates the initial density profile

