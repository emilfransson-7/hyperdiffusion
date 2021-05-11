

Script features
===============

In this section we will introduce different features of the hyperdiffusion script. First, we will discuss the more general information the script needs to work properly. 
The script has three major inputs which the script loop over all iterations of. The three inputs are:

* rho_pos_array
	Define the number of points in the radial grid. Can use multiple values as input, the script loop through all input. 

* time_step_array 
	Timestep in [s]. Can use multiple values as input, the script loop through all input. 

* hyper_diff_array 
	Value for the hyperdiffusion if one uses it. Can use multiple values as input, the script loop through all input. 

The script will do a simulation with each instant of these three inputs. It will later save the output (only the last time step) as a simple .txt-file.
Be aware if you have done a simulation with the same input it will be overwritten. The script will also print the time evolution of the density profile. A typical output can look like this:

.. image:: https://am3pap007files.storage.live.com/y4mOLfMRlg0-83WthxD5bNOrti0QJEx1RzrR0Y0GX8IHd1uMsz91mCL4O93SPoFh31e6pOJC8k6WGhLTnA4bHpmHoIT-zzt3YXjkhEBgh5vHWgUzcC-Ax-62KKg4mfI5rO1UdXuOGiGiRifUXIzRgvLrYuy1uXJzdxjFvAnQJbv7RyO360LUMTJHNceoDMZWxZr?width=640&height=480&cropmode=none

In the figure we can see the initial and the evolution at 5 equidistant timepoints. We can notice that the profile converges towards the models 
quasi steady state.

The script also print a couple of things in the terminal which can be seen in the figure below:

.. image:: https://am3pap007files.storage.live.com/y4m60t0LxnAnJwSPUlLJq2G7ov7z9E5oJDFylAhBjuiQVoqVlL-h_j51FQ2z1MP3915T2PqIC4qrXf0w_xtt7fjopecayOX499paOrPotDi0YhFR1OTQucHBS0DoixEjZ2gTlvs6mm9ZCJD7_DXMeSIklL1Hc-J9VY6wsMFtISLzVM5l2WyIXb7a4TyLPTyqSB6?width=525&height=200&cropmode=none

This was a run with three simulation and the script indicates when each of them is done. The output also denotes when if it overwrites previously stored data. 
It also indicates of long the run took.

There are a couple of other inputs that the user might want to change. We will give a short description of them here:

* t_end 
	Define the end time of the simulation. (The simulation starts at 0). An easy way to see if t_end is large enough is to look if the profiles have saturated to a steady-state, 
	as the flux is dependent on the gradient and the source. At a certain gradient the flux and the source will be equal, and the profiles have reach a so called "quasi steady state". 
	It might take some time to reach there, however.  

* t_save
	Define of often the script store the density profile.

* rho_tor_cut_of 
	Defines the normalized radial position outside where we do not evolve the profiles

* d_hyperdiff_chooser
	Define how the hyperdiffusion is added or not. 0 means that we do not use hyperdiffusion. 1 we add hyperdiffusion as diffusion and pinch, 2 as a pinch and a "source" term. 
	We discussed both these methods in the previous section.

* choose_d
	Choose of how the diffusion is calculated. There are four different type of transport "models" to calculate the diffusion which the script uses. The first is just a constant diffusion, the second is a model quadratically proportional to the
	the normalized density gradient. The third is like the second one with a constant diffusion. The fourth, and the one recommended to use, is a mimic of the IFSPPL model. This gives the effect of the 
	Ion Temperature Gradient - mode (which is usually the dominant mode) and it’s dependency on a critical (density) gradient. Below this value we have a very small diffusion but above the diffusion is 
	rapidly increasing with the normalized gradient.    

* choose_convergence_factor 
	Defines the accuracy for the convergence. 0.0001 is a reasonable value. The simulations usually need more iteration in the beginning, for later time steps the simulations usually 
	only need one iteration. It is an average over the whole radial profile. 

* max_iterations 
	Determines the maximum number of iterations in the convergence loop. A value of 10 is usually a high enough value. 
	The number of iterations is usually largest for the first time steps.

* b_pos, b_height, b_sol, b_width, b_slope
	Defines the density initial profile. 

A word of the pinch we use in the script. We only use a constant value and the default is 0.5 [m/s]. This is a very crude "model" but as we only evolve the density and the pinch is not dependent of the density gradient it is 
sufficient for our script.

	






