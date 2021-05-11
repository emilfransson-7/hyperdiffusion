 

Governing equations
===================

In this section we will derive the governing equation for the particle density for an integrated model. We will only derive the equations for an implicit solver as explicit solvers are not used
in modern integrated models due to their higher numerical instabilities. 

We will start with the continuity equation which describes the density evolution in the plasma:

.. math::
	\frac{\partial}{\partial t} <n> + \frac{1}{v'}\frac{\partial}{\partial \rho} ( v' <\Gamma^{\rho}>) = <S>

n is the particle density, v' is te, :math:`\rho` is the flux coordinate, :math:`\Gamma` is the particle flux and S is the internal particle source.
The brackets denote flux surface averaging. 

The particle flux can be divided into a diffusive and a pinch part.

.. math::
	<\Gamma^{\rho}>=-D \frac{\partial <n>}{\partial \rho} <|\vec{\nabla} \rho_t|^2>+ <n> <|\vec{\nabla} \rho_t|> V

We can put this in the continuity equation:

.. math::
	v' \frac{\partial <n>}{\partial t} + \frac{\partial}{\partial \rho} \left( - v' D \frac{\partial <n>}{\partial \rho} <|\vec{\nabla} \rho_t|^2>+ v'<n> <|\vec{\nabla} \rho_t|> V\right) =v' <S>
	
By integrating and using the finite difference approach we can get an equation system which we can solve. We will from now on also drop the brackets for flux surface averaging, but they
are still there.

.. math::
	&\int_{i-1/2}^{i+1/2} v' \frac{\partial n}{\partial t}  d \rho= - \left(- v' D \frac{\partial n}{\partial \rho} |\vec{\nabla} \rho_t|^2+ v'n |\vec{\nabla} \rho_t| V \right)_{i+1/2} 
	+ \left( - v' D \frac{\partial n}{\partial \rho} |\vec{\nabla} \rho_t|^2+ v'n |\vec{\nabla} \rho_t| V \right)_{i-1/2} \\
	&+ \int_{i-1/2}^{i+1/2} v' S d\rho

Let's approximate the integral on the right hand side as the volume (:math'\Delta V') multiply with the value at the centre of the volume. And use finite difference method on the derivative on the density:

.. math::
	\frac{\partial n_i}{\partial t} \Delta V_i=&D_{i+1/2} v'_{i+1/2} |\vec{\nabla} \rho|^2_{i+1/2}\frac{n_{i+1}-n_{i}}{\rho_{i+1}-\rho_{i}} - |\vec{\nabla} \rho|_{i+1/2} V_{i+1/2} v'_{i+1/2} n_{i+1/2} \\
	+& D_{i-1/2} v'_{i-1/2} |\vec{\nabla} \rho|^2_{i-1/2} \frac{n_{i}-n_{i-1}}{\rho_{i}-\rho_{i-1}}+ |\vec{\nabla} \rho|_{i-1/2} V_{i-1/2} v'_{i-1/2} n_{i-1/2} \\
	+& \int_{i-1/2}^{i+1/2} v' S d\rho

Assume equidistant grid in :math:`\rho`. :math:`h=\rho_{i+1}-\rho_{i}=\rho_{i}-\rho_{i-1}`. And that the densities between two points are an average of the points closes: 
:math:`n_{i+1/2} = (n_{i+1}+n_{i})/2`, :math:`n_{i-1/2} = (n_{i}+n_{i-1})/2`

.. math::
	\frac{\partial n_i}{\partial t} \Delta V_i=& n_{i-1} \left( \frac{D_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|^2_{i-1/2}}{h} +\frac{V_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|_{i-1/2}}{2} \right) \\
	+&n_{i} \left(-\frac{D_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|^2_{i+1/2}}{h}-\frac{V_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|_{i+1/2}}{2} -\frac{D_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|^2_{i-1/2}}{h} \right. \\
	+& \left. \frac{V_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|_{i-1/2}}{2} \right) \\
	+&n_{i+1} \left( \frac{D_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|^2_{i+1/2}}{h} -\frac{V_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|_{i+1/2}}{2} \right) \\
	+& \int_{i-1/2}^{i+1/2} v' S d\rho

The integral RHS can be approximate as:

.. math::
	\int_{i-1/2}^{i+1/2} v' S d\rho= \Delta V_i S_i    

Now we will work towards an implicit solver. This means that the gradient on the RHS is of the next timestep. This makes an explicit solver much simpler but much more unstable as well.
We now denote :math:`\hat{n}` describe the next time step, describe the current time step. math:`\tau` is the time step.

.. math::
	& \hat{n}_{i-1} \left( \frac{D_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|^2_{i-1/2}}{h} +\frac{V_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|_{i-1/2}}{2} \right) \\
	+&\hat{n}_{i} \left(-\frac{\Delta V_i}{\tau}-\frac{D_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|^2_{i+1/2}}{h}-\frac{V_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|_{i+1/2}}{2} -\frac{D_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|^2_{i-1/2}}{h} \right. \\
	+& \left. \frac{V_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|_{i-1/2}}{2} \right) \\
	+&\hat{n}_{i+1} \left( \frac{D_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|^2_{i+1/2}}{h} -\frac{V_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|_{i+1/2}}{2} \right) \\
	=& -\frac{\Delta V_i}{\tau}n_{i}- \int_{i-1/2}^{i+1/2} v' S d\rho

The equation system is the basis for an implicit solver. It doesn't include any hyperdiffussion which we will introduce in the next section.




