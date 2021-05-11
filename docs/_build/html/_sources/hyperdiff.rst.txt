

Application of hyperdiffusion
=============================

In this section we will introduce hyperdiffusion into the equation system which we derived in the last section. We will add hyperdiffusion in two different ways, firstly as an additional diffusion and pinch and
secondly as an additional diffusion and pinch. Both these options are available in the script and they (should) yield the same results. 

First, we start with the additional diffusion and pinch. We use the flux split into two parts:

.. math::
	\Gamma^{\rho}=-D \frac{\partial n}{\partial \rho} |\vec{\nabla} \rho_t|^2+ n |\vec{\nabla} \rho_t| V


We add a hyperdiffusion by switching the diffusion: math:`\bar{D}`. i.e. :math:`D \rightarrow D +\bar{D}`. We need to also add a "hyperdiffusion" pinch so that we add zero flux:

.. math::
	\bar{V}=\bar{D}\frac{|\vec{\nabla} \rho_t|^2}{n |\vec{\nabla} \rho_t|}  \frac{\partial n}{\partial \rho}


To add hyperdiffussion in solver just: :math:`\bar{D}`. i.e. :math:`D \rightarrow D +\bar{D}` and :math:`\bar{D}`. i.e. :math:`V \rightarrow V +\bar{V}`

The other way how one can add hyperdiffusion is a bit different. This is the hyperdiffusion and an extra “source” approach. 


Now in the same way as before we add a hyperdiffusion :math:`\bar{D}`. i.e. :math:`D \rightarrow D +\bar{D}` and then add a “source” term. The “source” term is added to counter the additional flux 
from the added hyperdiffusion, but instead of using the gradient implicitly (as for the hyperdiffusion) we subtract explicitly calculated gradient with the same hyperdiffusion value. For an implicit solver:

.. math::
	\int_{i-1/2}^{i+1/2} v' \frac{\partial n}{\partial t}  d \rho &= - \left(- v' (D+\bar{D}) \frac{\partial n}{\partial \rho} |\vec{\nabla} \rho_t|^2+ v'n |\vec{\nabla} \rho_t| V \right)_{i+1/2, hat} \\
	&+\left( - v'(D+\bar{D}) \frac{\partial n}{\partial \rho} |\vec{\nabla} \rho_t|^2+v'n|\vec{\nabla} \rho_t| V \right)_{i-1/2, hat} \\
	&+ \left(- v' \bar{D} \frac{\partial n}{\partial \rho} |\vec{\nabla} \rho_t|^2 \right)_{i+1/2} -\left(- v' \bar{D} \frac{\partial n}{\partial \rho} |\vec{\nabla} \rho_t|^2 \right)_{i-1/2}+ \int_{i-1/2}^{i+1/2} v' S d\rho

"hat" denotes next time step and it is calculated at time step t + :math:`\tau` and everything else is calculated at time t. If we assume again assume equidistant radial positions we end up with this expression:

.. math::
	&\hat{n}_{i-1} \left( \frac{(D+\bar{D})_{i-1/2} v'_{i-1/2}\|\vec{\nabla} \rho|^2_{i-1/2}}{h} +\frac{V_{i-1/2} v'|\vec{\nabla} \rho|_{i-1/2}}{2} \right) \\
	+&\hat{n}_{i} \left(-\frac{(D+\bar{D})_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|^2_{i+1/2}}{h}-\frac{V_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|_{i+1/2}}{2} -\frac{(D+\bar{D})_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|^2_{i-1/2}}{h} \right. \\
	+& \left. \frac{V_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|_{i-1/2}}{2} \right) \\
	+&\hat{n}_{i+1} \left( \frac{(D+\bar{D})_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|^2_{i+1/2}}{h} -\frac{V_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|_{i+1/2}}{2} \right) = \\
	& n_{i-1} \left( \frac{\bar{D}_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|^2_{i-1/2}}{h}  \right) \\
	+& n_{i} \left( -\frac{\bar{D}_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|^2_{i+1/2}}{h}  - \frac{\bar{D}_{i-1/2} v'_{i-1/2}|\vec{\nabla} \rho|^2_{i-1/2}}{h}-\frac{\Delta V_i}{\tau} \right) \\
	+& n_{i+1} \left( \frac{\bar{D}_{i+1/2} v'_{i+1/2}|\vec{\nabla} \rho|^2_{i+1/2}}{h}  \right)  \\
	-&\int_{i-1/2}^{i+1/2} v' S d\rho
	
We can see that line 5,6 and 7 is the hyperdiffusion "source".
	






