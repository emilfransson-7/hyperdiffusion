.. Hyperdiffusion documentation master file, created by
   sphinx-quickstart on Mon Apr 12 15:12:22 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Hyperdiffusion proof of principle's documentation!
=============================================================

An interesting future power source is the process of fusion. The process of fusion is where two light elements are merged to create a heavier element and the difference in binding energy are converted to 
kinetic energy. This kinetic energy can be used to heat water which easily can be converted to electric power. Fusion demand the elements to be extremely energetic to overcome the Coulomb repulsion, that is why the elements be in a so called plasma state. A plasma is a hot
ionized gas which is quasi-neutral and displays collective behaviour. As the plasma is ionized the particles can be confined in a powerful magnetic field. 
One promising device to contain the plasma is the Tokamak, which is torus shape vacuum-chambers with a strong magnetic field. There are several active Tokamak experimental facilities today, such as JET, AUG, DIII-D, JT-60SA etc. 

Integrated modelling is a powerful tool to analyse todays tokamak experiments. There are several different modelling tolls in use today, ETS, ASTRA, JINTRAC among others. 

This small project with corresponding python file is a proof-of-principle for the concept of hyperdiffusion in the integrated modelling environment. This work is based on the article by Pereverzev and Corrigan.  
We will briefly present the benefits with hyperduffusion and the different ways it can be implemented in an integrated framework. This project only focuses on the particle transport and does not include
the evolution of the temperatures. In the coming sections we will describe the functionality and use of the python script. 

If you have question feel free to contact us at emil.fransson@chalmers.se

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   equations
   hyperdiff
   codefeatures
   thecode	
   	
.. #automodule:: spdec.GRF.FEM

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
