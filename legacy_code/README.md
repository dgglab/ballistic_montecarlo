# ballistic_montecarlo legacy code
This subrepo contains the legacy matlab code for ballistic montecarlo. It runs off of two libraries: 
'fab_feature_code' generates the device geometry, 'stable_bsweep_bandstructure' models electron transport in 
this device according to the bandstructure of the material.

# Getting started
We will walk through an example of a simple two-terminal bar geometry. Everything is run through a main file, in 
this case delafossiteBarMain.m

## stable_bsweep_bandstructure/delafossiteBarMain.m
This is the main file to run. It provides a number of inputs for various other functions that will be run to generate
various objects necessary for running the simulation. For example, we need to input a few parameters for the device
geometry, what magnetic fields to simulate, properties of the device edge (currently scattering seems to work correctly but
specularly reflecting sometimes causes issues), the number of carriers to inject, and parameters for discretizing the Fermi surface.

### Generating the device geometry
First, the device geometry is generated using 'runDelaffositeBar.' This will return a 'delafossiteBar' object which we will
later convert to a more generic causticframe. When we do this, a plot show show up of the device as will be simulated. Red boxes
will circle each ohmic contact with a number indicating which contact they are (multiple different contacts can be linked together later).
Each ohmic will be assigned an edgestyle, with 0 being a generic wall. We can for example manually set the edge style of side walls where 
we don't want to inject contacts by setting the appropriate edgestyle to 0. Finally we will set which contact we want to inject electrons from
and which contacts we want to be ground (either of which can be multiple contacts). 

### Generating the bandstructure
We will no generate two instances of the bandstructure, one for each sign of the field which require the Fermi surface to be discritized
in opposite directions. We call 'delaffositeBandstructure' which generates an instance of the bandstructure which has the necessary
methods for propagating electrons according to the bandstructure. We also rotate the crystal lattice by 'phi' relative to the device
geometry. Finally, we calculate the probability of each state on the discritized Fermi surface being injected from each edge of the device.

### Running the simulation
We are now ready to run the simulation. Everything will be passed to 'runBsweepBS'. 

### Returned data
The simulation will return a few generic things for the device geometry as well as two files for each simulated magnetic field

#### arcdata.txt
This txt file contains the (x_initial, y_initial, initial Fermi surface index, final Fermi surface index), for each of the arcs
that the electrons follow

#### classes.mat
This file will contain a few things. Namely it will contain: 'arcclass' which is the relevant bandstructure for reploting any of 
the electron trajectories, 'frmgrp' which contains the device geometry, and 'edgenum' which contains the number of times an electron
hit each ohmic of the device. 


## stable_bsweep_bandstructure/bsweep_code/calc_ohmstats.m
We have finished a run of the simulation and we now want to visualize the data. We will run calc_ohmstats on the data (use 
calc_ohmstats_directory for combining data across multiple runs). Calc_ohmstats will plot the electron flux (a proxy for voltage) 
seen by each of the contacts as a function of magnetic field. 