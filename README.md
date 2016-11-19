# microhybrid

The matlab folder contains all the scripts and supporting functions required to run microHybrid.m

microHybrid.m simulates the motor through an explicit timestepping method.

At each iteration we:
    Calculate the mass flow through the motor and chamber pressure through an iterative scheme
    Calculate the required power to bring the nitrous up to the preheat temperature
    Calculate the momentum and pressure contributions to the thrust
    
Simulation is stopped when the assumption that sonic flow can be obtained at the throat fails to hold
