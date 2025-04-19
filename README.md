# NVSim

A program to simulate the resonance of nitrogen-vacancy centers and compare them with experimental ODMR measurements.
It can find the orientation of the NV centers and the magnetic field strength using a single microwave frequency sweep.

The repository includes a web app that accepts input data in JSON format and returns the reconstructed magnetic field and other results (uncertainty, trajectory, fit error).

The project also includes a number of Pluto Notebooks for interactive use, the most important one being [Iterative.jl](Notebooks/Iterative/Iterative.jl) which contains the iterative algorithm.
