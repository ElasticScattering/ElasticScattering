# Overview

This repository serves two purposes.

The first purpose is to share the C++/OpenCL code used to compute the real space impurity scattering. Here, we share both the full code (see the README how to build the code and the dependencies required) as well as an executable version (see XXXXXX)

The second purpose is to show the python code used to generate the graphs throughout the paper 'B^2 to B-linear magnetoresistance due to impeded orbital motion'. This includes calculations for the simple model shown in the main text, as well as transforming the data exported by the elastic impurity scattering program into figures for the supplementary material. The python files are themselves executables, the dependencies are given below.

# Elastic Scattering

![T3 MF0](https://user-images.githubusercontent.com/25907608/170357566-d476233d-febd-4996-9374-0275faa2c3e2.png)

Program to computes the lifetime of particles scattering around in a field of impurities.

The project can be run by giving it a config file, for example: `cl-es.exe examples/small.config`. The config allows you to tune various parameters of the simulation, such as the magnetic field strength and temperature, and the scale of the simulation. 

The results of a simulation are put inside a folder, with optional visualizations of intermediate results (see image above) and performance metrics.

See ```examples/explanation.config``` for a config file with comments for every possible parameter.

## Setup Elastic Scattering

This project was developed and tested on Windows 10, it might not work on another OS. with Visual Studio 2019. You should get OpenCL lib/header files (version >2.0) from the [Intel SDK for OpenCL](https://software.intel.com/en-us/intel-opencl). Open the solution file and press F5 to build and run.

## Reproducing Graphs

Once setup, you can run any simulations you would like. To create the figures shown in the paper, we used settings and code which are documented in the ```reproduce_graphs``` folder. In particular, the SM2 folder contains the config file used for the computation as well as the output of the program.

The python code was executed with the following dependencies:
- Python 3.8.2
- NumPy 1.17.2
- SciPy 1.4.1
- Numba 0.49.0
- Matplotlib 3.5.1


