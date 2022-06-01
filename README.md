# Elastic Scattering

![T3 MF0](https://user-images.githubusercontent.com/25907608/170357566-d476233d-febd-4996-9374-0275faa2c3e2.png)

Program to computes the lifetime of particles scattering around in a field of impurities.

The project can be run by giving it a config file, for example: `cl-es.exe examples/small.config`. The config allows you to tune various parameters of the simulation, such as the magnetic field strength and temperature, and the scale of the simulation. 

The results of a simulation are put inside a folder, with optional visualizations of intermediate results (see image above) and performance metrics.

See ```examples/explanation.config``` for a config file with comments for every possible parameter.

## Setup

This project was developed and tested on Windows 10, it might not work on another OS. with Visual Studio 2019. You should get OpenCL lib/header files (version >2.0) from the [Intel SDK for OpenCL](https://software.intel.com/en-us/intel-opencl). Open the solution file and press F5 to build and run.

## Reproduction

Once setup, you can run any simulations you would like. To create the figures shown in the paper, we used settings and code which are documented in the ```reproduce_graphs``` folder. This folder also contains the code used to generate the graphs in the main text.
