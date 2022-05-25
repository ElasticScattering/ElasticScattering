# Elastic Scattering

![T3 MF0](https://user-images.githubusercontent.com/25907608/170357566-d476233d-febd-4996-9374-0275faa2c3e2.png)

Computes the lifetime of particles scattering around in a field of impurities.

The project can be run by giving it a config file, for example: `cl-es.exe small.config`. The config allows you to configure a field of impurities and the number of particles in them. This field can be used multiple times for different magnetic fields and temperatures, making it easy to batch a large number of simulations.

The results of a simulation are put inside a folder, with individual and combined results for all samples used. Furthermore, it includes visualizations of all intermediate results (see image above) and performance metrics

![image](https://user-images.githubusercontent.com/25907608/170358570-70f2730e-2f03-472f-802a-aaa32f88a4d1.png)

## Setup

This project was developed and tested on Windows 10, with Visual Studio 2019. You should get OpenCL lib/header files (version >2.0) from the [Intel SDK for OpenCL](https://software.intel.com/en-us/intel-opencl). Open the solution file and press F5 to build and run.
