# OOFDTD
Parallel Object-Oriented Finite Difference Time Domain

This is an experiment with a FDTD code in a purely object oriented framework. This is done for code maintainability and extensibility. Obviously, performance is sacrificed for this cause but the code is parallelized with MPI to provide the ability to speedup the 

- Every Yee cell is an object. It is easy to derive new cell objects from the base class to add frequency dispersion, etc..
- Everything is query-able. ANY fields available at each cell (e.g. Ex, Hy, Susceptibility, Polarization...) can be probed and/or Fourier transformed for spectral content. 

The code in this project is built with CMake
- mkdir build
- cd build
- cmake ..
- make

NOTE: see oofdtd.pdf for a deeper explanation of the code and parallel performance. 
