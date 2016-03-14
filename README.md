# OOFDTD
Object-Oriented Finite Difference Time Domain

This is an experiment with a FDTD code in a purely object oriented framework. This is done for code maintainability and extensibility. Obviously, performance is sacrificed for this cause but hopefully the sacrifice is minimized by providing a fully parallel implementation. 

- Everything is query-able. ANY fields available at each cell (e.g. Ex, Hy, Susceptibility, Polarization...) can be probed and/or Fourier transformed for spectral content. 
