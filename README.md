# Purpose
3-PG2Py is a pure Python version of 3-PG2 (a modified version of 3-PG with respect to water balance prediction) to facilitate its extension and application to broader communities. Except the single plot simulation same as the original 3-PG2, 3-PG2Py includes two global sensitivity analyses algorithms, i.e., the variance-based sensitivity analysis method and Fourier amplitude sensitivity test, and the state-parameter estimation using ensemble Kalman filter algorithm. Additionally, an interface for spatial simulation was also implemented. 3-PG2Py is compatible with Python2.7+. With 3-PG2Py, the users could adapt the model easily to more diversified applications, especially computationally intensive ones.
# Usage
usage: python main.py -m <running mode> -a <end age the simulation> 
  
Arguments:
  
-m: running mode and could be one of ['r', 'spatial', 'VB', 'FAST', 'EnKF'], in which ‘r’ means single plot simulation, ‘spatial’ means the spatial simulation, ‘VB’ means the variance-based GSA, ‘FAST’ means the FAST GSA, and ‘EnKF’ means the state-parameter estimation using ensemble Kalman Filter.
  
-a: end age of the simulation.
  
Optional arguments:
  
-h: help information.
  
# Installation
3-PG2Py is compatible with Python 2.7+. 
  
# Author and contact
Xiaodong Song (xdsongy@gmail.com)
  
# Citation
Yu Song (2021). Introducing 3-PG2Py, an open-source forest growth model in Python. Environmental Modelling and Software (in review).
