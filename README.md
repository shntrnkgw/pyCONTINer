# pyCONTINer
A Python/numpy/scipy implementation of the famous inverse Laplace transformation algorithm CONTIN

## Overview
This is a re-implementation of the famous inverse Laplace transformation (ILT) 
algorithm _CONTIN_ using Python 2.7, numpy, and scipy. 

The original _CONTIN_ program, written in FORTRAN, was developed by 
[Stephen W. Provencher](http://s-provencher.com/) 
and have been used by countless scientists in various fields 
thanks to its excellent reliability and performance. 
However, although _CONTIN_ is regarded as a _de facto_ standard for ILT 
over the decades,  
the difficulty in understanding the original FORTRAN source code 
makes it a black box, which I, as a scientist, hate more than anything. 

That is why I decided to write my implementation of _CONTIN_ algorithm that is 
(hopefully) comprehensive to the users. 

In developing this program, I intensively referenced a paper 
by Scotti and coworkers (ref 1), 
in which the essence of the _CONTIN_ algorithm is explained 
in an extremely comprehensive manner. 

## Requirement
- Python 2.7
- numpy >= 1.13.0
- scipy >= 0.19.0

## Usage
1. Download the file `contin.py`. 
2. There you go. Use the method `contin.CONTIN()` as you like. 

Suppose that you have numpy arrays of the correlation function `g1` 
and corresponding delay time `tau`. 
You can then obtain the numpy arrays of the decay time `gamma` 
and its population `x` by doing something like this;  
```python
import contin

gamma, x = contin.CONTIN(tau, g1, N_gamma=64, range_gamma=[1e-4, 1e4], alpha=1.0)
```
`N_gamma` is the length of the output decay rate array. 
`range_gamma` specifies the lower and upper limit of the output decay rate. 
`alpha` is the parameter called _regularizer_ (see ref 1 for detail). 

## Caution & disclaimer
- This is my personal project and therefore has nothing to do with any institution. 
- The code is not thoroughly tested. Use it at your own risk.  

## How it works


## ToDo
- make it more object-oriented
- compare result with that from Provencher's CONTIN

## References
1. Scotti, A. et al "The CONTIN algorithm and its application to determine the size distribution of microgel suspensions", J. Chem. Phys. **2015**, _142_, 234905. doi: [10.1063/1.4921686](http://dx.doi.org/10.1063/1.4921686)

2. [pyCONTIN](https://github.com/kanhua/pyCONTIN) - a Python wrapper to the original Provencher's _CONTIN_