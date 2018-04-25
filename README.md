### Code to accompany [A formalism for steering with local quantum measurements](https://arxiv.org/abs/1708.00756)  
#### A. B. Sainz, L. Aolita, M. Piani, M. J. Hoban, P. Skrzypczyk

This repository provides a small collection of code which reproduces the numerical results presented in   
[A formalism for steering with local quantum measurements](https://arxiv.org/abs/1708.00756).  
A. B. Sainz, L. Aolita, M. Piani, M. J. Hoban, P. Skrzypczyk
arXiv:1708.00756

The main file is a Jupyter notebook that presents all of the codes (written in MATLAB) and performs all of the calculations required to obtain
the numerical results presented in this work. The notebook can be most conveniently viewed 
[here](https://nbviewer.jupyter.org/github/paulskrzypczyk/o-formalism-steering/blob/master/post-quantum-steering.ipynb).

In order to run the Jupyter notebook, it is necessary to have installed
- [Jupyter notebook](http://jupyter.org/)
- MATLAB and the [API for Python](https://uk.mathworks.com/help/matlab/matlab-engine-for-python.html)
- The [MATLAB kernal](https://github.com/Calysto/matlab_kernel) for Jupyter.

The MATLAB codes require in addition:
- [CVX](http://cvxr.com/) - a Matlab-based convex modeling framework
- [QETLAB](http://www.qetlab.com/) - A MATLAB Toolbox for Quantum Entanglement

Everything has been tested on Matlab R2017a, and CVX 2.1 

The MATLAB code comprises the following:

  - [GenerateAssem.m](https://github.com/paulskrzypczyk/o-formalism-steering/blob/master/GenerateAssem.m): generates an assemblage
  given a shared quantum state and measurements for Alice and Bob
  - [BreuerMapOnAssem.m](https://github.com/paulskrzypczyk/o-formalism-steering/blob/master/BreuerMapOnAssem.m): Applies the
  Breuer map to the elements of an assemblage
  - [IsAQAssemblage.m](https://github.com/paulskrzypczyk/o-formalism-steering/blob/master/IsAQAssemblage.m): Determines whether
  an assemblage belongs to the set of 'almost-quantum' assemblages
  
