# coding=utf-8
'''
Created on 2017/06/15

@author: snakagawa
'''

import numpy as np
from scipy import optimize
import time

__author__ = "Shintaro NAKAGAWA"
__email__ = "snaka@iis.u-tokyo.ac.jp"
__version__ = "1.0"

def CONTIN(tau, g1, N_gamma, range_gamma, alpha, verbose=False, reconst=False, full_result=False):
    '''
    A Python/numpy/scipy implementation of the CONTIN algorithm 
    for inverse Laplace transformation of DLS correlation functions. 
    
    Reference: 
    Scotti, A. et al J. Chem. Phys. 2015, 142, 234905. 
    doi: 10.1063/1.4921686
    
    @param tau: time axis of the input correlation function. 
    @param g1: input correlation function. 
    @param N_gamma: length of the solution vector. 
    @param range_gamma: lower and upper limit of the solution vector. 
    @param alpha: regularizer. 
    @param verbose: if True, message will be printed when the minimization is finished. 
    @param reconst: if True, a reconstructed g1 will be returned in addition to gammma and G(gamma). 
    @param full_result: if True, returns a scipy.optimize.OptimizeResult object. RuntimeError will not be raised. 
    
    @return: gamma, G(gamma)
    
    @raise RuntimeError: will be raised when the minimization does not converge. 
    '''
    assert len(tau) == len(g1)
    
    N_tau = len(tau)
    
    if verbose:
        tm = time.time()
        print "Preparing arrays...", 
    
    # gamma-axis of the solution vector
    gamma = np.logspace(np.log10(range_gamma[0]), 
                        np.log10(range_gamma[1]), 
                        N_gamma)
    # initial solution vector, G(gamma)
    x = np.zeros_like(gamma)
    
    # matrix that transforms the solution to the correlation function
    # i.e., g1 = A.x
    A = np.stack([np.exp(-gamma*t) for t in tau])
    
    if verbose:
        print "done. (ca. {0:f} sec). ".format(time.time() - tm)
    
    # define the function to minimize. 
    def V(u):
        # prevent any element of the solution vector from being negative. 
        ua = np.abs(u)
        # function to minimize. 
        U = np.power(np.linalg.norm(A.dot(ua) - g1), 2.0)
        # regularization term. 
        R = np.power(alpha * np.linalg.norm(np.diff(np.diff(ua))), 2.0)
        
        return U + R
    
    if verbose:
        print "Minimizing..."
        tm = time.time()

    # minimize V(x)
    # so far, the Powell algorithm is the fastest & most reliable method. 
    res = optimize.minimize(V, x, 
                            method="Powell", 
                            options={"disp": verbose, "maxiter": N_gamma*128})
    
    if verbose:
        print "Minimization done (ca. {0:f} sec). ".format(time.time() - tm)
    
    xa = np.abs(res.x)
    res.x = xa
    
    if full_result:
        return res
    elif not res.success:
        raise RuntimeError(res.message)
    elif reconst:
        return gamma, xa, A.dot(xa)
    else:
        return gamma, xa

    
if __name__ == '__main__':
    
    from matplotlib import pyplot
    
    tau, g1 = np.loadtxt("test.txt")
    
    gamma, G, g1_re = CONTIN(tau, g1, 
                             N_gamma=64, 
                             range_gamma=[1e-4, 1e4], 
                             alpha=0.1, 
                             verbose=True, 
                             reconst=True)
        
    pyplot.semilogx()
    
    axl = pyplot.gca()
    axr = pyplot.twinx()
    
    axl.set_yscale("log")
    axl.set_xlim(1e-4, 1e4)
    axl.set_ylim(1e-3, 1e0)
    axl.set_xlabel(ur"$\tau$ [ms]")
    axl.set_ylabel(ur"$g^{(1)}(\tau)$ [-]")
    
    axr.set_xlim(1e-4, 1e4)
    axr.set_ylim(0.0, 1.0)
    axr.set_ylabel(ur"$G(\tau)$ [-]")
    
    axl.plot(tau, g1,    ls="", marker="o", mec="k", mfc="w", label="Experimental")
    axl.plot(tau, g1_re, ls="-", marker="", color="b", label="Reconstructed")
    axr.plot(1.0/gamma, G/np.max(G), ls="-", marker="o", color="g", mec="g", mfc="w", label="Distribution")
    
    pyplot.show()

    