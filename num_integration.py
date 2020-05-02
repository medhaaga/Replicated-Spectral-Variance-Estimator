import numpy as np
import scipy
import scipy.integrate
from scipy.integrate import quad, dblquad

def density(A,C,x,y):
    return np.exp(-0.5*(A*(x*y)**2 + x**2 + y**2 - 2*C*x - 2*C*y))

def num_integration(A,C):
    return dblquad(lambda x,y: density(A,C,x,y), -np.inf, np.inf, lambda x: -np.inf, lambda x: np.inf)

def mean_fn(A,C,normalise_const,x,y):
    return x*np.exp(-0.5*(A*(x*y)**2 + x**2 + y**2 - 2*C*x - 2*C*y))/normalise_const

def expectation_x(A,C,normalise_const):
    return dblquad(lambda x,y: mean_fn(A,C,normalise_const,x,y), -np.inf, np.inf, lambda x: -np.inf, lambda x: np.inf)
