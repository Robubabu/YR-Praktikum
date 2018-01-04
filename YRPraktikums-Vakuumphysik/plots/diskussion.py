import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)

def mittel(x):              #the real mean()-ing of life
    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100
def loga(p,pe,po):
    a=unp.log(((p-pe)/(po-pe)))
    return a
def f(x,a,b):
    return (a*np.exp((-x)*b))+0.012
def g(x,m,n):
    return m*x+n

Stdreh=ufloat(1.1,0)
Stturbo=ufloat(77,0)
Sd=ufloat(1.05,0.08)
St=ufloat(25,3)

fd=relf(Stdreh,Sd)
ft=relf(Stturbo,St)
print(fd)
print(ft)
