import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)

d=ufloat(12,0.05)
l=ufloat(61,1)
Vzyeigen=((d/2)**2)*l*np.pi

ta =ufloat(9.5,0.8)
lS=ufloat(0.8,0.1)
kS=ufloat(0.087,0.011)
kT=ufloat(0.013,0.002)
gT=ufloat(0.25,0.01)
kK=ufloat(0.016 , 0.002)
gK=ufloat(0.177 , 0.09)
V=ufloat(0.0069 , 0.0001)
Q=ufloat(0.067 , 0.004)
oKu=ufloat(0.015 , 0.002)
zKu=ufloat(0.005 , 0.001)
V1=ufloat(0.044 , 0.004)
V4=ufloat(0.025,0.005)
VgesamtET=ta+gT+gK+Q+2*zKu+V1
VgesamtLT=ta+zKu+0.5*V1+gT+gK
VgesamtED=ta+lS+kS+kT+gT+kK+gK+V+oKu+2*zKu+0.5*V1+V4
VgesamtLD=ta+lS+kS+kT+gT+kK+gK+V+2*oKu+zKu+0.5*V1+0.5*V4
print ('3 significant digits on the uncertainty: {:.3u}'.format(VgesamtET))
print ('3 significant digits on the uncertainty: {:.3u}'.format(VgesamtLT))
print ('3 significant digits on the uncertainty: {:.3u}'.format(VgesamtED))
print ('3 significant digits on the uncertainty: {:.3u}'.format(VgesamtLD))
