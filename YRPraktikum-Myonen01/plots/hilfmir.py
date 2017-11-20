import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
x = np.linspace(1, 512, 512)
y = np.genfromtxt('Spektrum.Spe', unpack=True)
np.savetxt('Spektrum2.txt',np.column_stack([x,y]))
