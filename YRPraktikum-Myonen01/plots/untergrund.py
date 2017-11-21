# 0.000280495500997 p
# 589.756657109 Nfehl
# 1.15186847092 Unter
# 26+/-19
import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.umath import *
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
x = np.linspace(0, 10, 1000)
mhub = const.value('Bohr magneton') #das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt
def mittel(x):              #the real mean()-ing of life
    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100
def poisson(a): #nur wenn k=1 aber passt hier
    return ((a)*unp.exp(-a))
Ufit=ufloat(0.848559707879, 0.221185434551)
startim=ufloat(2102553,np.sqrt(2102553))
stopim=ufloat(15711, np.sqrt(15711))
gesamtzeit=83930*(10**6) #µs
suchzeit=11.2 #µs
kanal=512
startimpulsrate= startim/gesamtzeit
P=poisson(startimpulsrate*suchzeit)
Nfehl=P*startim
U=Nfehl/kanal
fehler=relf(Ufit,U)
print(P)
print(Nfehl)
print(U)
print(fehler)
