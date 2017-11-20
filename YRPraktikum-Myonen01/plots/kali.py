#a = 0.0223839634048 ± 1.42384531816e-05
#b = -0.0141610316421 ± 0.00365692550549
import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
x= np.linspace(0,500)
k,t= np.genfromtxt('Calib.txt', unpack=True) #k kanäle, t in µs
#mhub = const.value('Bohr magneton') #das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt
# def mittel(x):              #the real mean()-ing of life
#     return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
# def relf(l,m):  #in Prozent
#     return (np.absolute(l-m)/l)*100
def f(x,a,b):
    return a*x+b
#Fit
params, covariance = curve_fit(f, k, t)
errors = np.sqrt(np.diag(covariance))
print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])

# Tabelle
# np.savetxt('tab.txt',np.column_stack([x,y]), delimiter=' & ',newline= r'\\'+'\n' )
#plt.subplot(1, 2, 1)
plt.plot(k, t,'kx', label='Messdaten')
plt.plot(x, f(x,*params),'r-',linewidth=2, label='Fit')
plt.ylabel(r'$t \:/\: \mu s$')
plt.xlabel(r'Kanal')
plt.legend(loc='best')
plt.savefig('plotKanal.pdf')
plt.clf()
#plt.subplot(1, 2, 2)
# plt.plot(x, y, label='Kurve')
# plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
# plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
# plt.legend(loc='best')

# in matplotlibrc leider (noch) nicht möglich
#plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
#plt.savefig('build/plot2.pdf')
