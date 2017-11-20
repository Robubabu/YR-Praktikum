import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
# x = np.linspace(0, 10, 1000)
# mhub = const.value('Bohr magneton') #das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt
# def mittel(x):              #the real mean()-ing of life
    # return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
# def relf(l,m):  #in Prozent
    #return (np.absolute(l-m)/l)*100

def poisson(k,a): #nur wenn k=1 aber passt hier
    return ((a**k)*np.exp(-a))

def f(t,a,U,N):
    return N*np.exp(-a*t)+U


x,y = np.genfromtxt('Spektrum2.txt', unpack=True)
erry=np.sqrt(y)

#Fit
# params, covariance = curve_fit(f, x, y)
# errors = np.sqrt(np.diag(covariance))
# print('a =', params[0], '±', errors[0])
# print('b =', params[1], '±', errors[1])


#Tabelle
# np.savetxt('tab.txt',np.column_stack([x,y]), delimiter=' & ',newline= r'\\'+'\n' )
#plt.subplot(1, 2, 1)
plt.errorbar(x, y, yerr=erry,fmt='rx', label='Messdaten')
#plt.plot(x, y, 'kx',label='Messwerte')
#plt.plot(x, f(x,*params),'r-',linewidth=2, label='Fit')
plt.xlabel(r'Kanal')
plt.ylabel(r'Counts')
plt.legend(loc='best')
plt.savefig('KanalCounts.pdf')
plt.clf()
#plt.subplot(1, 2, 2)
#plt.plot(x, y, label='Kurve')
#plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
#plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
#plt.legend(loc='best')

# in matplotlibrc leider (noch) nicht möglich
#plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
#plt.savefig('build/plot2.pdf')
