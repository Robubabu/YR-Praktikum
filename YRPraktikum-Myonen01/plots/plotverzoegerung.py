import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)

t,N=np.genfromtxt('verzoegerung.txt', unpack=True)
errN=np.sqrt(N)
N2=N[5:26]

x = np.linspace(0, 10, 1000)
#mhub = const.value('Bohr magneton') #das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt
def mittel(x):              #the real mean()-ing of life
    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100

#Fit
#params , cov = curve_fit(f , x ,y )
#params = correlated_values(params, cov)
#for p in params:
#    print(p)

Nmeant=np.mean(N2) #207.09
Nmean=205 #fehlerbalken
#Tabelle

# np.savetxt('tab.txt',np.column_stack([x,y]), delimiter=' & ',newline= r'\\'+'\n' )
#plt.subplot(1, 2, 1)
plt.errorbar(t, N,yerr=errN, fmt='kx', label='Messdaten')
plt.xlabel(r'$T \:/\: s$')
plt.ylabel(r'Counts$ \:/\: 10s$')
plt.xlim(-18,25)
plt.plot((-6, 16), (Nmean,Nmean), 'r-', label='Plateau')
plt.plot((4, 4), (0, 300), 'k-')
plt.plot((-9, 17), ((Nmean/2),(Nmean/2)), 'r:', label='Halbwertsbreite')
plt.plot((-9, -9), (0, 300), 'k.-.')
plt.plot((17, 17), (0, 300), 'k.-.')
plt.legend(loc='best')
plt.savefig('plotVerzoegerung.pdf')
plt.clf()
#plt.subplot(1, 2, 2)
#plt.plot(x,y, label='Kurve')
#plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
#plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
#plt.legend(loc='best')

# in matplotlibrc leider (noch) nicht möglich
#plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
#plt.savefig('build/plot2.pdf')
