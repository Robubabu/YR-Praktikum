import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
x = np.linspace(-0.03,0.03 , 1000)
mhub = const.value('Bohr magneton') #das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt
def mittel(x):              #the real mean()-ing of life
    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100
def fite(x,o,A):
	return A*np.exp(-0.5*(1/o**2)*x**2) #funktioniert nicht :(
def fitf(x,o,l):
	return (l/(np.sqrt(2*np.pi)*o))*np.exp(-0.5*(1/o**2)*x**2) #funktioniert
z , B = np.genfromtxt('BFeld.txt', unpack = True)
z=z-127 #Zentrum des Magfeldes 
z*=10**-3
B*=10**-3

#Fit
params , cov = curve_fit(fitf, z ,B)
params = correlated_values(params, cov)
o = params[0]
l = params[1]
#params2 , cov2 = curve_fit(fite, z ,B)
#params2 = correlated_values(params2, cov2)
#a = params2[0]
#b = params2[1]
print((l/(np.sqrt(2*np.pi)*o)))
print(o,l)
#Tabelle
np.savetxt('BFeldtab.txt',np.column_stack([B,z]), delimiter=' & ',newline= r'\\'+'\n' )

#plt.subplot(1, 2, 1)
plt.plot(z, B,'ro', label='Mag. Feldstärke')
plt.plot(x, fitf(x,noms(o),noms(l)), 'b-', label='Gaußfunktion')
#plt.plot(x, fite(x,noms(a),noms(b)), 'c-', label='Gaußfunktion2')
plt.xlabel(r'$z \:/\: m$')
plt.ylabel(r'$B \:/\: T$')
plt.xlim(-0.03,0.03)
plt.grid()
plt.legend(loc='best')
#plt.show()
plt.savefig('BFeld.pdf')

#plt.clf()

#plt.subplot(1, 2, 2)
#plt.plot(x, y, label='Kurve')
#plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
#plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
#plt.legend(loc='best')

# in matplotlibrc leider (noch) nicht möglich
#plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
#plt.savefig('build/plot2.pdf')
