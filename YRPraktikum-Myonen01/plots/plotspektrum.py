# a = 0.495128318768 ± 0.00891667974057
# U = 0.848559707879 ± 0.221185434551
# N = 163.07150846 ± 2.8239401661
# 0.495+/-0.009
# 2.02+/-0.04s
#fehler 8.1+/-1.7
import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)

xlin = np.linspace(0, 12, 1000)
# mhub = const.value('Bohr magneton') #das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt
def mittel(x):              #the real mean()-ing of life
     return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100

tautheo=ufloat(2.1969811,0.0000022)

katot1o=(0.0223839634048)
katot2o=(-0.0141610316421)
def poisson(k,a): #nur wenn k=1 aber passt hier
    return ((a**k)*np.exp(-a))

def f(t,a,U,N):
    return (N*np.exp(-a*t))+U

def g(x,a,b):
    return a*x+b

def katot(x): # kanal to time
    return katot1o*x+katot2o

x,y = np.genfromtxt('Spektrumfusch.txt', unpack=True)
erry=np.sqrt(y)
ylin=np.log(y)
errylin=np.log(erry)
x2= katot(x)
#Fit
params, covariance = curve_fit(f, x2, y,sigma=erry)
errors = np.sqrt(np.diag(covariance))
print('a =', params[0], '±', errors[0])
print('U =', params[1], '±', errors[1])
print('N =', params[2], '±', errors[2])

params2, covariance2 = curve_fit(g, x2, ylin)
errors2 = np.sqrt(np.diag(covariance2))
print('a =', params2[0], '±', errors2[0])
print('b =', params2[1], '±', errors2[1])

#lam= ufloat(params[0], errors[0])
#Zeitbla= 1/lam
# print(lam)
# print(Zeitbla)
# fehler=relf(tautheo,Zeitbla)
# print(fehler)
#Tabelle
# np.savetxt('tab.txt',np.column_stack([x,y]), delimiter=' & ',newline= r'\\'+'\n' )

# plt.errorbar(x, y, yerr=erry,fmt='rx', label='Messdaten')
# plt.xlabel(r'Kanal')
# plt.ylabel(r'Counts')
# plt.legend(loc='best')
# plt.savefig('KanalCounts.pdf')
# plt.clf()
# plt.errorbar(x2, y, yerr=erry,fmt='rx',linewidth=1, label='Messdaten')
# plt.plot(xlin, f(xlin,*params),'k-',linewidth=2, label='Fit')
# plt.xlabel(r'$t\:/\:\upmu s$')
# plt.ylabel(r'Counts')
# plt.legend(loc='best')
# plt.savefig('ZeitCounts.pdf')
# plt.clf()

plt.errorbar(x2, y, yerr=(erry),fmt='rx',linewidth=1, label='Messdaten')
plt.plot(xlin,f(xlin,*params),'k-',linewidth=2, label='Fit')
plt.xlabel(r'$t\:/\:\mu s$')
plt.yscale('log')
plt.ylabel(r'Counts')
plt.legend(loc='best')
plt.savefig('ZeitCountslinear.pdf')
plt.clf()
#plt.plot(x, y, 'kx',label='Messwerte')
#
