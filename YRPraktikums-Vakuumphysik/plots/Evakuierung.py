import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from uncertainties.umath import *
from astropy.io import ascii

def mittel(x):              #the real mean()-ing of life
    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100
def loga(p,pe,po):
    a= unp.log(((p-pe)/(po-pe)))
    return a
def f(x,a,b):
    return (a*np.exp(b*(-x)))+(1*10**(-5))
def g(x,m,n):
    return m*x+n
Ve=ufloat(11,0.8)
p,t1r,t2r,t3r,t4r,t5r = np.genfromtxt('Turbo_p.txt', unpack = True)
t0=[0.99, 1.18, 1.05, 0.99, 0.92]
t1=[ 2.42, 2.70, 2.55, 2.51, 2.52]
t2=[ 2.83, 3.33, 3.11, 3.25, 3.15]
t3=[ 4.70, 4.93, 4.94, 4.84, 4.83]
t4=[ 5.23, 5.55, 5.44, 5.49, 5.35]
t5=[ 6.62, 6.97, 6.85, 6.91, 6.59]
pE=ufloat(10**(-5), 0.5*10**(-5))
p0=ufloat(7.8*10**(-3),7.8*10**(-4))
perr=p*0.1
pmf=unp.uarray(p,perr)
Logwerte=loga(pmf,pE,p0)
t0m=mittel(t0)
t1m=mittel(t1)
t2m=mittel(t2)
t3m=mittel(t3)
t4m=mittel(t4)
t5m=mittel(t5)
tm=[t0m,t1m,t2m,t3m,t4m,t5m]
xlin=np.linspace(0.1,10 , 1000)
xlin2=np.linspace(0.1,4.5 , 1000)
xlin3=np.linspace(3.5,7.5 , 1000)
xlin4=np.linspace(4,7.5 , 1000)
#expfit
params, covariance = curve_fit(f=f, xdata=noms(tm), ydata=p)
errors = np.sqrt(np.diag(covariance))
print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])
para=ufloat(params[1], errors[1])
Sgemittelt= para*Ve
print(Sgemittelt)
#geradenfit
params2, covariance2 = curve_fit(f=g, xdata=noms(tm[0:3]), ydata=noms(Logwerte[0:3]))
errors2 = np.sqrt(np.diag(covariance2))
print('a2 =', params2[0], '±', errors2[0])
print('b2 =', params2[1], '±', errors2[1])
params3, covariance3 = curve_fit(f=g, xdata=noms(tm[3:6]), ydata=noms(Logwerte[3:6]))
errors3 = np.sqrt(np.diag(covariance3))
print('a3 =', params3[0], '±', errors3[0])
print('b3 =', params3[1], '±', errors3[1])
# params4, covariance4 = curve_fit(f=g, xdata=noms(tm[4:6]), ydata=noms(Logwerte[4:6]))
# errors4 = np.sqrt(np.diag(covariance4))
# print('a4 =', params4[0], '±', errors4[0])
# print('b4 =', params4[1], '±', errors4[1])
para2=ufloat(params2[0], errors2[0])
S1= para2*Ve
print(S1)
para3=ufloat(params3[0], errors3[0])
S2= para3*Ve
print(S2)
# para4=ufloat(params4[0], errors4[0])
# S3= para4*Ve
# print(S3)
#
print(p[0:3])
print(p[3:6])
plt.errorbar(noms(tm),noms(pmf), xerr=stds(tm), yerr=stds(pmf),fmt='rx',linewidth=1, label='Messdaten')
plt.plot(xlin, f(xlin,*params),'k-',linewidth=2, label='Fit')
plt.plot((0,7), (noms(pE),noms(pE)), 'r-', label=r'$p_{E}$')
plt.xlim(0.5,7)
plt.ylim(-0.001,0.004)
plt.grid(True)
plt.ylabel(r'$p\,/\,$mbar')
plt.xlabel(r'$t\,/\,$s')
plt.legend(loc='best')
plt.axes([0.55, 0.45, 0.3, 0.2])
ax = plt.gca()
ax.errorbar(noms(tm[4:12]),noms(pmf[4:12]), xerr=stds(tm[4:12]), yerr=stds(pmf[4:12]),fmt='rx',linewidth=1, label='Messdaten')
ax.plot(xlin, f(xlin,*params),'k-',linewidth=2, label='Fit')
ax.plot((0,10), (noms(pE),noms(pE)), 'r-', label=r'$p_{E}$')
plt.grid(True)
plt.xlim(5,8)
plt.ylim(-0.00001,0.00003)
plt.savefig('EvakuierungTurboExp.pdf')
plt.clf()
plt.errorbar(noms(tm),noms(Logwerte), xerr=stds(tm), yerr=stds(Logwerte),fmt='rx',linewidth=1, label='Messdaten')
plt.plot(xlin2, g(xlin2,*params2),'r-',linewidth=2, label='Druckbereich 1')
plt.plot(xlin3, g(xlin3,*params3),'g-',linewidth=2, label='Druckbereich 2')
#plt.plot(xlin4, g(xlin4,*params4),'b-',linewidth=2, label='Druckbereich 3')
plt.xlim(0,7)
plt.ylim(-8,0)
plt.grid(True)
plt.ylabel(r'$ln\left(\frac{p(t)-p_{E}}{p_0-p_{E}}\right)$')
plt.xlabel(r'$t\,/\,$s')
plt.legend(loc='best')
plt.savefig('EvakuierungTurbolin.pdf')

# print(t0)
# print(t1)
# print(t2)
# print(t3)
# print(t4)
# print(t5)
# print(tm)
# [0.002+/-0.0002 0.0004+/-4e-05 0.0002+/-2e-05 6e-05+/-6e-06
#  4e-05+/-4.000000000000001e-06 2e-05+/-2.0000000000000003e-06]
# [-1.361348464830098+/-0.1414659206383651
#  -2.9727893824406344+/-0.14161273423013987
#  -3.6684459746060054+/-0.14180731484845902
#  -4.884213355424778+/-0.14288161643133884
#  -5.298189153200852+/-0.14382171436697727
#  -6.017311820164057+/-0.14755354160362774]
# [0.99, 1.18, 1.05, 0.99, 0.92]
# [2.42, 2.7, 2.55, 2.51, 2.52]
# [2.83, 3.33, 3.11, 3.25, 3.15]
# [4.7, 4.93, 4.94, 4.84, 4.83]
# [5.23, 5.55, 5.44, 5.49, 5.35]
# [6.62, 6.97, 6.85, 6.91, 6.59]
# [1.026+/-0.04365775990588613, 2.54+/-0.045497252664309346, 3.134+/-0.08518215775618741, 4.848000000000001+/-0.043289721643826704, 5.412000000000001+/-0.056071383075504695, 6.787999999999999+/-0.07722693830523127]
