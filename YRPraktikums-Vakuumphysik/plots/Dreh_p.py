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

Ve=ufloat(10.9,0.9)

tr1,tr2,tr3,tr4,tr5= np.genfromtxt('Dreh_p_nurt.txt', unpack=True)
p0=ufloat(1000,200)
pE=ufloat(1.2*(10**(-2)),0.2*(1.2*(10**(-2))))
p=[100, 40 ,20 ,10 ,4 ,2,1 ,0.4 ,0.2,0.1 ,0.04 ,0.02]
fp=[100*0.2, 40*0.2 ,20*0.2 ,10*0.2 ,4*0.2 ,2*0.2,1*0.2 ,0.4*0.2 ,0.2*0.2,0.1*0.2 ,0.04*0.2 ,0.02*0.2]
pmf=unp.uarray(p,fp)
t1=[ 9.57, 11.69, 12.21, 15.35, 13.91]
t2=[22.16, 23.58, 23.64, 28.5 ,23.91]
t3=[27.84, 30.94, 31.25, 36.1 ,31.06]
t4=[35.27, 38.08, 38.22, 43.11, 38.52]
t5=[44.03, 47.02, 47.05, 51.73, 47.18]
t6=[50.45, 53.3 ,53.54 ,58.66 ,54.23]
t7=[57.57, 60.86, 60.74, 65.82, 60.86]
t8=[69.94, 73.72, 73.69, 79.68, 73.72]
t9=[82.48, 85.89, 86.01, 91.06, 85.95]
t10=[93.18, 96.60, 96.79, 102.09, 100.52]
t11=[132.1, 134.55, 133.69, 140.88, 133.88]
t12=[199.5, 204.34, 202.08, 219.14, 202.06]
t1m=mittel(t1)
t2m=mittel(t2)
t3m=mittel(t3)
t4m=mittel(t4)
t5m=mittel(t5)
t6m=mittel(t6)
t7m=mittel(t7)
t8m=mittel(t8)
t9m=mittel(t9)
t10m=mittel(t10)
t11m=mittel(t11)
t12m=mittel(t12)
tm=[t1m,t2m,t3m,t4m,t5m,t6m,t7m,t8m,t9m,t10m,t11m,t12m]
Logwerte=loga(pmf,pE,p0)
#print(loga(p0,pE,p0))???
xlin=np.linspace(5,250 , 1000)
xlin2=np.linspace(5,100 , 1000)
xlin3=np.linspace(50,140 , 1000)
xlin4=np.linspace(90,250 , 1000)
params, covariance = curve_fit(f=f, xdata=noms(tm), ydata=p)
errors = np.sqrt(np.diag(covariance))
# print('a =', params[0], '±', errors[0])
# print('b =', params[1], '±', errors[1])
para=ufloat(params[1], errors[1])
Sgemittelt= para*Ve
#print('sg =', Sgemittelt)
#fitsgeradeneva
params2, covariance2 = curve_fit(f=g, xdata=noms(tm[0:6]), ydata=noms(Logwerte[0:6]))
errors2 = np.sqrt(np.diag(covariance2))
#print('a2 =', params2[0], '±', errors2[0])
#print('b2 =', params2[1], '±', errors2[1])
params3, covariance3 = curve_fit(f=g, xdata=noms(tm[6:9]), ydata=noms(Logwerte[6:9]))
errors3 = np.sqrt(np.diag(covariance3))
#print('a3 =', params3[0], '±', errors3[0])
#print('b3 =', params3[1], '±', errors3[1])
params4, covariance4 = curve_fit(f=g, xdata=noms(tm[9:12]), ydata=noms(Logwerte[9:12]))
errors4 = np.sqrt(np.diag(covariance4))
#print('a4 =', params4[0], '±', errors4[0])
#print('b4 =', params4[1], '±', errors4[1])
para2=ufloat(params2[0], errors2[0])
S1= para2*Ve
#print('S1=',S1)
para3=ufloat(params3[0], errors3[0])
S2= para3*Ve
#print('S2=',S2)
para4=ufloat(params4[0], errors4[0])
sneu=ufloat(0.043,0.007)
S3= sneu*Ve
#print('S4=',S3)
# plt.axes([0.1, 0.1, 0.8, 0.8])
# ax = plt.gca()
# ax.errorbar(noms(tm),noms(pmf), xerr=stds(tm), yerr=stds(pmf),fmt='rx',linewidth=1, label='Messdaten')
# ax.plot(xlin, f(xlin,*params),'k-',linewidth=2, label='Fit')
# ax.plot((5,220), (noms(pE),noms(pE)), 'r-', label=r'$p_{E}$')
# plt.xlim(5,220)
# plt.ylim(-10,120)
# plt.grid(True)
# plt.ylabel(r'$p\,/\,$mbar')
# plt.xlabel(r'$t\,/\,$s')
# plt.legend(loc='best')
# plt.axes([0.55, 0.45, 0.3, 0.2])
# ax = plt.gca()
# ax.errorbar(noms(tm[6:12]),noms(pmf[6:12]), xerr=stds(tm[6:12]), yerr=stds(pmf[6:12]),fmt='rx',linewidth=1, label='Messdaten')
# ax.plot(xlin, f(xlin,*params),'k-',linewidth=2, label='Fit')
# ax.plot((0,250), (noms(pE),noms(pE)), 'r-', label=r'$p_{E}$')
# plt.grid(True)
# plt.xlim(100,250)
# plt.ylim(-0.05,0.05)
# plt.savefig('EvakuierungDrehExp.pdf')
# plt.clf()
# plt.errorbar(noms(tm),noms(Logwerte), xerr=stds(tm), yerr=stds(Logwerte),fmt='rx',linewidth=1, label='Messdaten')
# plt.plot(xlin2, g(xlin2,*params2),'r-',linewidth=2, label='Druckbereich1')
# plt.plot(xlin3, g(xlin3,*params3),'g-',linewidth=2, label='Druckbereich2')
# plt.plot(xlin4, g(xlin4,*params4),'b-',linewidth=2, label='Druckbereich3')
# plt.xlim(5,220)
# plt.ylim(-20,0)
# plt.grid(True)
# plt.ylabel(r'$ln\left(\frac{p(t)-p_{E}}{p_0-p_{E}}\right)$')
# plt.xlabel(r'$t\,/\,$s')
# plt.legend(loc='best')
# plt.savefig('EvakuierungDrehlin.pdf')
# plt.clf()
