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
def g(x,m,n):
    return m*x+n

Ve=ufloat(10.9,0.9)

xlin=np.linspace(-10,150 , 1000)
#p=1mbar
p1=[1,2,3,4,5,6,7,8,9,10]
p1err=[0.2,0.2*2,0.2*3,0.2*4,0.2*5,0.2*6,0.2*7,0.2*8,0.2*9,0.2*10]
p1mf=unp.uarray(p1,p1err)
t10=[0,0,0]
t11=[10.33, 10.45, 9.99]
t12=[19.98, 18.53, 19.65]
t13=[34.64, 33.15, 30.51]
t14=[44.01, 44.24, 42.24]
t15=[57.81, 55.82, 54.51]
t16=[75.91, 75.57, 74.43]
t17=[89.78, 86.33, 87.44]
t18=[103.05, 104.70, 95.94]
t19=[126.90, 122.17, 133.25]
t10m=mittel(t10)
t11m=mittel(t11)
t12m=mittel(t12)
t13m=mittel(t13)
t14m=mittel(t14)
t15m=mittel(t15)
t16m=mittel(t16)
t17m=mittel(t17)
t18m=mittel(t18)
t19m=mittel(t19)
t1m=[t10m,t11m,t12m,t13m,t14m,t15m,t16m,t17m,t18m,t19m]
print(t1m)
#p=01mbar
p2=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
p2err=[0.02,0.2*0.2,0.2*0.3,0.2*0.4,0.2*0.5,0.2*0.6,0.2*0.7,0.2*0.8]
p2mf=unp.uarray(p2,p2err)
t20=[0,0,0]
t21=[ 10.97,  13.59,  14.97]
t22=[  27.9,  29.63,  29.97]
t23=[ 49.15,  49.91,  49.45]
t24=[ 63.49,  65.36,  65.27]
t25=[ 80.91,  83.10,  84.01]
t26=[104.43, 105.74, 105.24]
t27=[120.23, 124.02, 123.86]
t20m=mittel(t20)
t21m=mittel(t21)
t22m=mittel(t22)
t23m=mittel(t23)
t24m=mittel(t24)
t25m=mittel(t25)
t26m=mittel(t26)
t27m=mittel(t27)
t2m=[t20m,t21m,t22m,t23m,t24m,t25m,t26m,t27m]
print(t2m)
#p=04mbar
p3=[0.4,1  ,1.5,2  ,2.5,3  ,3.5,4  ]
p3err=[0.4*0.2, 0.2*1  ,0.2*1.5,0.2*2  ,0.2*2.5,0.2*3  ,0.2*3.5,0.2*4]
p3mf=unp.uarray(p3,p3err)
t30=[0,0,0]
t31=[20.18,  20.39,  19.48]
t32=[31.77,  33.25,  32.43]
t33=[50.87,  52.66,  51.41]
t34=[64.74,  65.97,  66.19]
t35=[80.20,  80.11,  79.38]
t36=[99.70,  99.94,  99.94]
t37=[127.1, 124.27, 127.43]
t30m=mittel(t30)
t31m=mittel(t31)
t32m=mittel(t32)
t33m=mittel(t33)
t34m=mittel(t34)
t35m=mittel(t35)
t36m=mittel(t36)
t37m=mittel(t37)
t3m=[t30m,t31m,t32m,t33m,t34m,t35m,t36m,t37m]
print(t3m)
#p=08mbar
p4=[0.8,2,3,4,5,6,7,8]
p4err=[0.2*0.8,0.2*2,0.2*3,0.2*4,0.2*5,0.2*6,0.2*7,0.2*8]
p4mf=unp.uarray(p4,p4err)
t40=[0,0,0]
t41=[ 17.83,  18.23,  17.7]
t42=[ 29.13,  29.56,  28.82]
t43=[ 51.04,  48.94,  48.42]
t44=[ 63.62,  63.57,  62.97]
t45=[ 79.68,  81.27,  79.22]
t46=[106.97, 105.83, 107.01]
t47=[123.16, 127.22, 123.07]
t40m=mittel(t40)
t41m=mittel(t41)
t42m=mittel(t42)
t43m=mittel(t43)
t44m=mittel(t44)
t45m=mittel(t45)
t46m=mittel(t46)
t47m=mittel(t47)
t4m=[t40m,t41m,t42m,t43m,t44m,t45m,t46m,t47m]
print(t1m)
params1, covariance1 = curve_fit(f=g, xdata=noms(t1m), ydata=noms(p1mf))
errors1 = np.sqrt(np.diag(covariance1))
#print('a2 =', params1[0], '±', errors1[0])
#print('b2 =', params1[1], '±', errors1[1])
params2, covariance2 = curve_fit(f=g, xdata=noms(t2m), ydata=noms(p2mf))
errors2 = np.sqrt(np.diag(covariance2))
#print('a2 =', params2[0], '±', errors2[0])
#print('b2 =', params2[1], '±', errors2[1])
params3, covariance3 = curve_fit(f=g, xdata=noms(t3m), ydata=noms(p3mf))
errors3 = np.sqrt(np.diag(covariance3))
#print('a3 =', params3[0], '±', errors3[0])
#print('b3 =', params3[1], '±', errors3[1])
params4, covariance4 = curve_fit(f=g, xdata=noms(t4m), ydata=noms(p4mf))
errors4 = np.sqrt(np.diag(covariance4))
#print('a4 =', params4[0], '±', errors4[0])
#print('b4 =', params4[1], '±', errors4[1])
a1 =ufloat(params1[0], errors1[0])
a2 =ufloat(params2[0], errors2[0])
a3 =ufloat(params3[0], errors3[0])
a4 =ufloat(params4[0], errors4[0])
#Saugleistungen
S1=(Ve/p1mf[0])*a1
S2=(Ve/p2mf[0])*a2
S3=(Ve/p3mf[0])*a3
S4=(Ve/p4mf[0])*a4
print(S1)
print(S2)
print(S3)
print(S4)
Stheo=1.1
Sexp=ufloat(0.92 , 0.07)
Slin1=ufloat(1.05, 0.08) #2-100
Slin2=ufloat(0.73, 0.07) #0,2-1
Slin3=ufloat(0.24, 0.04) #0,02-0,1

plt.errorbar(noms(t1m),noms(p1mf), xerr=stds(t1m), yerr=stds(p1mf),fmt='rx',linewidth=1, label='Messdaten')
plt.plot(xlin, g(xlin,*params1),'k-',linewidth=2, label='Fit')
plt.xlim(-10,135)
#plt.ylim(-0.001,0.004)
plt.grid(True)
plt.ylabel(r'$p\,/\,$mbar')
plt.xlabel(r'$t\,/\,$s')
plt.legend(loc='best')
plt.savefig('LeckrateDreh1.pdf')
plt.clf()
plt.errorbar(noms(t2m),noms(p2mf), xerr=stds(t2m), yerr=stds(p2mf),fmt='rx',linewidth=1, label='Messdaten')
plt.plot(xlin, g(xlin,*params2),'k-',linewidth=2, label='Fit')
plt.xlim(-10,125)
#plt.ylim(-0.001,0.004)
plt.grid(True)
plt.ylabel(r'$p\,/\,$mbar')
plt.xlabel(r'$t\,/\,$s')
plt.legend(loc='best')
plt.savefig('LeckrateDreh0_1.pdf')
plt.clf()
plt.errorbar(noms(t3m),noms(p3mf), xerr=stds(t3m), yerr=stds(p3mf),fmt='rx',linewidth=1, label='Messdaten')
plt.plot(xlin, g(xlin,*params3),'k-',linewidth=2, label='Fit')
plt.xlim(-10,130)
#plt.ylim(-0.001,0.004)
plt.grid(True)
plt.ylabel(r'$p\,/\,$mbar')
plt.xlabel(r'$t\,/\,$s')
plt.legend(loc='best')
plt.savefig('LeckrateDreh0_4.pdf')
plt.clf()
plt.errorbar(noms(t4m),noms(p4mf), xerr=stds(t4m), yerr=stds(p4mf),fmt='rx',linewidth=1, label='Messdaten')
plt.plot(xlin, g(xlin,*params4),'k-',linewidth=2, label='Fit')
plt.xlim(-10,130)
#plt.ylim(-0.001,0.004)
plt.grid(True)
plt.ylabel(r'$p\,/\,$mbar')
plt.xlabel(r'$t\,/\,$s')
plt.legend(loc='best')
plt.savefig('LeckrateDreh0_8.pdf')
plt.clf()

#plt.plot((0,250), (Stheo,Stheo), 'k-', label=r'$S_{theo}$')
#plt.plot((0,250), (noms(Sexp),noms(Sexp)), 'r-', label=r'$S_{g}$')
plt.errorbar(50,noms(Stheo), xerr=50, yerr=0,fmt='kx',linewidth=1, label=r'$S_{theo}$')
plt.errorbar(50,noms(Sexp), xerr=50, yerr=stds(Sexp),fmt='rx',linewidth=1, label=r'$S_{g}$')
##plt.errorbar(125,noms(Sexp), xerr=0, yerr=stds(Sexp),fmt='rx',linewidth=1, label=r'$S_{g}$')
#plt.plot((2,100), (noms(Slin1),noms(Slin1)), 'b-', label=r'$S_{lin1}$')
plt.errorbar(49,noms(Slin1), xerr=47, yerr=stds(Slin1),fmt='bx',linewidth=1, label=r'$S_{lin1}$')
#plt.plot((0.2,1), (noms(Slin2),noms(Slin2)), 'g-', label=r'$S_{lin2}$')
plt.errorbar(0.6,noms(Slin2), xerr=0.4, yerr=stds(Slin2),fmt='gx',linewidth=1, label=r'$S_{lin2}$')
#plt.plot((0.02,0.1), (noms(Slin3),noms(Slin3)), 'm-', label=r'$S_{lin3}$')
plt.errorbar(0.06,noms(Slin3), xerr=0.4, yerr=stds(Slin3),fmt='mx',linewidth=1, label=r'$S_{lin3}$')

#plt.plot((0.1,0.8), (noms(S2),noms(S2)), 'y-', label=r'$S_{leck1}$')
plt.errorbar(0.45,noms(S2), xerr=0.35, yerr=stds(S2),fmt='yx',linewidth=1, label=r'$S_{leck1}$')
#plt.plot((1,10), (noms(S1),noms(S1)), 'c-', label=r'$S_{leck2}$')
plt.errorbar(5.5,noms(S1), xerr=4.5, yerr=stds(S1),fmt='cx',linewidth=1, label=r'$S_{leck2}$')
#plt.plot((0.4,4), (noms(S3),noms(S3)), 'b.', label=r'$S_{leck3}$')
plt.errorbar(2.2,noms(S3), xerr=1.8, yerr=stds(S3),fmt='bx',linewidth=1, label=r'$S_{leck3}$')
#plt.plot((0.8,8), (noms(S4),noms(S4)), 'r.', label=r'$S_{leck4}$')
plt.errorbar(4.4,noms(S4), xerr=3.6, yerr=stds(S4),fmt='rx',linewidth=1, label=r'$S_{leck4}$')
plt.xlim(-2,100)
plt.grid(True)
plt.xlabel(r'$p\,/\,$mbar')
plt.ylabel(r'$S\,/\,\frac{l}{s}$')
plt.legend(loc=4)
plt.savefig('SaugverDreh.pdf')
plt.clf()
