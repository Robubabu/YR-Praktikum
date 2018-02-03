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

Ve=ufloat(11,0.8)
xlin=np.linspace(-10,100 , 1000)
#p2_4
p1=[2*10**(-4),6*10**(-4),8*10**(-4),2*10**(-3),4*10**(-3),6*10**(-3),8*10**(-3)]
p1err=[2*10**(-5),6*10**(-5),8*10**(-5),2*10**(-4),4*10**(-4),6*10**(-4),8*10**(-4)]
p1mf=unp.uarray(p1,p1err)
t10=[0,0,0]
t11=[0.92, 1.05 ,1.07 ,1.06, 1.05]
t12=[1.55, 1.77 ,1.51 ,1.75, 1.48]
t13=[4.43, 5.83 ,4.26 ,5.80, 4.45]
t14=[8.81, 11.08, 8.50, 10.97, 8.56]
t15=[12.52, 15.83, 12.19, 15.59, 12.26]
t16=[15.96, 20.19, 15.45, 19.76,15.2]
t10m=mittel(t10)
t11m=mittel(t11)
t12m=mittel(t12)
t13m=mittel(t13)
t14m=mittel(t14)
t15m=mittel(t15)
t16m=mittel(t16)
t1m=[t10m,t11m,t12m,t13m,t14m,t15m,t16m]
#print(p1mf)
#print(t1m)
#1_4
p2=[1*10**(-4),4*10**(-4) ,6*10**(-4) ,8*10**(-4) ,2*10**(-3) ,4*10**(-3) ,6*10**(-3) ,8*10**(-3)]
p2err=[1*10**(-5),0.1*4*10**(-4) ,0.1*6*10**(-4) ,0.1*8*10**(-3) ,0.1*2*10**(-3) ,0.1*4*10**(-3) ,0.1*6*10**(-3) ,0.1*8*10**(-3)]
p2mf=unp.uarray(p2,p2err)
t20=[0,0,0]
t21=[1.64, 1.7 ,1.32]
t22=[2.74, 2.83, 2.46]
t23=[3.89, 3.76, 3.48]
t24=[9.85, 9.72, 9.32]
t25=[18.40, 18.13, 17.66]
t26=[25.21, 25.21, 24.86]
t27=[32.16, 31.96, 31.68]
t20m=mittel(t20)
t21m=mittel(t21)
t22m=mittel(t22)
t23m=mittel(t23)
t24m=mittel(t24)
t25m=mittel(t25)
t26m=mittel(t26)
t27m=mittel(t27)
t2m=[t20m,t21m,t22m,t23m,t24m,t25m,t26m,t27m]
#print(t2m)
#p3_5
p3=[3*10**(-5),8*10**(-5),2*10**(-4),4*10**(-4),6*10**(-4),8*10**(-4),2*10**(-3)]
p3err=[3*10**(-6),0.8*10**(-5),0.2*10**(-4),0.4*10**(-4),0.6*10**(-4),0.8*10**(-4),0.2*10**(-3)]
p3mf=unp.uarray(p3,p3err)
t30=[0,0,0]
t31=[1.77, 0.6, 0.7 ,0.98]
t32=[5.26, 4.0, 3.48, 4.58]
t33=[11.68, 9.6, 8.23 ,10.73]
t34=[17.56, 14.78, 12.79, 17.01]
t35=[22.93, 20.24, 17.01, 23.12]
t36=[58.08, 60.2 ,41.35 ,68.07]
t30m=mittel(t30)
t31m=mittel(t31)
t32m=mittel(t32)
t33m=mittel(t33)
t34m=mittel(t34)
t35m=mittel(t35)
t36m=mittel(t36)
t3m=[t30m,t31m,t32m,t33m,t34m,t35m,t36m]
# hilf=unp.uarray(noms(t3m),noms(p3mf))
# print(hilf)
#print(t3m)
#p8_5
p4=[8*10**(-5),4*10**(-4),6*10**(-4),8*10**(-4),2*10**(-3),4*10**(-3),6*10**(-3),8*10**(-3)]
p4err=[8*10**(-6),0.4*10**(-4),0.6*10**(-4),0.8*10**(-4),0.2*10**(-3),0.4*10**(-3),0.6*10**(-3),0.8*10**(-3)]
p4mf=unp.uarray(p4,p4err)
t40=[0,0,0]
t41=[2.29 ,  2.69, 2.40]
t42=[3.89 ,  4.35, 4.34]
t43=[5.28 ,  5.74, 6.35]
t44=[13.07, 13.95, 18.15]
t45=[24.26, 25.19, 35.31]
t46=[33.88, 35.25, 49.92]
t47=[42.54, 44.08, 60.02]
t40m=mittel(t40)
t41m=mittel(t41)
t42m=mittel(t42)
t43m=mittel(t43)
t44m=mittel(t44)
t45m=mittel(t45)
t46m=mittel(t46)
t47m=mittel(t47)
t4m=[t40m,t41m,t42m,t43m,t44m,t45m,t46m,t47m]
#print(t4m)
params1, covariance1 = curve_fit(f=g, xdata=noms(t1m), ydata=noms(p1mf))
errors1 = np.sqrt(np.diag(covariance1))
print('a2 =', params1[0], '±', errors1[0])
print('b2 =', params1[1], '±', errors1[1])
params2, covariance2 = curve_fit(f=g, xdata=noms(t2m), ydata=noms(p2mf))
errors2 = np.sqrt(np.diag(covariance2))
print('a2 =', params2[0], '±', errors2[0])
print('b2 =', params2[1], '±', errors2[1])
params3, covariance3 = curve_fit(f=g, xdata=noms(t3m), ydata=noms(p3mf))
errors3 = np.sqrt(np.diag(covariance3))
print('a3 =', params3[0], '±', errors3[0])
print('b3 =', params3[1], '±', errors3[1])
params4, covariance4 = curve_fit(f=g, xdata=noms(t4m), ydata=noms(p4mf))
errors4 = np.sqrt(np.diag(covariance4))
print('a4 =', params4[0], '±', errors4[0])
print('b4 =', params4[1], '±', errors4[1])
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
Stheo=77
Sg=ufloat(11.9 , 0.9)
Slin1=ufloat(12.2 , 0.9)
Slin2=ufloat(9.1 ,0.7)
h1=(((2*10**(-2))+(2*10**(-5)))/2)
h2=h1-(2*10**(-5))
h3=((2*10**(-4)+8*10**(-3))/2) #s1 2.4
h4=h3-(2*10**(-4))
h5=((10**(-4))+(8*10**(-3)))/2 #s2 1.4
h6=h5-(10**(-4))
h7=((3*10**(-5))+(2*10**(-3)))/2
h8=h7-(3*10**(-5))
h9=((8*10**(-5))+(8*10**(-3)))/2
h10=h9-(8*10**(-5))

plt.errorbar(noms(t1m),noms(p1mf), xerr=stds(t1m), yerr=stds(p1mf),fmt='rx',linewidth=1, label='Messdaten')
plt.plot(xlin, g(xlin,*params1),'k-',linewidth=2, label='Fit')
plt.xlim(-1.5,19)
#plt.ylim(-0.001,0.004)
plt.grid(True)
plt.ylabel(r'$p\,/\,$mbar')
plt.xlabel(r'$t\,/\,$s')
plt.legend(loc='best')
plt.savefig('LeckrateTurbo2_4.pdf')
plt.clf()
plt.errorbar(noms(t2m),noms(p2mf), xerr=stds(t2m), yerr=stds(p2mf),fmt='rx',linewidth=1, label='Messdaten')
plt.plot(xlin, g(xlin,*params2),'k-',linewidth=2, label='Fit')
plt.xlim(-2.5,34)
#plt.ylim(-0.001,0.004)
plt.grid(True)
plt.ylabel(r'$p\,/\,$mbar')
plt.xlabel(r'$t\,/\,$s')
plt.legend(loc='best')
plt.savefig('LeckrateTurbo1_4.pdf')
plt.clf()
plt.errorbar(noms(t3m),noms(p3mf), xerr=stds(t3m), yerr=stds(p3mf),fmt='rx',linewidth=1, label='Messdaten')
plt.plot(xlin, g(xlin,*params3),'k-',linewidth=2, label='Fit')
plt.xlim(-5,70)
#plt.ylim(-0.001,0.004)
plt.grid(True)
plt.ylabel(r'$p\,/\,$mbar')
plt.xlabel(r'$t\,/\,$s')
plt.legend(loc='best')
plt.savefig('LeckrateTurbo3_5.pdf')
plt.clf()
plt.errorbar(noms(t4m),noms(p4mf), xerr=stds(t4m), yerr=stds(p4mf),fmt='rx',linewidth=1, label='Messdaten')
plt.plot(xlin, g(xlin,*params4),'k-',linewidth=2, label='Fit')
plt.xlim(-5,60)
#plt.ylim(-0.001,0.004)
plt.grid(True)
plt.ylabel(r'$p\,/\,$mbar')
plt.xlabel(r'$t\,/\,$s')
plt.legend(loc='best')
plt.savefig('LeckrateTurbo8_5.pdf')
plt.clf()

plt.errorbar(h1,noms(Stheo), xerr=h2, yerr=0,fmt='kx',linewidth=1, label=r'$S_{theo}$')
plt.errorbar(h1,noms(Sg), xerr=h2, yerr=stds(Sg),fmt='rx',linewidth=1, label=r'$S_{g}$')
#plt.errorbar(125,noms(Sexp), xerr=0, yerr=stds(Sexp),fmt='rx',linewidth=1, label=r'$S_{g}$')
plt.errorbar(0.0011,noms(Slin1), xerr=0.0009, yerr=stds(Slin1),fmt='bx',linewidth=1, label=r'$S_{lin1}$')
#plt.plot((0.2,1), (noms(Slin2),noms(Slin2)), 'g-', label=r'$S_{lin2}$')
plt.errorbar(0.00004,noms(Slin2), xerr=0.00002, yerr=stds(Slin2),fmt='gx',linewidth=1, label=r'$S_{lin2}$')
#plt.plot((0.02,0.1), (noms(Slin3),noms(Slin3)), 'm-', label=r'$S_{lin3}$')
#plt.plot((0.1,0.8), (noms(S2),noms(S2)), 'y-', label=r'$S_{leck1}$')
plt.errorbar(h3,noms(S2), xerr=h4,yerr=stds(S2),fmt='yx',linewidth=1, label=r'$S_{leck1}$')
plt.errorbar(h5,noms(S1), xerr=h6, yerr=stds(S1),fmt='cx',linewidth=1, label=r'$S_{leck2}$')
#plt.plot((0.4,4), (noms(S3),noms(S3)), 'b.', label=r'$S_{leck3}$')
plt.errorbar(h7,noms(S3), xerr=h8, yerr=stds(S3),fmt='mx',linewidth=1, label=r'$S_{leck3}$')
#plt.plot((0.8,8), (noms(S4),noms(S4)), 'r.', label=r'$S_{leck4}$')
plt.errorbar(h9,noms(S4), xerr=h10, yerr=stds(S4),fmt='rx',linewidth=1, label=r'$S_{leck4}$')
plt.xlim(0,0.02)
plt.grid(True)
plt.xlabel(r'$p\,/\,$mbar')
plt.ylabel(r'$S\,/\,\frac{l}{s}$')
plt.legend(bbox_to_anchor=(1,0.8))
plt.savefig('SaugverTurbo.pdf')
plt.clf()
