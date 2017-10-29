import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
x = np.linspace(0, 10, 1000)
mhub = const.value('Bohr magneton') #das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt
n=3.3543 # Brechungsindex für GaAs bei einer Wellenlänge von 1771.14  
A= (const.e**3 * 0.43) / ( 8*(np.pi**2) * const.epsilon_0*(const.c**3) *n)  #konst zum bestimmen der eff Masse 
def mittel(x):              #the real mean()-ing of life
    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100

# d = dicke= Länge im Strahlengang
hd = 5.11*10**-3
d1 = 1.296*10**-3
d2 = 1.36*10*10**-3
#Werte einelsen l = Wellenläng. tb= teta mit bfeld , to= teta ohne (b/o)bm= bogenmin mit und ohne bfeld
#hochreines GaAg
hl , htb , hbbm , hto , hobm = np.genfromtxt('hochreinesGaAs.txt', unpack = True)
#1. n-dotiertes GaAs 
l1 , tb1, bbm1 , to1 , obm1 = np.genfromtxt('1.n-dotiertesGaAs.txt', unpack = True) 
#2. n-dotiertes GaAs
l2  , tb2 , bbm2 , to2 , obm2 = np.genfromtxt('2.n-dotiertesGaAs.txt', unpack = True) 
#jetz noch die Winkel um und zusammenrechnen. dann die ohne bFeld von den mit BFeld abziehen. Alle Messwerte in ein Plot mit 'o-'
#Winkel zusammenrechnen eine bogenmin = (1°/60) Grad 
u = 1/60 #u der kleine Umrechnungsfaktor, weil faul
htb+= hbbm*u
hto+= hobm*u
tb1+= bbm1*u
to1+= obm1*u
tb2+= bbm2*u
to2+= obm2*u
# The True Teta ausrechnen t mit BFeld - t ohne
ht = htb - hto 
t1 = tb1 - to1 
t2 = tb2 - to2 
ht = abs(ht)
t1 = abs(t1)
t2 = abs(t2)
#Umrechnen der Wellenlängen
hl*=10**-6
l1*=10**-6
l2*=10**-6
#Differenzzwischen der hochreinen und den zwei Proben
D1 = abs(ht/hd-t1/d1) 
D2 = abs(ht/hd-t2/d2)
D = np.array([D1[0:2]])
D = np.append(D,D1[4:])
L = np.array([hl[0:2]])
L = np.append(L,hl[4:])
print('Mittelwert der Wellenlängen:' , np.mean(hl))
#Fit
C1= A*2.8*10**16 #weil umgerechnet auf Meter 
C2= A*1.2*10**16
def fitf1(x,a,b):
	return x*a + b
def fitf2(l,m,b):
	return C2*(1/m)**2 * l**2 +b 
params , cov = curve_fit(fitf1 ,hl**2,D2)
params = correlated_values(params, cov)
#params2 , cov2 = curve_fit(fitf2 ,hl,D2 )
#params2 = correlated_values(params2, cov2)
for p in params:
	print(p)  
a =  params[0]
b = params[1]

#Tabelle
# np.savetxt('tab.txt',np.column_stack([x,y]), delimiter=' & ',newline= r'\\'+'\n' )
plt.plot(L**2 , D, 'ro', label= 'Differenz mit dem 1.n-dotiertem')
#plt.plot(hl**2,fitf1(hl**2,noms(a)),'b--')
plt.plot(hl**2, D2, 'bo', label='Differenz zwischen dem 2.n-dotiertem')
plt.xlabel(r'Wellenlänge zum Quadrat $\lambda^2 \:/\: m^2 $')
plt.ylabel(r'$\frac{\Delta\theta}{L} / \degree m^{-1} $')
#plt.legend(loc='best')
plt.show()
#plt.savefig('DeltaTheta.pdf')
#plt.title('Farady-Rotation der GaAs-Proben im Vergleich')
#plt.plot((l2**2),t2/d2,'go',label='n-dotiertes GaAs mit $N=1,2\cdot10^{18}cm^{-3}$')
#plt.plot((l1**2),t1/d1,'bo',label=r'n-dotiertes GaAs mit $ N=2,8\cdot10^{18}cm^{-3}$')
#plt.plot((hl**2),ht/hd,'ro', label=r'Hochreines GaAs')
#plt.xlabel(r'Wellenlänge zum Quadrat $\lambda^2 \:/\: m^2 $')
#plt.ylabel(r'Winkelrotation normiert auf Probenlänge $\frac{\theta}{L} \:/\:\degree m^{-1}$')
#plt.legend(loc='best')
#plt.savefig('GaAsimVgl.pdf')
#plt.show()



#plt.title('Hochreines GaAs')
#plt.plot(hl**2,ht,'ro', label=r'Fraday-Rotation gg. Wellenlänge')
#plt.xlabel(r'$\lambda \:/\:\mu m $')
#plt.ylabel(r'$\Delta\theta \:/\:\degree $')
#plt.legend(loc='best')
#plt.savefig('hochreinesGaAsplot.pdf')
#plt.show()

#plt.clf()

#plt.title(r'n-dotiertes GaAs mit $ N=2,8\cdot10^{18}cm^{-3}$')
#plt.plot(l1**2,t1,'bo',label=r'Fraday-Rotation gg. Wellenlänge')
#plt.xlabel(r'$\lambda \:/\:\mu m$')
#plt.ylabel(r'$\Delta\theta \:/\:\degree $')
#plt.legend(loc=r'best')
#plt.savefig('1.nGaAsplot.pdf')
#plt.show()

#plt.clf()

#plt.title('n-dotiertes GaAs mit $N=1,2\cdot10^{18}cm^{-3}$')
#plt.plot(l2**2,t2,'go',label='Fradday Rotation gg. Wellenlänge')
#plt.xlabel(r'$\lambda \:/\:\mu m $')
#plt.ylabel(r'$\Delta\theta \:/\:\degree $')
#plt.legend(loc='best')
#plt.savefig('2.nGaAsplot.pdf')
#plt.show()


