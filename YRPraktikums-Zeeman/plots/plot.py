import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
x = np.linspace(0, 21, 1000)
mhub = const.value('Bohr magneton') #das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt
c = const.c 
h = const.h
def mittel(x):              #the real mean()-ing of life
	return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
	return (np.absolute(l-m)/l)*100
# BFeld rauf	
Ia,Ba = np.genfromtxt('B-Feld_rauf.txt', unpack = True)
#BFeld runter
Iu,Bu = np.genfromtxt('B-Feld_runter.txt', unpack = True)
# Fit Arrays bis (13,825) und Plot arrays 
Bf = Ba[Ba<825] 
Bp = Ba[Ba>825]
If = Ia[Ia<14]
Ip = Ia[Ia>13]
#Fit
def f(x,a,b):
	return a*x +b
params , cov = curve_fit(f ,If, Bf )
params = correlated_values(params, cov)
for p in params:
	print(p)
a = params[0]
b = params[1]

plt.plot(Iu,Bu, 'gx', label= 'Absteigende Hysteresekurve')
plt.plot(If,Bf ,'rx', label='Aufsteigende Hyysteresekurve')
plt.plot(Ip,Bp ,'kx', label='Nicht-lin. Teil der aufsteig. Hysteresekurve')
plt.plot(x, f(x,noms(a),noms(b)), 'b-', label= 'Ausgleichsgerade d. linearen Teils')
plt.xlim(0,21)
plt.ylim(0,1500)
plt.xlabel(r'$Feldstrom \; I \:/\: A$')
plt.ylabel(r'$Magnetische Feldstärke \; B \:/\: mT$')
plt.legend(loc='best')
#plt.show()
plt.savefig('BFeldplot.pdf')

# del lambda funktion:
def dl(ds,Ds,Dl):
	return 0.5 * (ds/Ds) * Dl
# gij funktion:
def gij(B,l,dl):
	return (h*c/l**2) * (1/(mhub*B)) * dl
# energie fkt:
def E(l,dl):
	return h *c * dl/l**2
#rote werte 
#ohne pol
r = np.genfromtxt('644nm_B=0_P=0.txt', unpack = True)
#mit pol = 0 aufspaltung 
r1, r2 = np.genfromtxt('644nm_I=9.5_P=0.txt' , unpack = True)
# BFeld der Aufspaltung rot I = 9.5 P=0
Bra = f(9.5,a,b)
#delta lambda d für rot 638,8 nm
Dlr = (643.8e-9**2/(2*4e-3))*np.sqrt(1/(1.4567**2 -1))
# Berechnung delta s
Dsr = np.empty(0)
q = range(len(r[1:]))
for i in q:
	Dsr= np.append(Dsr,[r[i+1]-r[i]])
DS1=Dsr
Dsr = mittel(Dsr)
# Berechnung del s
dsr = np.array(r2-r1)
dsr = mittel(dsr)
#Berechnung del lambda 
dlr = dl(dsr,Dsr,Dlr)
#Txt ausgaben 
print(" Für Rot 643.8 nm ")
print(" Delta Lambda D LG-Platte:" , Dlr)
print(" Gangunterschied:", Dsr)
print(" del s der Auspaltung", dsr)
print("BFeld der Aufspaltung:", Bra)
print(" del Lambda:", dlr)
print(" gij Faktor:", gij(Bra*10**-3, 643.8e-9,dlr))
print("Energie:",E(643.8e-9,dlr/const.e))
print("Relativer Fehler von Erwartungswert 1:", relf(1,gij(Bra*10**-3, 643.8e-9,dlr)))

#blau 480nm
#blau P = 0 Werte
b1 = np.genfromtxt('480nm_B=0_P=0.txt', unpack = True)
bn1,bn2 = np.genfromtxt('480nm_I=5.5_P=0.txt', unpack = True)
#blau P = 90 Werte 
b2 = np.genfromtxt('480nm_B=0_P=90.txt', unpack = True)
ba1, ba2 = np.genfromtxt('480nm_B=0.981mT_P=90.txt' ,unpack = True)
# BFeld der anregung bei P = 0
Bn = f(5.5,a,b)*10**-3 
# Delta lambda D für blau 480nm
Dlb = (480e-9**2/(2*4e-3))*np.sqrt(1/(1.4635**2 -1))
#Berechnung von Delta s für P = 0
Dsb1 = np.empty(0)
g = range(len(b1[:-1]))
for i in g:
	Dsb1 = np.append(Dsb1,b1[i+1]-b1[i])
DS2 = Dsb1
Dsb1 = mittel(Dsb1)
#Berechnung von Dekta s für P = 90
Dsb2 = np.empty(0)
g = range(len(b2[:-1]))
for i in g:
	Dsb2 = np.append(Dsb2, b2[i+1] - b2[i])
DS3= Dsb2
Dsb2 = mittel(Dsb2)
#Berechnung von del s für P = 0
dsb1 = bn2 -bn1
dsb1 = mittel(dsb1)
#Berechnung von del s für P = 90
dsb2 = ba2 - ba1
dsb2 = mittel(dsb2)
#Berechnung von del lambda für P = 0
dl1 = dl(dsb1,Dsb1,Dlb)
#Berechnung von del lambda für P = 90
dl2 = dl(dsb2,Dsb2,Dlb)
#Txt ausgaben 
print(" Für 480 nm ")
print(" Delta Lambda D LG-Platte:" , Dlb)
print(" Gangunterschied für P= 0:", Dsb1)
print(" del s der Auspaltung für P = 0", dsb1)
print("BFeld der Aufspaltung für P = 0:", Bn)
print(" del Lambda für P= 0:", dl1)
print("Energie für P=0:",E(480e-9,dl1)/const.e)
print(" gij Faktor für P = 0:", gij(Bn, 480e-9,dl1))
print("Relativer Fehler von Erwartungswert 1.75:", relf(1.75,gij(Bn, 480e-9,dl1)))
print(" Gangunterschied für P=90:", Dsb2)
print(" del s der Auspaltung für P =90", dsb2)
print("BFeld der Aufspaltung für P = 0:",0.981)
print(" del Lambda für P=90:", dl2)
print("Energie für P=0:",E(480e-9,dl2)/const.e)
print("Relativer Fehler von Erwartungswert 0.5:", relf(0.5, gij(0.981, 480e-9,dl2)))
print(" gij Faktor für P =90:", gij(0.981, 480e-9,dl2))

##Tabelle
np.savetxt('rotP0PPtab.txt',np.column_stack([r]), delimiter=' & ',newline= r'\\'+'\n' )
np.savetxt('rotP0DStab.txt',np.column_stack([DS1]), delimiter=' & ',newline= r'\\'+'\n' )
np.savetxt('rotP0PtPtab.txt',np.column_stack([r1,r2,r2-r1]), delimiter=' & ',newline= r'\\'+'\n' )
np.savetxt('blauP0PPtab.txt',np.column_stack([b1]), delimiter=' & ',newline= r'\\'+'\n' )
np.savetxt('blauP0DStab.txt',np.column_stack([DS2]), delimiter=' & ',newline= r'\\'+'\n' )
np.savetxt('blauP0PtPtab.txt',np.column_stack([bn1,bn2,bn2-bn1]), delimiter=' & ',newline= r'\\'+'\n' )
np.savetxt('blauP90PPtab.txt',np.column_stack([b2]), delimiter=' & ',newline= r'\\'+'\n' )
np.savetxt('blauP90DStab.txt',np.column_stack([DS3]), delimiter=' & ',newline= r'\\'+'\n' )
np.savetxt('blauP90PtPtab.txt',np.column_stack([ba1,ba2,ba2-ba1]), delimiter=' & ',newline= r'\\'+'\n' )

