import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
zero = const.zero_Celsius 
def mittel(x):              #the real mean()-ing of life
    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100
# Widerstände in Temperataturen umrechnen
def T(R):
    return 0.00134 * R**2 + 2.296 * R - 243.02
# cp berechnen
def cp(dT, U, I, M, dt, m):
    return (U * I * dt * M)/(dT * m)
# cv bestimmen
def cv(cp, alpha_T, kappa, V0, Tbar):
    return -9*alpha_T**2*kappa*V0*Tbar + cp
def alpha(T):
	alph = np.empty(0)
	for t in T:
		for a in alphaT:
			if(a>=t):
				alph=np.append(alph,alphaW[alphaT==a])
				break
		if(t>alphaT[-1]):
			alph=np.append(alph,alphaW[-1])

	return alph*10**(-6)	
#Werte einlesen 
t,RP,RG,UP,IP,UG,IG = np.genfromtxt('Messwerte.txt', unpack = True)
Masse = 0.342 #[kg]
MolMasse = 63.546*const.N_A*const.value("atomic mass constant") #[kg/mol] aus Quelle 
kappa = 137.8*10**9 #Pa Quelle
alphaT, alphaW = np.genfromtxt('ausdehnungskoeffizient.txt', unpack = True) # Anleitung
rho = 8.96*100**3/1000 # Dichte [kg/m³] bei 20 Celsius Quelle
MVol = MolMasse / rho
print("Masse in kg:", Masse)
print("Molmasse in kh /mol:", MolMasse)
print("Kompressionsmodul uas Quelle in Pa:",kappa)
print("Dichte aus Quelle in kg /m ³:" , rho)
print("Molvolumen in mol/m³:", MVol)

#Fehler
t*=60 #in Sekunden
t = unp.uarray(t,5) # 5 sec Fehler Faktor Mensch 
RP = unp.uarray(RP,RP*0.002) #bei beiden Wiederständen 2% Fehler 
RG = unp.uarray(RG,RG*0.002)
UP = unp.uarray(UP,0.001) #10% der Skalar xx,xx [Volt]
IP = unp.uarray(IP,0.01) #10% der Skalar  xxx,x [mA]
# Umrechnen
TP =  T(RP) + zero #in Kelvin
TG = T(RG) + zero
IP *= 10**(-3) # [mA] zu [A]

#Delta 
DeltaTP = np.empty(0)
errDeltaTP = np.empty(0)
Deltat = np.empty(0)
errDeltat= np.empty(0)
laufi = range(len(TP)-1) #laufi der nicht pythonische Laufindex
for i in laufi:
	Deltat= np.append(Deltat,noms(t[i+1] - t[i]))
	errDeltat = np.append(errDeltat,stds(t[i+1] - t[i]))
	DeltaTP = np.append(DeltaTP,noms(TP[i+1] - TP[i]))
	errDeltaTP = np.append(errDeltaTP, stds(TP[i+1] - TP[i]))

DeltaTP = unp.uarray(DeltaTP,errDeltaTP)
dt = unp.uarray(Deltat,errDeltat)
#Berechne Cp und Cv:
x = noms(TP[1:])
Cp =  cp(DeltaTP,UP[1:],IP[1:],MolMasse,dt, Masse)
alphax = alpha(x)
Cv = cv(Cp,alpha(x),kappa,MVol,TP[:-1])
mitCp = np.mean(Cp)
#print(alphax)
#Tabelle
np.savetxt('CPtab.txt',np.column_stack([noms(DeltaTP),stds(DeltaTP),noms(UP[1:]),stds(UP[1:]),noms(IP[1:]),stds(IP[1:]),noms(dt),stds(dt),noms(Cp), stds(Cp)]), delimiter=' & ',newline= r'\\'+'\n' )
np.savetxt('CVtab.txt',np.column_stack([x,stds(TP[1:]),noms(Cp), stds(Cp),alphax, noms(Cv), stds(Cv)]), delimiter=' & ',newline= r'\\'+'\n' )


#Plot des Temperatur verlaufes
#plt.plot(noms(TP),noms(TP), 'ro-', label='Probe')
#plt.plot(noms(TP),noms(TG), 'bo-', label ='Gefäß')
#plt.ylabel(r'$T \; / \; K$')
#plt.xlabel(r'$Temperatur \; der \, Probe \; T_P \; / \;K $')
#plt.title('Temperaturverlauf')
#plt.xlim(70,310)
#plt.grid()
#plt.legend(loc='best')
##plt.show() 
#plt.savefig("Tempv.pdf")
#Plot der spez Wärme 
#print(cv(Cp,alpha(x),kappa,MVol,TP[:-1]))
#print(Cp)
##plt.subplot(1, 2, 1)
#plt.errorbar(x,noms(cv(Cp,alpha(x),kappa,MVol,x)),stds(cv(Cp,alpha(x),kappa,MVol,TP[1:])), stds(TP[1:]) ,fmt='rx',ecolor="red",label=r'$ C_V $' )
#plt.errorbar(x,noms(Cp),stds(Cp), stds(TP[:-1]) ,fmt='bx',ecolor="blue", label= r'$C_p$')
#plt.axhline(y=3*const.R, xmin=0,xmax=350, color="k", label=r'$3 \cdot R$')
#plt.axvline(x=170, ymin=0,ymax=35,ls='--', color='g', label=r'$ 170 K$')
#plt.plot(x,noms(cv(Cp,alpha(x),kappa,MVol,x)) ,'rx', label='Kurve')
#plt.plot(x,noms(Cp),'bx')
#plt.xlabel(r'$T \; / \; K$')
#plt.ylabel(r'$spezifische \; Wärme \; C \; / \; \frac{J}{mol\cdot K}$')
#plt.legend(loc='best')
##plt.show()
#plt.savefig('Cplot.pdf')

#Teil c Debye Kurve:
Tthet = TP[1:] 
Tthet = Tthet[Tthet<170] # alle TP unter 170K
cvd = noms(Cv[TP[1:]<170])
#print(Cv[TP[1:]<170]) # alle  Cv werte unter 170K
thetaT = np.array([4.1, 3.0 ,2.9 , 3.1 , 2.6 , 2.2 , 2.2 , 2 , 2.2 , 2.4 , 1.4 ]) #ausgewählte Theta / T werte
theta = thetaT*Tthet
#print(cvd, thetaT)
#print(" theta aus Cv und T:", theta)
#Tabelle
np.savetxt('thetatab.txt',np.column_stack([noms(Cv[TP[1:]<170]), stds(Cv[TP[1:]<170]),thetaT,noms(Tthet),stds(Tthet),noms(theta),stds(theta)]), delimiter=' & ',newline= r'\\'+'\n' )

##Fit
#def fit(C,a,b):
#	return a*C +b 
#params , cov = curve_fit(fit , cvd ,thetaT )
#params = correlated_values(params, cov)
#af = params[0]
#bf = params[1]
#print(af,bf)
##plot theta / T in Abh von Cv
#y = np.linspace(11,25,1000)
#plt.plot(cvd,thetaT, 'rx')
#plt.plot(y, noms(fit(y,af,bf)),'b--', label= 'Ausgleichsgerade')
#plt.title(r'$ Debye-Temperatur durch T in Abhängigkeit von  C_V$')
#plt.ylabel(r'$\theta_D / T $')
#plt.xlabel(r'$C_V \; / \; \frac{J}{mol\cdot K} $')
##plt.legend(loc='best')
#plt.show() 
#

#Teil d berechne omega d und theta d aus dem intergral 
NL = (Masse/MolMasse) * const.N_A #Teilchenzhal in der Probe gleich NL 
print("Teilchenzahl in der Probe NL:",NL)
Vol = MVol * (Masse/MolMasse) # Volumen ist gleich L³
print("Volumen der Probe L³:",Vol)
omegaD = (((18*NL*np.pi**2)/Vol)*(1/((1/(4.7*1000)**3) + (2/(2.26*1000)**3))))**(1/3)
print("omegaD in Herz:",omegaD)
print("theta aus omega und Zw:",omegaD*(const.hbar/const.k))
print("Mittelwert der bestimmten theta:",np.mean(theta))
print("Rel. Fehler von dem berechneten omega:",relf(omegaD*(const.hbar/const.k),np.mean(theta)))



#Tabelle
# np.savetxt('tab.txt',np.column_stack([x,y]), delimiter=' & ',newline= r'\\'+'\n' )

