import scipy.integrate as integrate, numpy as np, math, cmath, matplotlib.pyplot as plt

fichier=open("valeurs.txt")
x,y=[],[]
t=0
for ligne in fichier:
    x.append(t)
    y.append(float(ligne.replace(",",".")))
    t+=0.1
fichier.close

fichier=open("Pulsations_37.txt")
w=[]
for ligne in fichier:
    w.append(float(ligne)*np.pi/180)
fichier.close

fichier=open("Valeurs2017+2018.txt")
x2,y2=[],[]
t=0
for ligne in fichier:
    x2.append(t)
    y2.append(float(ligne.replace(",",".")))
    t+=0.1
fichier.close

fichier=open("Dates.txt")
date=[]
for ligne in fichier:
    date.append(ligne)
fichier.close

def jour(n):
    """Transforme un numéro de jour n en une date"""
    return date[n]

def y_moy(y):
    res=0
    for i in y:
        res+=i
    return res/len(y)

def trapz(y,x,T):
    res=0
    h=x[1]-x[0]
    for k in range(T-1):
        res+=y[k]+y[k+1]
    return res*h/2

y0=y_moy(y)

def Complexe_liste(y,wi):
    T=len(x) #Donne l'indice du temps final t[-1]
    y0=y_moy(y)
    cos=[(y[t]-y0)*np.cos(wi*x[t]) for t in range(T)]
    sin=[(y[t]-y0)*np.sin(wi*x[t]) for t in range(T)]
    return [trapz(cos,x,T)*2/x[-1],trapz(sin,x,T)*2/x[-1]]  #x[-1]=temps total

def parametres_l():
    A,P=[],[]
    for wi in w:
        C=Complexe_liste(y,wi)
        a=np.sqrt(C[0]**2+C[1]**2)
        p=np.arctan(C[1]/C[0])
        A+=[a]
        P+=[p]
    return [A,P]

para=parametres_l()
        
def h(t):
    res=y0
    Amp,Phi=para[0],para[1]
    for i in range(len(w)):
        res+=Amp[i]*np.cos(w[i]*t-Phi[i])
    return res
    
def ecart_max(j1,j2):
    """Recherche l'écart maximal entre la courbe réelle et celle prédite sur l'intervalle [j1,j2] (en jour)"""
    p1,p2=j1*240,j2*240  #Correspond à la position de j1 et j2 dans x2 (on a 24*10 points par jours)
    x3=x2[p1:p2+1]
    Y=y2[p1:p2+1]
    H = [h(t) for t in x3]
    X = [t/24 for t in x3]
    tmax=0  
    ecart_max=abs(H[tmax]-Y[tmax])
    for t in range(len(x3)):
        ecart=abs(H[t]-Y[t])
        if ecart>ecart_max:
                  tmax=t
                  ecart_max=ecart
    return [x3[tmax]/24,ecart_max]
        
def nb_ecart(j1,j2,eps):
    """Compte le nombre d'écarts dépassant l'ecart limite eps sur l'intervalle [j1,j2]"""
    p1,p2=j1*240,j2*240  #Correspond à la position de j1 et j2 dans x2 (on a 24*10 points par jours)
    x3=x2[p1:p2+1]
    Y=y2[p1:p2+1]
    H = [h(t) for t in x3]
    X = [t/24 for t in x3]
    N=0 #compteur
    for t in range(len(x3)):
        ecart=abs(H[t]-Y[t])
        if ecart>=eps:
                  N+=1
    return N  

def Trace(j1,j2):
    ''' Trace la fonction prédite et celle mesurée entre les jours j1 et j2'''
    H = [h(t) for t in x2]
    X = [t/24 for t in x2]
    
    plt.plot(X,H,'r',label='Courbe prédite')
    plt.plot(X,y2,'b',label='Courbe réelle')
    plt.xlim(j1,j2)
    plt.ylabel('Hauteur des marées (en m)')
    plt.xlabel('Jour')
    plt.legend()
    plt.grid()
    plt.show()
