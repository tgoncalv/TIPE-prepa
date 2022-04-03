import scipy.integrate as integrate, numpy as np, math, cmath, matplotlib.pyplot as plt

fichier=open("valeurs.txt")
t,H=[],[]
temps=0
for ligne in fichier:
    t.append(temps)
    H.append(float(ligne.replace(",",".")))
    temps+=0.1 #6 minutes=0.1 heure
fichier.close

fichier=open("Pulsations_37.txt")
w=[]
for ligne in fichier:
    w.append(float(ligne)*np.pi/180)
fichier.close

fichier=open("Valeurs2017+2018.txt")
x2,y2=[],[]
temps=0
for ligne in fichier:
    x2.append(temps)
    y2.append(float(ligne.replace(",",".")))
    temps+=0.1
fichier.close

fichier=open("Dates.txt")
date=[]
for ligne in fichier:
    date.append(ligne)
fichier.close

def jour(n):
    """Transforme un numéro de jour n en une date"""
    return date[n]

N=len(t)
M=len(w)

def C(k):
    res=0
    wk=w[k]
    for i in range(N):
        res+=np.cos(wk*t[i])
    return res
def S(k):
    res=0
    wk=w[k]
    for i in range(N):
        res+=np.sin(wk*t[i])
    return res
def CC(k,j):
    res=0
    wk,wj=w[k],w[j]
    for i in range(N):
        res+=np.cos(wk*t[i])*np.cos(wj*t[i])
    return res
def SS(k,j):
    res=0
    wk,wj=w[k],w[j]
    for i in range(N):
        res+=np.sin(wk*t[i])*np.sin(wj*t[i])
    return res
def CS(k,j):
    res=0
    wk,wj=w[k],w[j]
    for i in range(N):
        res+=np.cos(wk*t[i])*np.sin(wj*t[i])
    return res

def B11(): #matrice de taille (M+1)x(M+1)
    Mat=[N]+[C(j) for j in range(M)] #Matrice B11 constitué seulement de la première ligne
    for i in range(M):
        lig=[C(i)]+[CC(i,j) for j in range(M)]  #ligne i
        Mat+=lig #on rajoute la ligne i à Mat
    return np.reshape(Mat,(M+1,M+1))

def B12(): #matrice de taille (M+1)x(M)
    Mat=[S(j) for j in range(M)]
    for i in range(M):
        lig=[CS(i,j) for j in range(M)]
        Mat+=lig
    return np.reshape(Mat,(M+1,M))

def B22(): #matrice de taille (M)x(M)
    Mat=[]
    for i in range(M):
        lig=[SS(i,j) for j in range(M)]
        Mat+=lig
    return np.reshape(Mat,(M,M))

def B():
    A=np.concatenate((B11(),B12()),axis=1)
    B=np.concatenate((np.transpose(B12()),B22()),axis=1)
    return np.concatenate((A,B),axis=0)

def Y(H): #BX=Y avec Y dépendant des valeurs H mesurées
    res=0
    for i in range(N):
        res+=H[i]
    Mat=[res]
    for lig in range(M):
        res=0
        w_lig=w[lig]
        for i in range(N):
            res+=H[i]*np.cos(w_lig*t[i])
        Mat+=[res]
    for lig in range(M):
        res=0
        w_lig=w[lig]
        for i in range(N):
            res+=H[i]*np.sin(w_lig*t[i])
        Mat+=[res]
    return np.reshape(Mat,(2*M+1,1))

def parametres(H):
    X=np.linalg.solve(B(),Y(H)).tolist()
    C0=X[0] #Valeur moyenne
    a,b=X[1:M+1],X[M+1:] #Renvoie la liste des a et b
    Amp,Phi=[],[]
    for i in range(M):
        Amp+=[np.sqrt(a[i][0]**2+b[i][0]**2)]
        Phi+=[np.arctan(b[i][0]/a[i][0])*180/np.pi]
    return [C0,Amp,Phi]

#Sol=parametres(H) (trop long, le résultat obtenu est celui ci-dessous)
Sol=[[2.120673713125022], [1.0992755526136821, 0.26844099600035809, 0.21999657185845872, 0.75123465306451209, 0.022701245238320813, 0.3901994193096599, 0.011125302036567123, 0.033628892307596525, 0.0020160904280276605, 0.0090825605634479688, 0.045792350315563869, 0.00052113388020196211, 0.037050056286104506, 0.029179263407002411, 0.015802197049013568, 0.016292449872538122, 0.01506044364931885, 0.038409581841495427, 0.032892022908085489, 0.043317708887760542, 0.029890175512935248, 0.088860213967523194, 0.022577006553855472, 0.025590578076850633, 0.0088930163393407587, 0.06614044665810348, 0.016608372434884604, 0.0023996783921557906, 0.0054593675582169177, 0.25543666228395967, 0.007579807337348414, 0.0032686904010345063, 0.050866346628601684, 0.030697136169299311, 0.059120987099488369, 0.00082620676278820978, 0.012166053923803273], [84.995395465130741, 36.403804358927317, -44.009579611320014, -87.483128822427773, -11.605401778215843, -27.856469785420305, -5.9322262558770591, -32.161761602292145, 65.190416727076737, 46.955727032679995, 59.489573471241329, -14.268562285114992, 22.35700243848375, 1.8050988755478148, 74.845533437314188, -47.481555923163235, 30.138619202904852, 40.416783984833565, 43.531903506231416, 74.125481647562779, -70.135340606431171, 46.721899693896916, 41.768332567017374, 18.464882803921544, -49.746787103678052, 46.306735388745793, 27.545757232852768, 44.519696233256809, -80.329688597514377, -72.187521595289581, 32.793253981877335, -9.1801649791273459, 50.291112459285799, 21.12236951890667, 25.628343834880802, -24.201240288075869, -53.17368993287689]]

def h_estime(t):
    res=Sol[0]
    Amp2=Sol[1]
    Phi2=[np.pi*p/180 for p in Sol[2]]  #problème de déphasage à étudier
    for i in range(len(w)):
        res+=Amp2[i]*np.cos(w[i]*t-Phi2[i])
    return res

def ecart_max(j1,j2):
    """Recherche l'écart maximal entre la courbe réelle et celle prédite sur l'intervalle [j1,j2] (en jour)"""
    p1,p2=j1*240,j2*240  #Correspond à la position de j1 et j2 dans x2 (on a 24*10 points par jours)
    x3=x2[p1:p2+1]
    Y=y2[p1:p2+1]
    H2 = [h_estime(t) for t in x3]
    X = [t/24 for t in x3]
    tmax=0  
    ecart_max=abs(H2[tmax]-Y[tmax])
    for t in range(len(x3)):
        ecart=abs(H2[t]-Y[t])
        if ecart>ecart_max:
                  tmax=t
                  ecart_max=ecart
    return [x3[tmax]/24,ecart_max]

def nb_ecart(j1,j2,eps):
    """Compte le nombre d'écarts dépassant l'ecart limite eps sur l'intervalle [j1,j2]"""
    p1,p2=j1*240,j2*240  #Correspond à la position de j1 et j2 dans x2 (on a 24*10 points par jours)
    x3=x2[p1:p2+1]
    Y=y2[p1:p2+1]
    H2 = [h_estime(t) for t in x3]
    X = [t/24 for t in x3]
    N=0 #compteur
    for t in range(len(x3)):
        ecart=abs(H2[t]-Y[t])
        if ecart>=eps:
                  N+=1
    return N

def Trace(j1,j2):
    ''' Trace la fonction prédite et celle mesurée entre les jours j1 et j2'''
    X = [tps/24 for tps in x2]
    H_estime = [h_estime(tps) for tps in x2]
    
    plt.plot(X,y2,'b',label='Courbe réelle')
    plt.plot(X,H_estime,'r',label='Courbe prédite')
    plt.xlim(j1,j2)
    plt.ylabel('Hauteur des marées (en m)')
    plt.xlabel('Jour')
    plt.legend()
    plt.grid()
    plt.show()



