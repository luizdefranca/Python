# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------
#      Fonctions Partie I
#------------------------------------
def step_Euler(y,t,h,f):
    """Paramètres :
    y : le résultat de l'estimation du pas précédent
    t : le point où est évaluée y
    h : la valeur du pas
    f : la fonction de l'équa diff
    Retourne l'estimation de la fonction f pour le pas suivant"""
    return y + h*f(y,t)

def step_middle_point(y,t,h,f):
    """Paramètres :
    y : le résultat de l'estimation du pas précédent
    t : le point où est évaluée y
    h : la valeur du pas
    f : la fonction de l'équa diff
    Retourne l'estimation de la fonction f pour le pas suivant"""
    yn2 = y + (h/2)*f(y,t)
    pn = f(yn2,t+(h/2))
    return y + h*pn

def step_Heun(y,t,h,f):
    """Paramètres :
    y : le résultat de l'estimation du pas précédent
    t : le point où est évaluée y
    h : la valeur du pas
    f : la fonction de l'équa diff
    Retourne l'estimation de la fonction f pour le pas suivant"""
    pn1 = f(y,t)
    yn1 = y + h * pn1
    pn2 = f(yn1, t + h)
    return y + (h/2) * (pn1 + pn2)

def step_Runge_Kutta_4(y,t,h,f):
    """Paramètres :
    y : le résultat de l'estimation du pas précédent
    t : le point où est évaluée y
    h : la valeur du pas
    f : la fonction de l'équa diff
    Retourne l'estimation de la fonction f pour le pas suivant"""
    pn1 = f(y,t)
    yn1 = y + (h/2) * pn1
    pn2 = f(yn1, t + (h/2))
    yn2 = y + h * pn2 / 2
    pn3 = f(yn2, t + (h/2))
    yn3 = y + h * pn3
    pn4 = f(yn3, t + h)
    return y + (h/6) * (pn1 + 2*pn2 + 2*pn3 + pn4)
    
def meth_n_step(y0,t0,N,h,f,meth):
    """Paramètres :
    y0 : la valeur de la condition initiale de l'équa diff
    t0 : le point où est évaluée y0
    h : la valeur du pas
    N : le nombre de valeurs à estimer
    f : la fonction de l'équa diff
    meth : la méthode utilisée : Euler, point du milieu, Heun ou RK4
    Retourne un tableau contenant l'estimation des N pas calculés de la fonction f"""
    y = []
    y.append(y0)
    t = t0
    i = 1
    while( i < N ):
        t = t + h
        ynew = meth(y[i-1],t,h,f)
        y.append(ynew)
        i = i + 1
    return y    

def meth_epsilon(y0,t0,tf,eps,f,meth):
    """Paramètres :
    y0 : la valeur de la condition initiale de l'équa diff
    t0 : le point où est évaluée y0
    tf : la dernier point d'évaluation de l'équa diff
    eps : condition d'arrêt, précision du calcul
    f : équa diff
    meth : méthode utilisée : Euler, point du milieu, Heun ou RK4
    Retourne l'estimation de la fonction solution de l'équa diff sur tous les points entre t0 et tf, espacés par un pas h calculé"""
    N = 1
    h = float(tf - t0)/float(N)
    yold = meth_n_step(y0,t0,N,h,f,meth)
    while( True ):
        N = N * 2
        h = h / 2
        y = meth_n_step(y0,t0,N,h,f,meth)
        #Recherche de l'écart max entre les valeurs obtenues et valeurs précédentes => condition d'arrêt
        size = np.shape(yold)[0]
        if( size > 1 ):
            maxi = yold[1] - y[2]
            for i in range(2,size):
                if( maxi < (yold[i] - y[2*i]) ):
                    maxi = yold[i] - y[2*i]
            if( maxi < eps ):
                return np.array([y,h,N])
        #------------------------------------
        yold = y

def meth_epsilon_dim(y0,t0,tf,eps,f,meth):
    """Paramètres :
    y0 : la valeur de la condition initiale de l'équa diff
    t0 : le point où est évaluée y0
    tf : la dernier point d'évaluation de l'équa diff
    eps : condition d'arrêt, précision du calcul
    f : équa diff
    meth : méthode utilisée : Euler, point du milieu, Heun ou RK4
    Retourne l'estimation de la fonction solution de l'équa diff sur tous les points entre t0 et tf, espacés par un pas h calculé"""
    N = 1
    h = float(tf - t0)/float(N)
    yold = meth_n_step(y0,t0,N,h,f,meth)
    while( True ):
        N = N * 2
        h = h / 2
        y = meth_n_step(y0,t0,N,h,f,meth)
        #Recherche de l'écart max entre les valeurs obtenues et valeurs précédentes => condition d'arrêt
        size = np.shape(yold)[0]
        if( size > 1 ):
            maxi = np.linalg.norm(yold[1] - y[2])
            for i in range(2,size):
                if( maxi < np.linalg.norm(yold[i] - y[2*i]) ):
                    maxi = np.linalg.norm(yold[i] - y[2*i])
                    print maxi
            if( maxi < eps ):
                return np.array([y,h,N])
        #------------------------------------
        yold = y

def tangentes(min, max, preci, f):
    X = np.arange(min,max,preci)
    for i in X:
        for j in X:
            plt.quiver(i,j,1,f(j,i))
    plt.savefig('Test')
    plt.show()

#------------------------------------
#      Tests Partie I
#------------------------------------
#TESTS Dimension 1
y0 = 1.0
y1 = lambda y,t: y/(1.0 + t*t)
#La solution exacte est y(t)=exp(arctan(t))

#Tracé des différentes courbes obtenues
t0 = 0
N = 10
preci = 0.5

A = meth_n_step(y0,t0,N,preci,y1,step_Euler)
B = meth_n_step(y0,t0,N,preci,y1,step_middle_point)
C = meth_n_step(y0,t0,N,preci,y1,step_Runge_Kutta_4)
D = meth_n_step(y0,t0,N,preci,y1,step_Heun)

f = lambda x: np.exp(np.arctan(x))
E = []
for i in range(0,N):
    E.append(f(t0+i*preci))

X = np.arange(t0,N,1)

plt.plot(X,A,'r',label='Euler')
plt.plot(X,B,'g',label='Point du milieu')
plt.plot(X,C,'b',label='Runge Kutta 4')
plt.plot(X,D,'y',label='Heun')
plt.plot(X,E,'c',label='Solution')

plt.xlabel('x')
plt.ylabel('Estimation de f(x)')
plt.legend()
plt.title('Comparaison des differentes methodes de resolution d\'equa diff')
plt.grid(True)
plt.savefig('Test2')
plt.show()

mini = -10
maxi = 10
preci = 1
tangentes(mini,maxi,preci,y1)

#Tests Dimension 2
y0 = np.mat('[1. ; 0.]')
y1 = lambda y,t: np.mat('[-y[1](t) ; y[0](t)]')

#Tracé des cercles obtenus
t0 = 0
N = 5
preci = 0.02
A = meth_n_step(y0,t0,N,preci,y1,step_Euler)
B = meth_n_step(y0,t0,N,preci,y1,step_middle_point)
C = meth_n_step(y0,t0,N,preci,y1,step_Runge_Kutta_4)
D = meth_n_step(y0,t0,N,preci,y1,step_Heun)
#Résultats identiques => bug. Vient de la définition du problème de Cauchy.
#Les fonctions ne sont pas clairement spécifiées (fonction de y,t, comment manipuler y et t), d'où le terme constant
