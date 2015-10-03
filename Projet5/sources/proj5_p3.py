import proj5 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import pdb

def Find_Min (f,t=1000):
    """Permet de trouver le min d'une fonction dans l'intervalle [0,1]"""
    x=0
    h=f(x)
    for x in range (1,t+1):
        tmp=f(float(x)/t)
        if ( tmp < h):
            h=tmp
    return h

#Test Find_Min
Find_Min(lambda x : (x-0.25)**2)
    
def Find_Max ( f , t=1000) :
    """Permet de trouver le max d'une fonction dans l'intervalle [0,1]"""
    x=0
    h=f(x)
    for x in range (1,t+1):
        tmp=f(float(x)/t)
        if ( tmp > h):
            h=tmp
    return h

 
#Test Find_Min
Find_Max(np.cos)


def define_curves ( airfoil_up ,airfoil_down , lambd=40  ) :
    """ Permet de renvoyer une liste de deux tableaux contenant respectivement les lambda courbes d'air  au-dessus(t_up) et au-dessous(t_down) de l'aile d'avion """
    t_up=range(lambd+1)
    t_down=range(lambd+1)
    hmin=Find_Min( airfoil_down)
    hmax=Find_Max (airfoil_up )
    lamb= float(1)/lambd
    t_up= [(lambda y : (lambda x :(1- y*lamb ) * airfoil_up(x) + 3 * y*lamb * hmax))(i) for i in range(lambd+1)] #liste de fonction (pointeur de pointeur)
    t_down= [(lambda y : (lambda x :(1-y*lamb) * airfoil_down(x) + 3 * y*lamb * hmin))(i) for i in range(lambd+1)]
    return (t_up,t_down)



def plot_curves(t_up,t_down,msg='Courants d\'air autour de l\'aile',lambd=40,t=1000):
   "" "Permet de tracer les courants d'air autour de l'aile mise en entree"""
   x_up_p=[float(i) for i in range(t+1)]
   x_down_p=[float(i) for i in range(t+1)]
   x_down=[[]]
   x_up=[[]]
   x=np.linspace(0,1,t+1)
   for l in range (0,lambd):
       for j in range (0 , t+1):
           x_up_p[j]=t_up[l](x[j])
           x_down_p[j]=t_down[l](x[j])
       x_up.append(x_up_p)
       x_up_p=[]
       x_up_p=[float(i) for i in range(t+1)]
       x_down.append(x_down_p)
       x_down_p=[]
       x_down_p=[float(i) for i in range(t+1)]
   plt.plot(x_up[1],'r') #pour avoir un effet mettant en valeur l'ecoulement dans l'ail 
   plt.plot(x_down[1],'r')
   for l in range (1,lambd+1):
       plt.plot(x_up[l],'k')
   for l in range (1,lambd+1):
       plt.plot(x_down[l],'k')
   plt.ylabel('hauteur de l\'aile')
   plt.title(msg)
   plt.show()
   return (x_up[1],x)
    
#Test define_curves by using y: x=>x(x-1)
k_up,k_down=define_curves((lambda x: -x*(x-1)) , (lambda x: x*(x-1) ))
l,x=plot_curves(k_up,k_down,'Test en utilisant y: x=>x(x-1)')

#Test define_curves en utilisant les deux fonction de l'ailfoil
ye = lambda x: splint(ex,ey,y2e,ne,x)
yi = lambda x: splint(ix,iy,y2i,ni,x)

k_up,k_down=define_curves(ye ,yi)
l,x=plot_curves(k_up,k_down,'Courants d\'air autour de l\'aile')


def static_pressure ( P , distance , time ,r ):
    """Calcule la pression static \(non utilise dans l'algo\)"""
    V=distance/time
    return P-(0.5 * r * V*V)

def speed(distance,time):
    """calcul la vitesse"""
    return float(distance)/time

def dynamic_pressure(distance, time,r):
    """permet de trouver la pression dynamic tel que r est la densite de l'air"""
    V=speed(distance,time);
    return 0.5*r*V*V

        

def appropriate_length(T_up,T_down,y,x,lambd=40):
    """permet de renvoyer la longeur de la vague de l'air.
    Parametres: T_up la liste de fonctions des courants de l'air au-dessus de l'aile
                T_down celles en bas de l'aile
                y la coordonee dans l'abscisse y du point dans le plan
                x la coordonee dans l'abscisse x du point dans le plan
                lamnbda: le nombre de fonctions contenues dans les T_up et T_down"""
    if y>=T_up[0](x) :
        i=0
        while y>T_up[i](x) and i<lambd :
            i=i+1
        return length(0,1,Simpson,0.7,50,T_up[i-1])[0]
    elif y<=T_down[0](x) :
        i=0
        while y<T_down[i](x) and i<lambd :
            i=i+1
        return length(0,1,Simpson,0.7,50,T_down[i-1])[0]
    else :           #si on est dans l'aile
        return 0.99  #meilleur compromis trouve pour avoir un bon graph de pression
    

def map_pressure(T_up,T_down,time=10,r=0.5,t=200):
    """permet de generer la map de pression """
    y=np.linspace(-0.5,0.5,t) #abscisse y
    x=np.linspace(0,1,t)      #abscisse x
    xmax=1                    #/* les limites des abscisses 
    ymax=0.25                 # 
    ymin=-0.25                #
    xmin=0                    #*/
    T = np.zeros([ t, t ])    #on initialise notre matrice t*t a 0
    for j in range(0,t):  
        print "Generation de map pressure:",(float(j)/t)*100,"%"
        for i in range(0,t):
            T[i][j]=dynamic_pressure(appropriate_length(T_up,T_down,y[j],x[i]),time,r) #on remplie notre matrice
    cmap=cm.get_cmap('afmhot') #on choisit le degrade voulut 
    pd=plt.pcolor(x,y,T.T,cmap=cmap)   #equivalent a plot ou contourf ou contour
    cd=plt.colorbar(pd,orientation='horizontal') #modification de l'emplacement du colrbar
    cd.set_label('Dynamic Pressure') #permet de modifier le label de colorbar
    plt.xlim(xmin,xmax)       #mis en place des limites en x
    plt.ylim(ymin,ymax)       #mis en place des limites en y
    plt.show()
    print "----------------------------------------------------------------- "
    print "!!!!!!!!!!!!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
    print "Pour generer une map avec la qualite utilisee dans le rapport,  il"
    print "faut mettre comme en parametre t=1000 mais cela va prendre 30 min "
    print "il faut aussi changer la couleur (cm.get_cmap()) en 'jet'."
    return T

#Application de la map_pressure  
T= map_pressure(k_up,k_down)
