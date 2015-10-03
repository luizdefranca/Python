import numpy as np
import matplotlib.pyplot as plt
import pdb
import matplotlib.cm as cm

def frequence_borda(Theta_O,l=5,g=10):
    """Calcul la frequence pour n'importe quel amplitude"""
    wO=np.sqrt((float(l)/g))
    return 1/float(wO*(1+float(Theta_O**2)/16))


def plot_frequence(l=5,g=10,t=300):
    O=np.linspace(-np.pi,np.pi,t)
    w=np.linspace(-np.pi,np.pi,t)   
    for i in range(0,t):
        w[i]=frequence_borda(O[i])
    plt.plot(O,w,'g',label='$\Theta_0$')
    C=np.linspace(np.sqrt(float(g)/l),np.sqrt(float(g)/l),t)
    plt.plot(O,C, 'r',label='$\sqrt{l/g}$')
    plt.xlabel('frequence')
    plt.ylabel('$\Theta_0$ en radian')
    plt.legend(loc='best')
    plt.title('Comparaison entre les frequences et $\sqrt{l}/g$')
    plt.savefig('un_maillon')
    plt.show()


#Recherche du graphe de deux pendules
m1 = 42
m2 = 61
L1 = 1002
L2 = 135
g = 9.81

# y0 = lambda t :np.array([0. , 0. , 0. , 0.])

# y1 = lambda t,y: np.array([y(t)[2], y(t)[3] , (-g*(2.0*m1+m2)*np.sin(y(t)[0]) - m2*g*np.sin(y(t)[0] - 2.0*y(t)[1]) - 2.0*np.sin(y(t)[0]-y(t)[1])*m2*(y(t)[3]*y(t)[3]*L2 + y(t)[2]*y(t)[2]*L1*np.cos(y(t)[0]-y(t)[1]))) / (L1*(2.0*m1 + m2 - m2 * np.cos(2.0*y(t)[0] - 2.0*y(t)[1]))) ,(2.0*np.sin(y(t)[0] - y(t)[1])*(y(t)[2]*y(t)[2]*L1*(m1+m2)+g*(m1+m2)*np.cos(y(t)[0]) + y(t)[3]*y(t)[3]*L2*m2*np.cos(y(t)[0]-y(t)[1]))) / (L2*(2.0*m1+m2-m2*np.cos(2.0*y(t)[0]-2.0*y(t)[1])))])

y1 = lambda t,y: np.array([y[2], y[3] , (-g*(2.0*m1+m2)*np.sin(y[0]) - m2*g*np.sin(y[0] - 2.0*y[1]) - 2.0*np.sin(y[0]-y[1])*m2*(y[3]*y[3]*L2 + y[2]*y [2]*L1*np.cos(y[0]-y[1]))) / (L1*(2.0*m1 + m2 - m2 * np.cos(2.0*y[0] - 2.0*y[1]))) ,(2.0*np.sin(y[0] - y[1])*(y[2]*y[2]*L1*(m1+m2)+g*(m1+m2)*np.cos(y[0]) + y[3]*y[3]*L2*m2*np.cos(y[0]-y[1]))) / (L2*(2.0*m1+m2-m2*np.cos(2.0*y[0]-2.0*y[1])))])

y0 = np.array([(np.pi-0.25)/2, (np.pi-0.25)/2,0.,0.])

def rk4(x, h, y, f):
    #pdb.set_trace()
    k1 = h * f(x, y)
    k2 = h * f(x + 0.5*h, y + 0.5*k1)
    k3 = h * f(x + 0.5*h, y + 0.5*k2)
    k4 = h * f(x + h, y + k3)
    return x + h, y + (k1 + 2*(k2 + k3) + k4)/6.0

def rk4_main(t,dt,y0,y1):
    # y0=y0(t)
    theta1=[]
    theta2=[]
    while t < 500:
        t, y0 = rk4( t, dt , y0, y1)
        theta1.append(y0[0])
        theta2.append(y0[1])
        print t,y0
    return theta1,theta2

#theta1,theta2=rk4_main(0,0.02,y0,y1)

def couple_1(theta1):
    return  L1*np.sin(theta1),-L1*np.cos(theta1) 
def couple_2(theta2,x1,y1):
    return  x1+L1*np.sin(theta2),y1-L1*np.cos(theta2) 

def plot_traj_deuxpendule(theta1,theta2):
    x1,y1=couple_1(theta1)
    x2,y2=couple_2(theta2,x1,y1)
    plt.plot(x1,y1,'b',label="trajectoire du premier maillon")
    plt.plot(x2,y2,'r',label="trajectoire du deuxieme maillon")
    #plt.plot(theta1,theta2)
    plt.legend(loc='best')
    plt.title('trajectoire pendant 500s de deux maillons avec $\Theta1$,$\Theta2$$=((\pi-0.25)/2,(\pi-0.25)/2$')
    plt.savefig('deux_maillon_p1')
    plt.show()

plot_traj_deuxpendule(theta1,theta2)


def time_reverse(t,dt,y0,y1,limite=500):
    #pdb.set_trace()
    # t,y0=rk4( t, dt , y0, y1)
    # x_1,y_1=couple_1(y0[0])
    # x_2,y_2=couple_2(y0[1],x_1,y_1)
    # x_prec=x_1
    # x_nv=x_1
    # while(x_nv<=x_prec):
    #     x_prec=x_nv
    #     t, y0 = rk4( t, dt , y0, y1)
    #     x_1,y_1=couple_1(y0[0])
    #     x_2,y_2=couple_2(y0[1],x_1,y_1)
    #     x_nv=x_1
    #     if (t>limite):
    #         return t
    # return t
    while(1):
        if y0[3]<=0:
            t, y0 = rk4( t, dt , y0, y1)
            if (y0[1]<-np.pi):
                return t 
        else:
            t, y0 = rk4( t, dt , y0, y1)
            if (y0[1]>np.pi):
                return t
        if (t>limite):
            return t
            


def map_time_reverse(y1,t=100):
    y0=np.array([0.,0.,0.,0.])
    y=np.linspace(-3,3,t) #abscisse y
    x=np.linspace(-3,3,t)      #abscisse x
    xmax=3
    xmin=-3
    ymax=3
    ymin=-3
    T = np.zeros([ t, t ])
    for j in range(0,t):
         print "Generation de map time_reverse:",(float(j)/t)*100,"%"
         for i in range(0,t):
             y0[0],y0[1]=x[i],y[j]
             T[i][j]=time_reverse(0,1,y0,y1)
    cmap=cm.get_cmap('jet') #on choisit le degrade voulut 
    pd=plt.pcolor(x,y,T.T,cmap=cmap)   #equivalent a plot ou contourf ou contour
    cd=plt.colorbar(pd,orientation='vertical') #modification de l'emplacement du colrbar
    cd.set_label('temps de renversement') #permet de modifier le label de colorbar
    plt.xlim(xmin,xmax)       #mis en place des limites en x
    plt.ylim(ymin,ymax)       #mis en place des limites en y
    plt.show()

map_time_reverse(y1)
