# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import load_foil as foil
from scipy import misc

#---------------------------------------------------
#---------------Airfoil refinment-------------------
#---------------------------------------------------
def spline(x, y, n, yp1, ypn):
    """Function from Numerical Recipes
    Parameters:
    -x is an array containing abscissae (of the wing in our case)
    -y is an array containing f(x) values y[i] = f(x[i])
    -n is the size of x (and so of y)
    -yp1 is the first derivative of f at x0 (f'(x[0]))
    -ypn is the first derivative of f at xn-1 (f'(x[n-1]))
    Returns: y2 an array containing the second derivatices of the interpolating function at the tabulated points xi
    """
    u = []
    y2 = []
    for i in range(0,n):
        u.append(0.0)
        y2.append(0.0)
    if( yp1 > 0.99e30 ):
        u[0] = 0.0
        y2[0] = 0.0
    else:
        y2[0] = -0.5
        u[0] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)
    for i in range(1,n-1): #Donc de 1 à n-2 inclus
        sig = (x[i]-x[i-1])/(x[i+1]-x[i-1])
        p = sig * y2[i-1]+2.0
        y2[i] = (sig-1.0)/p
        u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
        u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p
    if( ypn > 0.99e30 ):
        qn = 0.0
        un = 0.0
    else:
        qn = 0.5
        un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
    y2[n-1] = (un - qn*u[n-2])/(qn*y2[n-2]+1.0)
    for i in range(n-2,0): #De n-2 à 1 inclus
        y2[i] = y2[i]*y2[i+1] + u[k]
    return y2

def compute_all_derivate(xa, y, y2):
    """Parameters:
    -x is an array containing abscissae (of the wing in our case)
    -y is an array containing f(x) values y[i] = f(x[i])
    -y2 is an array containing f''(x) values y2[i] = f''(x[i])
    Returns deriv an array containing the derivatives of the function f on each interval"""
    n = xa.shape[0]
    deriv = []
    for i in range(0,n):
        deriv.append(lambda x: 0)
    for j in range(0,n-1): #Le bug de dérivée est ici
        deriv[j] = lambda x: (y[j+1] - y[j])/(xa[j+1] - xa[j]) - ((3*((xa[j+1] - x)/(xa[j+1] - xa[j]))**2 - 1)/6)*(xa[j+1] - xa[j])*y2[j] + ((3* ((x - xa[j])/(xa[j+1] - xa[j]))**2 - 1)/6)*(xa[j+1] - xa[j])*y2[j+1]
    return deriv

def derivate(xa, y1, n, x):
    """Parameters:
    -x is an array containing abscissae (of the wing in our case)
    -y1 is an array containing the derivatives of the function f on each interval
    -n is the size of the array x
    -x is the point where the derivate is calculated
    Returns the derivate function of f calculated at x"""
    h = 0
    for j in range(0,n-1):
        #print y1[j](0.25) #Notre dérivée est constante sur chaque intervalle
        if( (xa[j] < x) and (xa[j+1] >= x) ):
            h = j
            break
    return y1[h](x)

def splint(xa, ya, y2a, n, x):
    """Function from Numerical Recipes
    Parameters:
    -xa is an array containing abscissae (of the wing in our case)
    -ya is an array containing f(xa) values ya[i] = f(xa[i])
    -y2a is an array containing f''(xa) values (from the function spline above) y2a[i] = f''(xa[i])
    -n is the size of these arrays
    -x is the point where to interpolate f
    Returns: a cubic-spline interpolated value y = f(x)
    """
    klo = 0
    khi = n-1
    while(khi-klo > 1):
        k = (khi+klo)/2
        if( xa[k] > x ):
            khi = k
        else:
            klo = k
    h = xa[khi] - xa[klo]
    if( h == 0.0 ):
        print("Error.")
    a = (xa[khi]-x)/h
    b = (x-xa[klo])/h
    y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0
    return y

def draw_splints(x1,f1,n1,x2,f2,n2,precision):
    """Parameters:
    -x1 is an array containing abscissae (a part of the wing in our case)
    -f1 is the function interpolated to be drawn
    -n1 is the size of x1
    -x2, f2 and n2 are exactly the same for another function/abscissa
    -precision is related to the number of points drawn
    Draw the two functions
    """
    X1 = []
    Y1 = []
    X2 = []
    Y2 = []
    for i in range(0,n1-1):
        for j in range(0,precision):
            xi = x1[i]+float(j)/precision
            if( xi < x1[i+1] ):
                X1.append(xi)
                Y1.append(f1(xi))
    for i in range(0,n2-1):
        for j in range(0,precision):
            xi = x2[i]+float(j)/precision
            if( xi < x2[i+1] ):
                X2.append(xi)
                Y2.append(f2(xi))

    X1.append(x1[n1-1])
    Y1.append(f1(x1[n1-1]))
    X2.append(x2[n2-1])
    Y2.append(f2(x2[n2-1]))
    plt.plot(X1,Y1)
    plt.plot(X2,Y2)
    plt.title("Wing curve goe281 with spline")
    #plt.axis([-0.1,1.1,-0.1,1.1])
    plt.savefig("goe281_spline_preci_1000")
    plt.show()

def relative_err(x, y, n, f):
    """Parameters:
    -x is an array containing abscissae
    -y is an array containing the values of g(x). y[i] = g(x[i])
    -f is a cubic-spline interpolating g
    -n is the length of x
    Returns average of relatives error between f and y 
    """
    cpt = 0.0
    for i in range(0,n):
        if( y[i] > 1e-6 ):
            cpt = cpt + (f(x[i]) - y[i])/y[i]
    cpt = cpt / n
    return cpt
#---------------------------------------------------
#-------Airfoil refinment, application--------------
#---------------------------------------------------
(ex,ey,ix,iy) = foil.load_foil("goe281.dat")
plt.plot(ex,ey)
plt.plot(ix,iy)
plt.title("Wing curve goe281 without spline")
#plt.axis([-0.1,1.1,-0.1,1.1])
plt.savefig("goe281_no_spline")
plt.show()

#Taille des tableaux de points de l'aile
#e correspond à extrados
#i correspond à intrados
ne = ex.shape[0]
ni = ix.shape[0]

#On prend les dérivées premières comme nulles
yp1 = 0.0
ypn = 0.0

#Calcul du tableau des dérivées secondes aux points de l'aile
y2e = spline(ex,ey,ne,yp1,ypn)
y2i = spline(ix,iy,ni,yp1,ypn)

#Calcul des fonctions de la courbe de l'aile
ye = lambda x: splint(ex,ey,y2e,ne,x)
yi = lambda x: splint(ix,iy,y2i,ni,x)

#calcul du tableau des dérivées premières sur les intervalles de l'aile
y1e = compute_all_derivate(ex, ey, y2e)
y1i = compute_all_derivate(ix, iy, y2i)

#calcul des fonctions dérivées de la courbe de l'aile
ype = lambda x: derivate(ex, y1e, ne, x)
ypi = lambda x: derivate(ix, y1i, ni, x)

#Tests de la dérivée
print "-----Tests de la dérivée------"
print ype(0.9)
print ype(0.8)
print ypi(0.01)
print ypi(0.09)

#ye(0.0248400) #Quelques tests pour vérifier la validité de nos résultats/les chiffres significatifs et sert d'exemple d'utilisation
#yi(0.1500000)

precision = 1000 #Précision 1 correspond au recalcul des points donnés.
#A partir de 100, plus d'effet visuellement

draw_splints(ex,ye,ne,ix,yi,ni,precision)

err_ye = relative_err(ex,ey,ne,ye)
err_yi = relative_err(ix,iy,ni,yi)
#---------------------------------------------------
#-----Computing the length of plane curves----------
#---------------------------------------------------

def left_rect(a, b, f, h):
    """Parameters:
    -a,b define the interval on which f is evaluated
    -f is a function
    -h is the step for integration
    Returns left rectangle integral of f on [a,b] calculated with a step h
    """
    Irec = 0
    n = int((b-a)/h)
    for i in range(0,n):
        Irec = Irec + f(a + i*h)
    Irec = h*Irec
    return Irec

def right_rect(a, b, f, h):
    """Parameters:
    -a,b define the interval on which f is evaluated
    -f is a function
    -h is the step for integration
    Returns right rectangle integral of f on [a,b] calculated with a step h
    """
    Irec = 0
    n = int((b-a)/h)
    for i in range(1,n+1):
        Irec = Irec + f(a + i*h)
    Irec = h*Irec
    return Irec

def middle_point(a, b, f, h):
    """Parameters:
    -a,b define the interval on which f is evaluated
    -f is a function
    -h is the step for integration
    Returns middle point integral of f on [a,b] calculated with a step h
    """
    Imid = 0
    n = int((b-a)/h)
    for i in range(0,n):
        Imid = Imid + f(a+i*h+h/2)
    Imid = h*Imid
    return Imid

def trap(a, b, f, h):
    """Parameters:
    -a,b define the interval on which f is evaluated
    -f is a function
    -h is the step for integration
    Returns trapezic integral of f on [a,b] calculated with a step h
    """
    Itrap = 0
    n = int((b-a)/h)
    for i in range(0,n):
        Itrap = Itrap + f(a+i*h)
    Itrap = h*( (f(a)+f(b))/2 + Itrap )
    return Itrap

def Simpson(a, b, f, h):
    """Parameters:
    -a,b define the interval on which f is evaluated
    -f is a function
    -h is the step for integration
    Returns integral of f on [a,b] calculated with a step h calculated with the Simpson method
    """
    Isim1 = 0
    Isim2 = 0
    i = 1
    n = int((b-a)/h)
    while( i <= n/2 ):
        if( i <= (n/2 -1) ):
            Isim1 = Isim1 + f(a+2*i*h)
        Isim2 = Isim2 + f(a+(2*i - 1)*h)        
        i = i + 1
    Isim = (h/3)*(f(a) + 2*Isim1 + 4*Isim2 + f(b))
    return Isim

def apply_integrate(a, b, f, h, integrate):
    """Allows to choose the method for integration of f
    -integrate : Simpson, trap, middle point, left_rec or right rec
    """
    return integrate(a,b,f,h)

def integrate(a, b, f, integrate_method, epsilon, maxiter):
    h = float(b-a)/float(10)
    i = 0
    Sold = apply_integrate(a, b, f, h, integrate_method)
    while(i < maxiter):
        h = h/2
        Snew = apply_integrate(a, b, f, h, integrate_method) #Pas optimisé. On recalcule tous les termes à chaque fois. Il est possible de n'en recalculer que quelques-uns (voir cours)
        if( abs(Snew-Sold) < epsilon ):
            return np.array([Snew,float((b-a)/h)])
        Sold = Snew
        i = i + 1
    return np.array([Snew,float((b-a)/h)])
#Faire une fonction opti une fois la meilleure façon d'intégrer choisie ?

def length(a, b, integrate_method, epsilon, maxiter, spline_derivate):
    """Compute the length of SPLINE between A and B, using INTEGRATE_METHOD function to integrate"""
    return integrate(a,b,lambda x : np.sqrt(1 + misc.derivative(spline_derivate,x,dx=1e-6)**2), integrate_method,epsilon,maxiter)
#On utilise misc.derivative comme notre fonction de dérivation ne fonctionne pas

def integral_conv(a, b, abscissae, f, method):
    N = np.shape(abscissae)[0]
    Y = []
    for i in range(0,N):
        Y.append(method(a,b,f,float((b-a))/float(abscissae[i])))
    return Y

#---------------------------------------------------
#-Computing the length of plane curves, application-
#---------------------------------------------------

#TESTS de l'intégration
print "Test de la méthode d'intégration (pas encore optimisée) || Nombre d'itérations\n"
print "-----------Partie haute de l'aile------------------"
print "Rectangle gauche"
print integrate(0,1,ye,left_rect, 10e-8,1000)
print "Rectangle droit"
print integrate(0,1,ye,right_rect, 10e-8,1000)
print "Point milieu"
print integrate(0,1,ye,middle_point, 10e-8,1000)
print "Trapèze"
print integrate(0,1,ye,trap, 10e-8,1000)
print "Simpson"
print integrate(0,1,ye,Simpson, 10e-8,1000)
print "-----------Partie basse de l'aile------------------"
print "Rectangle gauche"
print integrate(0,1,yi,left_rect, 10e-8,1000)
print "Rectangle droit"
print integrate(0,1,yi,right_rect, 10e-8,1000)
print "Point milieu"
print integrate(0,1,yi,middle_point, 10e-8,1000)
print "Trapèze"
print integrate(0,1,yi,trap, 10e-8,1000)
print "Simpson"
print integrate(0,1,yi,Simpson, 10e-8,1000)

#TESTS Convergence des différentes méthodes d'intégrales
#Partie haute de l'aile, peu d'itérations
Xe = np.array([10,11,12,13,14,15,16,17,18,19,20,40,80,100])
Ye_rect = integral_conv(0,1,Xe,ye,left_rect)
Ye_mid = integral_conv(0,1,Xe,ye,middle_point)
Ye_trap = integral_conv(0,1,Xe,ye,trap)
Ye_Simpson = integral_conv(0,1,Xe,ye,Simpson)
plt.plot(Xe,Ye_rect,'r',label='Rectangle')
plt.plot(Xe,Ye_mid,'g',label='Point milieu')
plt.plot(Xe,Ye_trap,'b',label='Trapezes')
plt.plot(Xe,Ye_Simpson,'y',label='Simpson')
plt.xlabel('Nombre iterations (echelle log)')
plt.ylabel('Valeur integrale')
plt.legend()
plt.title("Convergence des differentes methodes d'integration Extrados")
#plt.axis([0,200,0.058,0.060])
plt.xscale('log')
plt.grid(True)
plt.savefig("integration_convergence_ex_few_it")
plt.show()

#Partie haute de l'aile, beaucoup d'itérations
Xe = np.array([1000,2000,3000,4000,5000,6000,7000,7100,7200,7300,7400,7500,7600,7700,7800,7900,8000,8001,8002,8003,8004,8005,8006,8007,8008,8009,8100,8200,8300,8400,8500,8600,8700,8800,8900,9000,10000])
Ye_rect = integral_conv(0,1,Xe,ye,left_rect)
Ye_mid = integral_conv(0,1,Xe,ye,middle_point)
Ye_trap = integral_conv(0,1,Xe,ye,trap)
Ye_Simpson = integral_conv(0,1,Xe,ye,Simpson)
plt.plot(Xe,Ye_rect,'r',label='Rectangle')
plt.plot(Xe,Ye_mid,'g',label='Point milieu')
plt.plot(Xe,Ye_trap,'b',label='Trapezes')
plt.plot(Xe,Ye_Simpson,'y',label='Simpson')
plt.xlabel('Nombre iterations (echelle log)')
plt.ylabel('Valeur integrale')
plt.legend()
plt.title("Convergence des differentes methodes d'integration Extrados")
#plt.axis([0,200,0.058,0.060])
plt.xscale('log')
plt.grid(True)
plt.savefig("integration_convergence_ex_lot_it")
plt.show()

#Partie basse de l'aile peu d'itérations
Xi = np.array([10,11,12,13,14,15,16,17,18,19,20,40,80,100])
Yi_rect = integral_conv(0,1,Xi,yi,left_rect)
Yi_mid = integral_conv(0,1,Xi,yi,middle_point)
Yi_trap = integral_conv(0,1,Xi,yi,trap)
Yi_Simpson = integral_conv(0,1,Xi,yi,Simpson)
plt.plot(Xi,Yi_rect,'r',label='Rectangle')
plt.plot(Xi,Yi_mid,'g',label='Point milieu')
plt.plot(Xi,Yi_trap,'b',label='Trapezes')
plt.plot(Xi,Yi_Simpson,'y',label='Simpson')
plt.xlabel('Nombre iterations (echelle log)')
plt.ylabel('Valeur integrale')
plt.legend()
plt.title("Convergence des differentes methodes d'integration Intrados")
#plt.axis([0,200,0.058,0.060])
plt.xscale('log')
plt.grid(True)
plt.savefig("integration_convergence_in_few_it")
plt.show()

#Partie basse de l'aile beaucoup d'itérations
Xi = np.array([1000,2000,3000,4000,5000,6000,7000,7100,7200,7300,7400,7500,7600,7700,7800,7900,8000,8001,8002,8003,8004,8005,8006,8007,8008,8009,8100,8200,8300,8400,8500,8600,8700,8800,8900,9000,10000])
Yi_rect = integral_conv(0,1,Xi,yi,left_rect)
Yi_mid = integral_conv(0,1,Xi,yi,middle_point)
Yi_trap = integral_conv(0,1,Xi,yi,trap)
Yi_Simpson = integral_conv(0,1,Xi,yi,Simpson)
plt.plot(Xi,Yi_rect,'r',label='Rectangle')
plt.plot(Xi,Yi_mid,'g',label='Point milieu')
plt.plot(Xi,Yi_trap,'b',label='Trapezes')
plt.plot(Xi,Yi_Simpson,'y',label='Simpson')
plt.xlabel('Nombre iterations (echelle log)')
plt.ylabel('Valeur integrale')
plt.legend()
plt.title("Convergence des differentes methodes d'integration Intrados")
#plt.axis([0,200,0.058,0.060])
plt.xscale('log')
plt.grid(True)
plt.savefig("integration_convergence_lot_it")
plt.show()


#TESTS longueur aile
print "----------------------------"
print "Longueur de l'aile avec Simpson"
print "Partie haute"
print length(0,1,Simpson,10e-6,1000,ye)
print "Partie basse"
print length(0,1,Simpson,10e-6,1000,yi)
print "Longueur de l'aile avec point du milieu"
print "Partie haute"
print length(0,1,middle_point,10e-6,1000,ye)
print "Partie basse"
print length(0,1,middle_point,10e-6,1000,yi)
