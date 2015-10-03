import numpy as np
import matplotlib.pyplot as mp
import random

print "Projet 1"

# FUNCTIONS RP _____________________________________________________________________________

# def order(x): #np.rint
#     # """ The function "order" takes one argument : 
#     # x is a real number    
#     # and returns the power of ten of the first significant number """
#     if( x < 1 ):
#         i = 1
#         x = x * 10.0
#         while( x < 1 ):
#             x = x * 10.0
#             i = i + 1
#         return -i+1
#     else:
#         i = 1
#         x = x / 10.0
#         while( x > 1 ):
#             x = x / 10.0
#             i = i + 1
#         return i

def rp(x,p):
    xp = abs(x)
    sign = 1
    i = 0
    j = 0
    if (x < 0) :
        sign = -1
    while (int(xp) != 0):
        xp = xp / 10
        j = j + 1
    while ((int(xp / pow(10.0,p-1)) < 1) and (xp != 0)):
        xp = xp*10.0
        i = i + 1
    xp = int(np.rint(xp))
    xp = xp / pow(10.0,i-j)
    return sign*xp

# def rp(x,p):
#     """ The function "rp" takes two arguments :
#     x is a real number
#     p is an integer
#     and returns the rounded number x with p significative numbers """
#     n = order(x)
#     num = x
#     num = num * pow(10,-n+p+1)
#     if( str(x)[-n+p] >= 5 ):
#         num = num + 5
#     num = num / 10
#     num = np.floor(num)
#     num = num * pow(10,n-p)
#     return num

# FUNCTIONS ERREUR ________________________________________________________________________

def erreur_rsum(x, y, p):
    """Calcule l'erreur relative pour la somme"""
    tmp = (rp(x,p)+rp(y,p))
    n = abs((x+y)-tmp)
    d = abs(tmp)
    return (n/d)

def erreur_rprod(x,y,p):
    """Calcule l'erreur relative pour le produit"""
    tmp = (rp(x,p)*rp(y,p))
    n = abs((x*y)-tmp) 
    d = abs(tmp) 
    return (n/d)

def trace_rsum(x,p):
    """Trace la courbe de la fonction erreur_rsum en fonction de y"""
    y = np.arange(0, 100.0, 0.1)
    s = []
    for i in range(0,len(y)):
        s.append(erreur_rsum(x,y[i],p))
    mp.plot(y, s, linewidth=1.0)
    mp.title('Graphe erreur relative de la somme')
    mp.xlabel('la valeur que nous fesons varier')
    mp.ylabel('erreur relative')
    mp.savefig("err_rs")    
    mp.show()

def trace_rprod(x,p):
    """Trace la courbe de la fonction erreur_rprod en fonction de y """
    y = np.arange(0.01, 100.0, 0.1)
    s = []
    for i in range(0,len(y)):
        s.append(erreur_rprod(x,y[i],p))
    mp.plot(y, s, linewidth=1.0)
    mp.xlabel('la valeur que nous fesons varier')
    mp.ylabel('erreur relative')
    mp.title('Graphe erreur relative du produit')
    mp.savefig("err_rp")    
    mp.show()
    
# FUNCTIONS APPROX LN(2) ___________________________________________________________________

def ap(p):
    """ The function "ap" takes one argument : 
    p is a real number    
    and returns the value of log(2) with p significant numbers """
    n=2
    l=1
    while n<=p:        
        l=float(l)+float(pow(-1,(n+1)))/n
        n=n+1
    return l

def ap2(p):
    """ The function "ap2" takes one argument : 
    p is a real number    
    and returns the value of log(2) with p significant numbers """
    i=1
    n=2
    l=1
    b=np.log(2)
    while i<=p:
        while str(b)[i]!=str(float(l))[i]:
            l=float(l)+float(pow(-1,(n+1)))/n
            n=n+1
        i=i+1        
    return l

# ALGORITHMES CORDIC _______________________________________________________________________

"""Precalcul des valeurs pour les algorithmes de CORDIC"""
L= [np.log(2), np.log(1.1), np.log(1.01), np.log(1.001), np.log(1.0001), np.log(1.00001), np.log(1.000001)]
A= [np.arctan(1), np.arctan(0.1), np.arctan(0.01), np.arctan(0.001), np.arctan(0.0001)]


def ln(x):
    """ The function "ln" takes one argument : 
    x is a real number (algorithm is optimized for numbers in [1, 10[)
    and returns the value of ln(x) with an error of 10^(-12) """ 
    k = 0
    y = 0
    p = 1
    while( k <= 6 ):
        while( (p + p * pow(10,-k)) <= x ):
            y = y + L[k]
            p = p + p * pow(10,-k);
        k = k + 1
    return ( y + ((x / p) - 1) )

def e(x):
    """ The function "e" takes one argument : 
    x is a real number (algorithm is optimized for numbers in [0, ln(10)[)
    and returns the value of e(x) with an error of 10^(-12) """
    k = 0
    y = 1
    while( k <= 6 ):
        while( x >= L[k] ):
            x = x - L[k]
            y = y + y * pow(10,-k)
        k = k + 1
    return( y + y * x )

def arctan(x):
    """ The function "arctan" takes one argument : 
    x is a real number (algorithm is optimized for numbers in [0, 1[)
    and returns the value of arctan(x) with an error of 10^(-12) """
    k = 0
    y = 1
    r = 0
    while( k <= 4 ):
        while( x < (y * pow(10,-k)) ):
            k = k + 1
        if( k > 4 ):
            break
        xp = x - y * pow(10,-k)
        y = y + x * pow(10,-k)
        x = xp
        r = r + A[k]
    return ( r + ( x / y ) )

def tan(x):
    """ The function "tan" takes one argument : 
    x is a real number (algorithm is optimized for numbers in [0, Pi/4[)
    and returns the value of tan(x) with an error of 10^(-12) """
    k = 0
    n = 0
    d = 1
    while( k <= 4 ):
        while( x >= A[k] ):
            x = x - A[k]
            np = n + d * pow(10,-k)
            d = d - n * pow(10,-k)
            n = np
        k = k + 1
    return ( (n + x * d) / (d - x * n) )

# GRAPHES ___________________________________________________________________________________

trace_rsum(np.pi,4)

trace_rprod(np.pi,4) 

# TESTS RP _________________________________________________________________________________
       
print "\nTests de rp(x,p) :"
print "\nPour  3.141592658"
print "Avec p = 4 : "
print rp(3.141592658,4)
print "Avec p = 6 : "
print rp(3.141592658,6)
print "\nPour 10507.1823"
print "Avec p = 4 : "
print rp(10507.1823,4)
print "Avec p = 6 : "
print rp(10507.1823,6)
print "\nPour 0.0001857563"
print "Avec p = 4 : "
print rp(0.0001857563,4)
print "Avec p = 6 : "
print rp(0.0001857563,6)

# TESTS ERREURS ___________________________________________________________________________

print "\nTests de erreur_rsum(x,y,p) :"
print "\nPour x = 3.141592658 et y = 10507.1823 :"
print "Avec p = 4 :"
print erreur_rsum(3.141592658,10507.1823,4)
print "Avec p = 6 :"
print erreur_rsum(3.141592658,10507.1823,6)
print "\nTests de erreur_rprod(x,y,p) :"
print "Pour x = 0.0001857563, et y = 10507.1823 :"
print "Avec p = 4 :"
print erreur_rsum(0.0001857563,10507.1823,4)
print "Avec p = 6 :"
print erreur_rsum(0.0001857563,10507.1823,6)

# TESTS APPROX LN(2) _______________________________________________________________________

print "\nApproximations de ln(2) avec ap(x) :"
print "Approximation :"
print ap(10000)
print "Approximation sur 3 chiffres :"
print rp(ap(10000),3)
print "Le vrai ln(2) :"
print np.log(2)
print "Erreur relative entre ln(2) et son approximation :"
print float(abs(np.log(2)-ap(1000000)/np.log(2)))

# TESTS APPROX2 LN(2) ______________________________________________________________________

print "Aproximation (2) de ln(2) avec ap2(x) :"
print ap2(7)
print "Erreur relative avec ln(2) :"
print float(abs(100*np.log(2)-ap2(7)))/np.log(2)

# TESTS CORDIC LOG _________________________________________________________________________

r = random.randint(1,100)
print "\nTest CORDIC ln(r) :"
print "r = %d" % r
print ln(r)
print "Le vrai ln(r)"
print np.log(r)

# TESTS CORDIC EXPONENTIELLE _______________________________________________________________

r = random.randint(0,100)
print "\nTest CORDIC e(r) :"
print "r = %d" % r
print e(r)
print "Le vrai e(r)"
print np.exp(r)

# TESTS CORDIC ARCTAN ______________________________________________________________________

r = random.randint(0,100)
print "\nTest CORDIC arctan(r) :"
print "r = %d" % r
print arctan(r)
print "Le vrai arctan(r)"
print np.arctan(r)

# TESTS CORDIC TAN _________________________________________________________________________

r = random.randint(0,100)
print "\nTest CORDIC tan(r) :"
print "r = %d" % r
print tan(r)
print "Le vrai tan(r)"
print np.tan(r)
