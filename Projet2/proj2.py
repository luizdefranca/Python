# -*- coding: utf-8 -*-

import numpy as np
import random
import scipy.linalg   # SciPy Linear Algebra Library
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# FONCTIONS CHOLESKY _______________________________________________________________________
def mSDP (n,p):
    """ Prend deux paramètres :
    n : la taille de la matrice carrée de sortie
    p : un entier pour borner les coefficients aléatoires de la matrice de sortie
    (la loi de probabilité est une loi uniforme entre 0 et p inclus
    Retourne une matrice symétrique définie positive dense dont les coefficients sont aléatoires"""
    
    A = np.zeros([n,n])
    
    for i in range(0,n):
        for j in range(0,n):
            T[i][j] = random.randint(0,p)
    
    return np.dot(A.T,A)

def Cholesky(A):
    """ Prend un paramètre :
    A : matrice (A supposée symétrique définie positive
    Retourne la matrice T triangulaire inférieure telle que A = T * transpose(T)"""
    n = len(A)
    T = np.zeros([n,n])

    for i in range(0,n):
        
        for j in range(i,n):
            
            somme = 0
            for k in range(0,i):
                somme = somme + T[i][k]*T[j][k]

            if( i == j ):
                T[i][i] = np.sqrt( A[i][i] - somme )
            else:
                T[j][i] = ( A[i][j] - somme ) / T[i][i]
  
    return T

def CholeskyIncomplete(A):
    """ Prend un paramètre :
    A : matrice (suppose A symétrique définie positive)
    Retourne la matrice T triangulaire inférieure telle que A = T * transpose(T)"""
    n = len(A)
    T = np.zeros([n,n])

    for i in range(0,n):
        
        for j in range(i,n):
            if( A[i][j] != 0 ):
                somme = 0
                for k in range(0,i):
                    somme = somme + T[i][k]*T[j][k]
                      
                if( i == j ):
                    T[i][i] = np.sqrt( A[i][i] - somme )
                else:
                    T[j][i] = ( A[i][j] - somme ) / T[i][i]
  
    return T


def MatSDP (n, p, rang, epsilon):
    """ Prend deux paramètres :
    n : taille de la matrice carrée (n supposé > 0)
    p : nombre de termes extra diagonaux non-nuls (p supposé pair et <= n*(n-1)
    Retourne la matrice carrée symétrique définie positive en fonction de n et p """
    A = np.zeros([n,n])

    while( p/2 != 0 ):
        x = random.randint(0,n-2)
        y = random.randint(x+1,n-1)

        while( A[x][y] != 0 ):
            x = random.randint(0,n-2)
            y = random.randint(x+1,n-1)

        tmp = np.random.uniform(-rang,rang)
        A[x][y] = tmp
        A[y][x] = tmp
        p = p - 2

    for i in range(0,n):
        somme = 0
        for j in range(0,n):
            somme = somme + abs(A[i][j])
        if( somme != 0 ):
            A[i][i] = somme + np.random.uniform(epsilon, rang)

    return A

# FONCTIONS GRADIENT _______________________________________________________________________
def conjgrad (A,b,x):
    """Prend trois paramètres :
    A : une matrice symétrique (de taille n*n) définie positive
    b : un vecteur (de taille n)
    x : un vecteur pouvant être préconditionné ou initialisé à 0
    Retourne le vecteur x proche de la solution de Ax=b"""
    r = b - np.dot( A, x )
    p = r
    rsold = np.dot( r.T, r )

    for i in range(1,(10^6)):
        Ap = np.dot( A, p)
        alpha = rsold / ( np.dot(p.T,Ap) )
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = np.dot( r.T, r )
        if( np.sqrt( rsnew ) - (1e-10) < 0 ):
            break
        p = r + ( rsnew/rsold ) * p
        rsold = rsnew
    return x

def conjgradPreCond (A,b,x,Minv):
    """ Prend trois paramètres :
    A : une matrice symétrique (de taille n*n) définie positive
    b : un vecteur (de taille n)
    x : un vecteur pouvant être préconditionné ou initialisé à 0
    Minv : une matrice préconditionnant la matrice A (calculée avec Cholesky, LU ...)
    Retourne le vecteur x proche de la solution de Ax=b """
    r = b - np.dot( A, x ) #r => r0    
    z = np.dot( Minv, r ) #z => z0
    p = z
    
    rsold = np.dot( r.T, r ) #rsold : produit scalaire rk
    zold = z #zold => zk
             #znew => z(k+1)
    rold = r #rold => rk
             #rnew => r(k+1)
             
    for i in range(1,(10^6)):
        Ap = np.dot( A, p )
        alpha = np.dot( rold.T, zold )  / ( np.dot(p.T,Ap) )
        x = x + alpha * p
        rnew = rold - alpha * Ap
        rsnew = np.dot( rnew.T, rnew )
        if( np.sqrt( rsnew ) - (1e-10) < 0 ):
            break
        znew = np.dot( Minv, rnew )
        p = znew + (np.dot( znew.T, rnew ) / np.dot( zold.T, rold )) * p
        rsold = rsnew
        rold = rnew
        zold = znew
    return x

# FONCTIONS CHALEUR ________________________________________________________________________
def McChaleur (n):
    """ Prend un paramètre :
    n : taille de la matrice carrée => n*n * n*n
    Retourne la matrice des dérivées partielles de l'équation de la chaleur"""
    N = n * n
    Mc = np.zeros([ N, N ])

    for i in range(0,n): #indice de ligne du carré courant
        for j in range(0,n): #indice de colonne du carré courant
            
            if( i == j ): #dans les carrés n*n de la diagonale de la matrice N*N
                for k in range(0,n): #indice de ligne dans le carré courant
                    for l in range(0,n): #indice de colonne dans le carré courant
                        if( k == l ): #dans la diagonale de ces carrés
                            Mc[ i * n + k ][ j * n + l ] = -4.0 / (N + 1)
                        if( abs(k - l) == 1 ): #juste autour de la diagonale
                            Mc[ i * n + k ][ j * n + l ] = 1.0 / (N + 1) 
                            
            if( abs(i - j) == 1 ): #dans les carrés n*n autour des carrés n*n de la diagonale de la matrice N*N
                for k in range(0,n): #indice de ligne dans le carré courant
                    for l in range(0,n): #indice de colonne dans le carré courant
                        if( k == l ): #dans la diagonale de ces carrés
                            Mc[ i * n + k ][ j * n + l ] = 1.0 / (N + 1) 

    return Mc

def VecToMat (x):
    """ Prend un paramètre :
    x : vecteur de taille n*n
    Retourne la matrice n*n construite à partir de ce vecteur"""
    N = len(x)
    n = int(np.sqrt(N))
    T = np.zeros([ n, n ])

    for i in range(0,n):
        for j in range(0,n):
            T[i][j] = x[ (i * n) + j ]
            
    return T

def MatToVec (F):
    """ Prend un paramètre :
    F : matrice de taille n*n
    Retourne le vecteur de taille n*n construit à partir de la matrice"""
    n = len(F)
    N = n * n
    x = np.zeros([ N, 1 ])

    for i in range(0,n):
        for j in range(0,n):
            x[ (i * n) + j ] = F[i][j]

    return x

def Radiateur (n):
    """ Prend un paramètre :
    n : taille de la matrice => n*n
    Retourne la matrice comportant un radiateur (un carré de nombres !=0) au centre"""
    Rad = np.zeros([ n, n ])

    for i in range( n/2 - n/4, n/2 + n/4 ): #taille du radiateur = n/2 * n/2, placé au centre
        for j in range( n/2 - n/4, n/2 + n/4 ):
            Rad[i][j] = -20 #Chaleur du radiateur

    return Rad

def MurChaud (n):
    """ Prend un paramètre :
    n : taille de la matrice => n*n
    Retourne la matrice comportant un mur chaud (une ligne de nombres !=0) au Nord"""
    Mur = np.zeros([ n, n ])

    i = n-1 #Mur Nord
    for j in range(0,n):
        Mur[i][j] = -20 #Chaleur du mur

    return Mur


# TESTS CHOLESKY ___________________________________________________________________________
print "\nTESTS CHOLESKY :\n"

print "\nTest Cholesky sur la matrice A :"
A = [[6, 3, 4, 8], [3, 6, 5, 1], [4, 5, 10, 7], [8, 1, 7, 25]]
print np.asmatrix(A)

print "\nCholesky(A) :"
T = Cholesky(A)
print T

print "\nPréconditioneur :"
Minv = np.dot( np.linalg.inv(T).T, np.linalg.inv(T) )
print Minv

print "\nTests de préconditionneur :"

print "\nMatrice inv(A) :"
print np.linalg.inv(A)

print "\nMatrice inv(M) :"
print Minv

print "\nConditionnement de A :"
print np.linalg.cond(A)

print "\nConditionnement de inv(M) * A :"
print np.linalg.cond( np.dot( Minv, A ) )
print "On a ici inv(M)=inv(A) d'où cond( inv(M) * A ) = 1."

print "\nLe vrai Cholesky(A) :"
T = np.linalg.cholesky(A)
print T

print "\nTest matrice symétrique définie positive aléatoire (dense, p=n*(n-1)) :"
A = MatSDP(5,20,20.0,0.001)
print A

print "\nTest Cholesky Incomplet sur la matrice A :"

print "\nCholesky(A) :"
T = CholeskyIncomplete(A)
print T

print "\nLe vrai Cholesky(A) :"
print np.linalg.cholesky(A)

print "\nTest de résolution avec Cholesky :"
print "\nMatrice A :"
A = MatSDP(5,20,20.0,0.001)
print A

print "\nVecteur b :"
b = [ [1], [2], [3], [4], [5] ]
print np.asmatrix(b)

print "\nMatrice T = Cholesky(A) :"
T = Cholesky(A)
print T

print "\nRésolution de T*y=b :"
y = np.linalg.solve(T,b)
print y

print "\nRésolution de transpose(T)*x=y :"
x = np.linalg.solve(T.T,y)
print x

print "\nRésultat bon ?"
print np.allclose(np.dot(A, x), b)

print "\nLa vraie résolution Ax=b :"
x = np.linalg.solve(A,b)
print x

print "\nTest de résolution avec Cholesky Incomplet :"
print "\nMatrice A creuse :"
A = MatSDP(5,8,20.0,0.001)
print A

print "\nVecteur b :"
b = [ [1], [2], [3], [4], [5] ]
print np.asmatrix(b)

""" print "\nMatrice T = CholeskyIncomplete(A) :"
T = CholeskyIncomplete(A)
print T

print "\nRésolution de T*y=b :"
y = np.linalg.solve(T,b)
print y
                                               #La résolution d'un tel système par np.linalg.solve peut engendrer des erreurs en cas de matrice particulière (souvent) => abandon
print "\nRésolution de transpose(T)*x=y :"
x = np.linalg.solve(T.T,y)
print x

print "\nRésultat bon ?"
print np.allclose(np.dot(A, x), b)

print "\nLa vraie résolution Ax=b :"
x = np.linalg.solve(A,b)
print x """

# TESTS CONJGRAD ___________________________________________________________________________
print "\nTESTS CONJGRAD :\n"

print "\nMatrice A :"
A = MatSDP(5,20,20.0,0.001)
print A

print "\nVecteur b :"
b = [ [1], [2], [3], [4], [5] ]
print np.asmatrix(b)

print "\nVecteur x (initialisé à 0) :"
x = [ [0], [0], [0], [0], [0] ]
print np.asmatrix(x)

print "\nTest conjgrad(A,b,x) :"
x = conjgrad(A,b,x)
print x

print "\nRésultat bon ?"
print np.allclose(np.dot(A, x), b)

print "\nLa vraie résolution Ax=b :"
x = np.linalg.solve(A,b)
print x

print "\nTESTS CONJGRAD AVEC PRECONDITIONNEUR :\n"

print "\nMatrice A :"
A = MatSDP(5,20,20.0,0.001)
print A

print "\nPréconditionneur :"
T = Cholesky(A)
Minv = np.dot( np.linalg.inv(T).T, np.linalg.inv(T) )
print Minv


print "\nVecteur b :"
b = [ [1], [2], [3], [4], [5] ]
print np.asmatrix(b)

print "\nVecteur x (initialisé à 0) :"
x = [ [0], [0], [0], [0], [0] ]
print np.asmatrix(x)

print "\nTest conjgradPreCond(A,b,x,Minv) :"
x = conjgradPreCond(A,b,x,Minv)
print x

print "\nRésultat bon ?"
print np.allclose(np.dot(A, x), b)

print "\nLa vraie résolution Ax=b :"
x = np.linalg.solve(A,b)
print x

# TESTS CHALEUR ____________________________________________________________________________
print "\nTESTS CHALEUR :\n"

print "\nTest McChaleur(n) :"
n = 100 #Erreur mémoire au-dessus d'une précision de 150 pour l'image. Et fait freeze la machine au-dessus de 100
A = McChaleur(n)
print A

#print "\nTest résolution chaleur radiateur :"
F = Radiateur(n)
Fx = MatToVec(F)
x = np.linalg.solve(A,Fx)
T = VecToMat(x)

plt.title("Radiateur")
plt.imshow(T.T, cmap=cm.jet, interpolation='nearest', origin='lower')
plt.axis('off')
plt.savefig("Radiateur")
plt.show()

#print "\nTest résolution chaleur mur :"
F = MurChaud(n)
Fx = MatToVec(F)
x = np.linalg.solve(A,Fx)
T = VecToMat(x)

plt.title("Mur chaud")
plt.imshow(T, cmap=cm.jet, interpolation='nearest', origin='lower')
plt.axis('off')
plt.savefig("Mur_Chaud")
plt.show()


#ERREUR RELATIVE ___________________________________________________________________________

def erreurRelatCholesky (n, p, rang, epsilon, it):
    erreurs = []
    conds = []

    b = np.zeros([n,1])
    for j in range(0,n):
        b[j][0] = random.randint(-rang,rang)
    
    for i in range (0,it):
        A = MatSDP (n,p,rang,epsilon)
    
        T = CholeskyIncomplete (A)
        y = np.linalg.solve(T,b)
        xexp = np.linalg.solve(T.T,y)
        xtheo = np.linalg.solve(A,b)

        erreurs.append(np.linalg.norm(np.dot(A,xexp) - np.linalg.norm(b)) / ( np.linalg.norm(np.dot(A,xexp)) + (np.linalg.norm(b)) ))
        conds.append(np.linalg.cond(A))
    plt.title("Erreur relative Cholesky")
    plt.plot(conds,erreurs,'ro')
    plt.savefig("erreur_relative_Cholesky")
    plt.show()
    





def erreurRelatConjGrad (n, p, rang, epsilon, it):
    erreurs = []
    conds = []


    x = np.zeros([n,1])
    b = np.zeros([n,1])
    for j in range(0,n):
        b[j][0] = random.randint(-rang,rang)
    
    for i in range (0,it):
        A = MatSDP (n,p,rang,epsilon)
       

        xexp = conjgrad(A,b,x)
        xtheo = np.linalg.solve(A,b)

        erreurs.append(np.linalg.norm(np.dot(A,xexp) - np.linalg.norm(b)) / ( np.linalg.norm(np.dot(A,xexp)) + (np.linalg.norm(b)) ))
        conds.append(np.linalg.cond(A))
    plt.title("Erreur relative ConjGrad")
    plt.plot(conds,erreurs,'ro')
    plt.savefig("erreur_relative_ConjGrad")
    plt.show()



#test Cholesky
erreurRelatCholesky(5,20,20.0,0.001,1000)

#test ConjGrad
erreurRelatConjGrad(5,20,20.0,0.001,1000)


