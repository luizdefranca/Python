# -*- coding: utf-8 -*-

import numpy as np
import random as rd
import matplotlib.image as img
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# 1. Housholder
"""
Fonction qui prend deux vecteurs X,Y de même norme et renvoie un vecteur U de norme euclidienne permettant d'obtenir la matrice de Householder
"""
def vector_U(X,Y):
    U = np.subtract( X, Y)
    if (np.linalg.norm(U)==0):
        n=len(U)
        return np.identity(n)
    else:
        U = U/(np.linalg.norm(U))
        return U

"""
Fonction qui prend les vecteurs X et Y et renvoie la matrice de Householder correspondante
"""
def Householder(X,Y):
    U = vector_U(X,Y)
    n = len(U)
    I = np.eye(n)
    #print np.outer(U,U.T)
    return I -( 2 * U * (U.T))
"""
Fonction calculant le produit d'une matrice de Housholder et d'un vecteur
"""
def opti_dot_vect (H, X):
    """Complexité n*(n-1)/2 au lieu de n*n"""
    n = len(H)
    V = np.zeros([n,1])
    for i in range(0,n):
        for j in range(i,n):
            V[i] = V[i] + H[i][j] * X[j]
            if( j>i ):
                V[j] = V[j] + H[j][i] * X[i]
    return V

"""
Fonction calculant le produit d'une matrice de Housholder et d'une autre matrice
"""
def opti_dot_mat (H,A):
    """Complexité de n*n*(n-1)/2 au lieu de n*n*n"""
    n = len(H)
    for i in range(0,n):
        Vold = opti_dot_vect( H, A[:,i] )
        if( i == 0 ):
            Vnew = Vold
        else:
            Vnew = np.vstack((np.array(Vnew),np.array(Vold)))
    return np.resize(Vnew,(n,n)).T  

# 2. Mise sous forme bidiagonale
"""
Algorithme de bidiagonalisation
Ne fonctionne pas. Une erreur part de la première colonne, et se propage, de plus en plus colonne par colonne. Mais le résultat reste quand même bidiagonal
"""
def bidiagonalization (A, Qleft, Qright):
    n = len(A)
    Qleft = np.eye(n)
    Qright = np.eye(n)
    BD = A
    for i in range(0,n):
        V = np.zeros([n-i,1])
        V[0] = 1 #np.linalg.norm(BD[i:n,i])
        Q1 = Householder( BD[i:n,i], V )
        Qleft[i:n,i:n] = np.dot( Qleft[i:n,i:n], Q1 )
        BD[i:n,i:n] = np.dot( Q1, BD[i:n,i:n] )
        for j in range(i+2,n):
            BD[i,j] = 0
        if( i <= (n-2) ):           
            V = np.zeros([1,n-i-1])
            V[0][0] = 1 #np.linalg.norm(BD[i,(i+1):n])
            Q2 = Householder( BD[i,(i+1):n], V )
            print Q2
            print Qright[(i+1):n,(i+1):n]
            Qright[(i+1):n,(i+1):n] = np.dot( Q2, Qright[(i+1):n,(i+1):n] )
            BD[(i+1):n,(i+1):n] = np.dot( BD[(i+1):n,(i+1):n], Q2 )        
        #print np.dot(Qleft, np.dot( BD, Qright) )
    return (Qleft,BD,Qright)

#3. Transformation QR
"""
Fonction incluant une portion de matrice carrée dans une autre.
Renvoie la matrice identité, modifiée par les coeficients de l'autre matrie à inclure.
"""
def include (A,n):
    m = len(A)
    I = np.eye(n)
    for i in range(0,m):
        for j in range(0,m):
            I[i+n-m,j+n-m] = A[i,j]
    return I

"""
Transformation QR, méthode itérative
"""
def QR (A):
    n = len(A)
    Rold = A
    for i in range(0,n):
        u = Rold[i:n,i]
        v = np.zeros([n-i,1])
        v[0] = np.linalg.norm(u)
        Qold = Householder( u, v )
        Qnew = include( Qold, n )
        if( i==0 ):
            Q = Qnew
        else:
            Q = np.dot(Q,Qnew)
        Rnew = np.dot( Qnew, Rold )
        Rold = Rnew
    return (Q,Rnew)        

def transUSV (B, nmax):
    """Application de la transformation USV à la matrice Bidiagonale B"""
    n = len(B)
    U = np.identity(n)
    V = np.identity(n)
    S = B
    for i in range(0, nmax+1):
        (Q1,R1) = np.linalg.qr(S.T)
        (Q2,R2) = np.linalg.qr(R1.T)
        S = R2
        U = np.dot(U,Q2)
        V = np.dot(Q1.T,V)
    return U, S, V

# Marche pas à cause de QR
# def transUSV (B, nmax):
#     """Application de la transformation USV à la matrice Bidiagonale B"""
#     n = len(B)
#     U = np.identity(n)
#     V = np.identity(n)
#     S = B
#     for i in range(0, nmax+1):
#         (Q1,R1) = QR(S.T)
#         (Q2,R2) = QR(R1.T)
#         S = R2
#         U = np.dot(U,Q2)
#         V = np.dot(Q1.T,V)
#     return U, S, V

def genMatrixBD (n,p,pos):
    """Génère des matrices bidiagonales"""
    A = np.zeros([n,n])
    for i in range (0, n):
        for j in range (0,n):
            if (i == j):
                A[i][j] = rd.randint(0,p)
            if (pos == 1):
                if (i == j-1):
                    A[i][j] = rd.randint(0,p)
            else:
                if (i == j+1):
                    A[i][j] = rd.randint(0,p)
    return A

def erreurRelMat(A, B):
    """Donne la moyenne de l'erreur relatives termes à termes des
    coefs des mats"""
    n = len(A)
    s = 0
    for i in range (0,n):
        for j in range (0,n):
            s = s + abs(A[i][j] - B[i][j])
    m = s /(n*n)
    return m

def erreurD(S) :
    """Donne l'erreur de S et de S matrice diagonale"""
    n = len(S)
    res = 0
    for i in range (0,n):
        for j in range (0,n):
            if (i!=j):
                res = res + S[i][j]
    return res/(n*n)

def erreurBD(S) :
    """Donne l'erreur de S et de S matrice Bidiagonale"""
    n = len(S)
    res = O
    for i in range (0,n):
        for j in range (0,n):
            if ((i!=j) and (i != j-1)):
                res = res + S[i][j]
    return res/(n*n)

def convergenceD(A, iter):
    """Prends une matrice bidiagonale et renvoie le graphe de la
    convergence de S vers une matrice diagonale en fonction
    du nombre d'itération"""
    y = []
    for i in range (0, iter+1):
        (B,C,D) = transUSV (A,i)
        y.append(erreurD(C))
    x = np.arange(0, iter+1, 1)
    err = y[0]
    #print y
    plt.plot(x,y)
    # plt.axis([0,10,err-0.000000000000001,err+0.000000000000001])
    plt.show()

def convergenceBD(A, iter):
    """Prends une matrice bidiagonale et renvoie le graphe de la
    convergence de S vers une matrice bidiagonale en fonction
    du nombre d'itération"""
    y = []
    for i in range (0, iter+1):
        (B,C,D) = transUSV (A,i)
        y.append(erreurBD(C))
    x = np.arange(0, iter+1, 1)
    err = y[0]
    #print y
    plt.plot(x,y)
    # plt.axis([0,10,err-0.000000000000001,err+0.000000000000001])
    plt.show()

def invariantUSV(A, iter):
    """Prends une matrice bidiagonale et renvoie un graphe en fonction
    du nombre d'itération"""
    y = []
    for i in range (0, iter+1):
        (B,C,D) = transUSV (A,i)
        E = np.dot(B, np.dot(C,D))
        y.append(erreurRelMat(A,E))
    x = np.arange(0, iter+1, 1)
    err = y[0]
    #print y
    plt.plot(x,y)
    # plt.axis([0,10,err-0.000000000000001,err+0.000000000000001])
    plt.show()

# 4 . Compression d'images

"""Décomposition d’une matrice A en matrices U,S,V avec U et V orthogonales carrées et S diagonale réelle """
def compression_rang_k(A, k):
    SVD = np.linalg.svd(A,full_matrices=False) 
    U = SVD[0]
    S = SVD[1]
    V = SVD[2]
    n = len(S)
    for i in range (k, n):
        S[i] = 0.0
    S1 = np.diag(S)
    A = np.dot(U, np.dot(S1,V))
    return A

"""Compression d'une image"""
def compression_img(M,k):
    M[:,:,0] = compression_rang_k(M[:,:,0],k) # On modifie la composante R
    M[:,:,1] = compression_rang_k(M[:,:,1],k) # On modifie la composante G
    M[:,:,2] = compression_rang_k(M[:,:,2],k) # On modifie la composante B
    imgplot = plt.imshow(M)
    plt.show()

"""Compression sans affichage"""
def compression_img_ws(M,k):
    M[:,:,0] = compression_rang_k(M[:,:,0],k)
    M[:,:,1] = compression_rang_k(M[:,:,1],k)
    M[:,:,2] = compression_rang_k(M[:,:,2],k)
    return M

"""Efficacité en fonction de k"""
def efficiency(M):
    Inv = M
    n = len(M)
    m = len(M[0])
    y = []
    for k in range (0,(n/2)):
        Mc = compression_img_ws(M,k)
        A = Inv-Mc
        res = np.linalg.norm(A)
        y.append(res)
    x = np.arange(0,n/2,1)
    plt.plot(x,y)
    plt.show()
    
# Tests 1. Housholder   
X = np.array([[3],[4],[0]])
Y = np.array([[0],[0],[5]])
print X
print Y
print "Vecteur U(X,Y)"
print vector_U( X, Y )
H = Householder(X,Y)
print "Matrice de Housholder envoyant X sur Y :"
print H
print "Vérification de la multiplication Housholder * Vecteur optimisée"
print np.dot(H,X)
print opti_dot_vect(H,X)

print "Définiton d'une matrice A"
A = np.matrix([[-4, 20, 35], [-4, -30, -15], [-8, 40, -80]])
print A
print "Vérification de la multiplication Housholder * Matrice"
print np.dot(H,A)
print opti_dot_mat(H,A)

# Tests 2. Mise sous forme bidiagonale
print "Définiton d'une matrice A"
# A = np.matrix([[-4, 20, 35, 5], [-4, -30, -15, 55], [-8, 40, -80, -65], [23, -15, 30, 15]])
A = np.matrix([[12,-51,4], [6, 167, -68], [-4, 24, -41]])
print A
n = len(A)
Qleft = np.eye(n)
Qright = np.eye(n)

#Qleft = bidiagonalization(A,Qleft,Qright)[0] #Qleft
#BD = bidiagonalization(A,Qleft,Qright)[1] #BD
#Qright = bidiagonalization(A,Qleft,Qright)[2] #Qright
#print "Tests de bidiagonalisation : affiche Qr, BD et Ql tels que A = Ql * BD * Qr"
#print Qleft
#print BD
#print Qright
#print "On vérifie qu'on retombe sur A :"
#print np.dot(Qleft, np.dot( BD, Qright) )

# Tests 3. Transformation QR
print "Tests de notre fonction QR"
Q = QR(A)[0]
R = QR(A)[1]
print Q
print R
print "Résultat avec l'algorithme QR de numpy (test si le même que précédemment)"
Q = np.linalg.qr(A)[0]
R = np.linalg.qr(A)[1]
print Q
print R

print "Génération d'une matricé bidiagonale aléatoire A"
A = genMatrixBD (5,54,1)
(B,C,D) = transUSV (A, 100)
print A
print "Test de la décomposition SVD, on regarde si on retombe sur A avec A = U*S*V"
E = np.dot(B,np.dot(C,D))
print E

print "Erreur relative entre la matrice A et sa décomposition en USV"
print (erreurRelMat(A,E))
print (erreurRelMat(A,A))
invariantUSV (A, 100)
convergenceD (A, 100)
M = img.imread("essai.jpeg")
compression_img(M, 50)

M = img.imread("essai.jpeg")
#efficiency(M)
