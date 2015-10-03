import numpy as np;



def Newton_Raphson(f, J, U0, N, epsilon):
    U = U0
    i = 0
    while( (i < N) and (np.linalg.norm( f(*U) ) >= epsilon) ):
        V = np.linalg.lstsq(J, -f(*U))[0]
        U = U + V
        i = i + 1
    return U



def Bairstow(P):
    k = 0; # compteur d'iteration;
    A = P.coeffs;
    B = 1; C = -8;
    #B = A[1]/A[0];
    #C = A[2]/A[0];
    epsilon = 10**-10;
    V = np.zeros(2);
    while ((abs(P(V[1])) > epsilon) & (abs(P(V[0])) > epsilon)):
        U = np.array([B, C]);
        T = np.poly1d([1, B, C]);
        
        Div = np.polydiv(P,T) #1ere div
        Reste = Div[1];
        Q = Div[0];       
        R = Reste[0];
        S = Reste[1];
        
        Div = np.polydiv(Q,T); #2nde div
        Reste = Div[1];
        G = Div[0];
        Rc = -Reste[0];
        Sc = -Reste[1];
        
        Rb = -B*Rc + Sc;
        Sb = -C*Rc;
        
        dv = 1.0/(Sb*Rc -Sc*Rb);
        delb = (R*Sc - S*Rc)*dv;
        delc = (-R*Sb + S*Rb)*dv;
        diff = np.array([delb, delc])
        
        
        B = B + diff[0];
        C = C + diff[1];        
        T = np.poly1d([1, B, C]);
        V = T.roots;
        k = k+ 1;
    return V,k;



def Bairstow2(P): # Ici la resolution se fait avec Newton_Raphson
    k = 0; 
    B = 100;
    C = -110;
    epsilon = 10**-10;
    V = np.zeros(2);
    while (abs(P(V[1])) > epsilon and abs(P(V[0])) > epsilon):
        U = np.array([B, C]);
        T = np.poly1d([1, B, C]);
        
        Div = np.polydiv(P,T) #1ere div
        Reste = Div[1];
        Q = Div[0];       
        R = Reste[0];
        S = Reste[1];
        
        Div = np.polydiv(Q,T); #2nde div
        Reste = Div[1];
        G = Div[0];
        Rc = -Reste[0];
        Sc = -Reste[1];
        
        Rb = -B*Rc + Sc;
        Sb = -C*Rc;
        
        f = lambda x,y: np.array([R,S]);
        J = np.array([[Rb, Rc], [Sb, Sc]]);
        diff = Newton_Raphson(f, J, U, 100, 10**-3);
        
        B = B + diff[0];
        C = C + diff[1];        
        T = np.poly1d([1, B, C]);
        V = T.roots;
        k = k+ 1;
    return V,k;


P = np.poly1d([6,11,-33,-33,11,6]);
#print "P(x) =", P 
print "les racines avec P.roots", P.roots
print "les racines avec Bairstow",Bairstow(P)[0];
print "\n";
print "evaluation de P(x+) \t",P(Bairstow(P)[0][0])
print "evaluation de P(x-) \t",P(Bairstow(P)[0][1]);
print "nombre d'iterations \t",Bairstow(P)[1];
print "\n";


#P2 = np.poly1d([3, -29, 57, 241, -700, 348]); 
#print "racine",Bairstow(P2);
#print P(Bairstow(P2)[0])
#print P(Bairstow(P2)[1]);




#_______________________________Comparaison avec une methode de newton raphson en dim 1______________________


def diff_P(P):
    n = P.order;
    R = np.zeros(n);
    A = P.coeffs;
    for i in range (0,n):
        R[i] = A[i]*(n-i);
    Q = np.poly1d(R);
    return Q;

print P;
print "\n";
print "test de derivee de P";
print diff_P(P);

def Newton_Solve(P):
    k = 0;
    x = 3;
    epsilon = 10**-100;
    while abs(P(x)) > epsilon:
        x = x - P(x)/diff_P(P)(x);
        k += 1; 
    return x,k;

print "test de resolution par N-R en dim1:";
print " ->racine trouvé \t",  Newton_Solve(P)[0];
print " ->nombre d'iterations \t",  Newton_Solve(P)[1];



