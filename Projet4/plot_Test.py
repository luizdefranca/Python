#__________________________________________Plot Test______________________________
import numpy as np;
import numpy as py;
import matplotlib as mp;
from pylab import *


def Bairstow(P,epsilon):
    k = 0; # compteur d'iteration;
    A = P.coeffs;
    #B = 1; C = -2;
    B = A[1]/A[0];
    C = A[2]/A[0];
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


def Bairstow2(P,epsilon,B,C):
    k = 0; # compteur d'iteration;
    A = P.coeffs;  
    C = A[2]/A[0];
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



def diff_P(P):
    n = P.order;
    R = np.zeros(n);
    A = P.coeffs;
    for i in range (0,n):
        R[i] = A[i]*(n-i);
    Q = np.poly1d(R);
    return Q;


def Newton_Solve(P,epsilon):
    k = 0;
    x = 3;
    while abs(P(x)) > epsilon:
        x = x - P(x)/diff_P(P)(x);
        k += 1; 
    return x,k;



def test():
    P = np.poly1d([6,11,-33,-33,11,6]);
    epsilon = np.array([1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,100]);
    iteration1 = np.zeros(17);
    iteration2 = np.zeros(17);
    iteration3 = np.zeros(17);
    for i in range(0,17):
        B = Bairstow(P,pow(10,-epsilon[i]));
        iteration1[i] = B[1]; 

        B2 = Bairstow2(P,pow(10,-epsilon[i]),10,-8);
        iteration3[i] = B2[1]; 

        N = Newton_Solve(P,10**-epsilon[i]);
        iteration2[i] = N[1];
    return iteration1, iteration2,iteration3;


epsilon = np.array([1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,100]);
y1 = test()[0];
y2 = test()[1];
y3 = test()[2];

#mp.plot(epsilon, y1, linewidth=1.0)
xlabel('-log(epsilon)')
ylabel('iterations')
legend('B = 1; C = -1; \n x0 = 1;')
title('comparison between Newton method and Bairstow method')

plot(epsilon, y3,'b-o', label='Bairstow without optimization', color = "red");
legend();
plot(epsilon, y1,'b-o', label='Bairstow with optimization');
legend();
plot(epsilon, y2,'b-o', label='Newton', color = 'green');
legend();
show();
#mp.savefig("figure");
