import numpy as np;
import matplotlib.pyplot as plt


def step_Euler(y0, t0, h, f):
    return y0 + h*f(y0,t0)


def step_Rung_Kutta_4(y0, t0, h, f):
    pn1 = f(y0,t0)
    yn1 = y0 + (h/2) * pn1

    pn2 = f(yn1, t0 + (h/2))
    yn2 = y0 + h * pn2 / 2

    pn3 = f(yn2, t0 + (h/2))
    yn3 = y0 + h * pn3

    pn4 = f(yn3, t0 + h)
    return y0 + (h/6) * (pn1 + 2*pn2 + 2*pn3 + pn4)
    


def meth_n_step(y0,t0,N,h,f,meth):
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


#test 

def test1():
    """un test sur l'exponentielle"""
    d = lambda y,t: y
    
    print step_Euler(1, 0.25, 0.25, d);
    
    print meth_n_step(1, 0, 10, 0.25, d, step_Euler);
    
    Y = meth_n_step(1, 0, 10, 0.25, d, step_Rung_Kutta_4);
    Z = meth_n_step(1, 0, 10, 0.25, d, step_Euler);
    print step_Rung_Kutta_4(1, 0, 0.25, d)
    
    x = np.linspace(0,2,10)
    
    print np.size(x)
    
    p1 = plt.plot(x, np.exp(x),label='exact')
    p2 = plt.plot(x, Y,label='runge')
    p3 = plt.plot(x, Z,label='euler')
    plt.legend();
    
    plt.legend();
    plt.show();


def test2():
    """un test sur exp(-x**2)"""
    d = lambda y,t: -2*t*y
    
    print step_Euler(1, 0.25, 0.25, d);
    
    print meth_n_step(1, 0, 10, 0.25, d, step_Euler);
    
    Y = meth_n_step(1, 0, 10, 0.25, d, step_Rung_Kutta_4);
    Z = meth_n_step(1, 0, 10, 0.25, d, step_Euler);
    print step_Rung_Kutta_4(1, 0, 0.25, d)
    
    x = np.linspace(0,2,10)
    
    print np.size(x)
    
    p1 = plt.plot(x, np.exp(-x**2),label='exact')
    p2 = plt.plot(x, Y,label='runge')
    p3 = plt.plot(x, Z,label='euler')
    plt.legend();
    
    plt.legend();
    plt.show();

test1();
test2();





#__________________________Modele de population_________________________

#Malthus


#Verrhulst

def test_verrhulst():
    d = lambda y,t: 3*y*(1-y/6)
    
    print step_Euler(2, 0.25, 0.25, d);
    
    print meth_n_step(2, 0, 10, 0.25, d, step_Euler);
    
    Y = meth_n_step(1, 0, 10, 0.25, d, step_Rung_Kutta_4);
    Z = meth_n_step(1, 0, 10, 0.25, d, step_Euler);
    print step_Rung_Kutta_4(1, 0, 0.25, d)
    
    x = np.linspace(0,2,10)
    
    print np.size(x)
    
    #p1 = plt.plot(x, 6/(1+ (1-1/3)*np.exp(3*x)),label='exact')
    p2 = plt.plot(x, Y,label='runge')
    p3 = plt.plot(x, Z,label='euler')
    plt.legend();
    
    plt.legend();
    plt.show();


test_verrhulst();


#Lotka Volterra


f = lambda P,N,t: N*(5 - 1*P); 
g = lambda P,N,t: P*(1*N - 1);


def step_Euler_system(N0, P0, t0, h, Fn, Fp):
    return N0 + h*Fn(P0,N0,t0), P0 + h*Fp(P0,N0,t0);

print step_Euler_system(100, 20, 0, 0.25, f, g)

def step_Rung_Kutta_4_system(N0, P0, t0, h, Fn, Fp):
    pn1N = Fn(P0, N0, t0)
    yn1N = N0 + (h/2) * pn1N
    pn1P = Fp(P0, N0, t0)
    yn1P = P0 + (h/2) * pn1P

    pn2N = Fn(yn1P,yn1N, t0 + (h/2))
    yn2N = N0 + h * pn2N / 2
    pn2P = Fp(yn1P,yn1N, t0 + (h/2))
    yn2P = P0 + h * pn2P / 2

    pn3N = Fn(yn2P, yn2N, t0 + (h/2))
    yn3N = N0 + h * pn3N
    pn3P = Fp(yn2P, yn2N, t0 + (h/2))
    yn3P = P0 + h * pn3P

    pn4N = Fn(yn3P, yn3N, t0 + h)
    pn4P = Fp(yn3P, yn3N, t0 + h)
  
    N =  N0 + (h/6) * (pn1N + 2*pn2N + 2*pn3N + pn4N)
    P =  P0 + (h/6) * (pn1P + 2*pn2P + 2*pn3P + pn4P)
    
    return N,P
    

print step_Rung_Kutta_4_system(100, 20, 0, 0.25, f, g)


def meth_n_step_system(N0, P0, t0, h, Fn, Fp, n, meth):
    N = []
    P = []   
    N.append(N0)
    P.append(P0)
    t = t0
    for i in range (0, n):
        t = t + h
        tab = meth(N[i-1], P[i-1], t, h, Fn, Fp)

        Nnew = tab[0]
        N.append(Nnew)

        Pnew = tab[1]
        P.append(Pnew)        
    return N,P   






#periodicity

def trapezes(t0,h,N):
    S = 0;
    n = np.size(N);
    for i in range(0,n-1):
        x1 = t0 + h*i;
        x2 = t0 + h*(i+1)
        S += ((N[i]+N[i+1])/2.0)*(x2-x1);
    return S;
         

def Average(N):
    n = 0;
    n = np.size(n);
    S = 0;
    for i in range (0,n):
        S += N[i];
    return S/n;


def period(t0,h,N):
    return trapezes(t0,h,N)*h/Average(N);


#ploting


Tab = meth_n_step_system(5, 3, 0, 0.05, f, g, 500, step_Euler_system)
Rtab = meth_n_step_system(5, 3, 0, 0.05, f, g, 500, step_Rung_Kutta_4_system)

print "periode", period(0,0.05,Rtab[1]);

x = np.linspace(0,501,501)

plt.xlabel('h*t')
p1 = plt.plot(x,Rtab[0],label='N(t) avec RK4')
p2 = plt.plot(x,Rtab[1],label='P(t) avec RK4')

p3 = plt.plot(x,Tab[0],label='N(t) avec Euler')
p4 = plt.plot(x,Tab[1],label='P(t) avec Euler')
plt.legend();
plt.show();


plt.xlabel('N(t)')
plt.ylabel('P(t)')
p4 = plt.plot(Tab[0],Tab[1],label='Euler')
p1 = plt.plot(Rtab[0],Rtab[1],label='Runge Kutta 4')
plt.legend();
plt.show();


#multiple plot

def multi_plot (f,meth,n):
    """trace n courbes pour n y0 differents"""
    plt.ylabel('f(t)')
    for i in range (0,n-1):
         Rtab = meth_n_step(i, 0, 500, 0.05, f, meth)
         plt.plot(Rtab)
    plt.legend();
    plt.show();



v = lambda y,t: y

multi_plot(v,step_Rung_Kutta_4,10)



def multi_plot_system (t0, f, g, meth,n):
    """trace n courbes pour n y0 differents"""
    plt.xlabel('N(t)')
    plt.ylabel('P(t)')
    h=0;
    for i in range (0,n-1):
        h += i*0.01 
        Rtab = meth_n_step_system(h, h+1, t0, 0.05, f, g, 500, meth)
        plt.plot(Rtab[0],Rtab[1])
    plt.legend();
    plt.show();

multi_plot_system(3,f,g,step_Rung_Kutta_4_system,10)




