# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math

#Other functions
def apply_Jacobian (J,U):
    n = len(U)
    V = np.zeros([n,n])
    for i in range(0,n):
        for j in range(0,n):
            V[i][j] = J[i][j](*U)
    return V 

#Question 1
#Equation linking H, U and V
#f(U)+H(U)*V = 0 <=> H(U)*V = -f(U)

#Question 2
def Newton_Raphson(f, J, U0, N, epsilon):
    # X = []
    # Y = []
    U = U0
    i = 0
    while( (i < N) and (np.linalg.norm( f(*U) ) >= epsilon) ):
        # X.append(i)
        # Y.append(np.linalg.norm(f(*U)))
        V = np.linalg.lstsq(apply_Jacobian(J,U), - f(*U))[0]
        U = U + V
        i = i + 1
    # plt.title("Newton-Raphson without backtracking")
    # plt.plot(X,Y,'ro')
    # plt.xlabel("Iterations")
    # plt.ylabel("||f(U)||")
    # plt.savefig("Newton-Raphson_no_backtracking")
    # plt.show()
    return U


#Question 3
print "-----------TESTS NEWTON-RAPHSON-----------------"
print "f1:"
f1 = lambda x: 2.0*x
print f1
print "\nf (R->R):"
function_test = lambda x: np.array([ f1(x) ])
print function_test
print "\nJ :"
J_test = [[lambda x: 2.0]]
print J_test
print "\nU0:"
U0_test = np.array([ 4.0 ])
print U0_test
N = 100
epsilon = 10e-6

print "\nSolution calculated with Newton-Raphson algorithm "
print Newton_Raphson(function_test, J_test, U0_test, N, epsilon)

#Question 4
def Newton_Raphson_bt(f, J, U0, N, epsilon):
    # X = []
    # Y = []
    U = U0
    i = 0
    while( (i < N) and (np.linalg.norm( f(*U) ) >= epsilon) ):
        # X.append(i)
        # Y.append(np.linalg.norm(f(*U)))
        V = np.linalg.lstsq(apply_Jacobian(J,U), - f(*U))[0]
        while( np.linalg.norm(f(*(U+V))) >= np.linalg.norm(f(*U)) ):
            V = V / 2.0 #comes back to a lower value of V in the same direction (a half)
        U = U + V
        i = i + 1
    # plt.title("Newton-Raphson with backtracking")
    # plt.plot(X,Y,'ro')
    # plt.xlabel("Iterations")
    # plt.ylabel("||f(U)||")
    # plt.savefig("Newton-Raphson_backtracking")
    # plt.show()
    return U

# LAGRANGIAN
#Question 1
def elastic_force (x0,k):
    f = lambda x: np.array([ lambda x : -k*(x-x0) ])
    U0 = np.array([x0])
    J = [[lambda x: -k]]
    return np.array([f,J,U0])

def centrifugal_force (x0,y0,k):
    f1 = lambda x: k*(x-x0)
    f2 = lambda y: k*(y-y0)
    f = lambda x,y: np.array([f1(x),f2(y)])
    U0 = np.array([x0,y0])
    J = [[lambda x,y: k, lambda x,y: 0],[lambda x,y :0, lambda x,y :k]]
    return np.array([f,J,U0])

def gravitational_force (x0,y0,k):
    f1 = lambda x,y: -k*(x-x0)/((x-x0)**2.0 + (y-y0)**2)**(3.0/2.0)
    f2 = lambda x,y: -k*(y-y0)/((x-x0)**2.0 + (y-y0)**2)**(3.0/2.0)
    f = lambda x,y: np.array([f1(x,y),f2(x,y)])
    U0 = np.array([x0+1, y0+1]) #avoid a 0 division
    J = [[lambda x,y: 0,lambda x,y: 0],[lambda x,y: 0,lambda x,y: 0]]
    J[0][0] = lambda x,y: (-k/((x-x0)**2.0+(y-y0)**2.0)**(3.0/2.0)+(3.0/2.0)*k*(x-x0)*(2.0*x-2.0*x0)/((x-x0)**2.0+(y-y0)**2.0)**(5.0/2.0))
    J[1][1] = lambda x,y: (-k/((x-x0)**2.0+(y-y0)**2.0)**(3.0/2.0)+(3.0/2.0)*k*(y-y0)*(2.0*y-2.0*y0)/((x-x0)**2.0+(y-y0)**2.0)**(5.0/2.0))
    J[0][1] = lambda x,y: ((3.0/2.0)*k*(x-x0)*(2.0*y-2.0*y0)/((x-x0)**2.0+(y-y0)**2.0)**(5.0/2.0))
    J[1][0] = lambda x,y: ((3.0/2.0)*k*(y-y0)*(2.0*x-2.0*x0)/((x-x0)**2.0+(y-y0)**2.0)**(5.0/2.0))
    return np.array([f,J,U0])

#Question 2
print "-----------TESTS LAGRANGIAN-----------------"
x1 = 0.0
y1 = 0.0
m1 = 1.0

x2 = 1.0
y2 = 0.0
m2 = 0.01

R = x2 - x1
print "First gravitational force originating from (0,0) with a 1 coefficient:" 
fg1 = gravitational_force(x1,y1,m1)
print fg1
print "\nSecond gravitational force originating from (1,0) with a 0.01 coefficient:"
fg2 = gravitational_force(x2,y2,m2)
print fg2

barycentre_x = (m1*x1 + m2*x2)/(m1 + m2)
barycentre_y = (m1*y1 + m2*y2)/(m1 + m2)

kc = 1

print "\nCentrifugal force centered on the barycenter of the two masses, with coefficient 1:"
fc = centrifugal_force(barycentre_x,barycentre_y,kc)
print fc

print "\nA test point U:"
U = [1.5,0.0]
print U

print "\nThe Jacobian matrix of the three forces, applied to U:"
J_apply = apply_Jacobian(fg1[1],U) + apply_Jacobian(fg2[1],U) + apply_Jacobian(fc[1],U)
print J_apply

print "\nThe three force sum at point U:"
print fc[0](*U) + fg1[0](*U) + fg2[0](*U)

#TESTS Newton
print "-----------TESTS NEWTON-RAPHSON-----------------"
epsilon = 10e-12
max_iter = 1000000

print "What Newton-Raphson algorithm obtains for conditions above (three forces):"
print Newton_Raphson(fg1[0],fg1[1],fg1[2],max_iter,epsilon)

print "\nWhat the backtracking algorithm obtains:"
print Newton_Raphson_bt(fg1[0],fg1[1],fg1[2],max_iter,epsilon)

#Question 3
print "-----------LAGRANGIAN POINTS-----------------"
# def addFunctions(f, g, h):
#     return lambda x,y: f(x,y) + g(x,y) + h(x,y)

# U0_L = np.array([ 0.15, 0.15 ])
# J = fg1[1] + fg2[1] + fc[1]
# f = addFunctions( fg1[0], fg2[0], fc[0] )

# R123 = Newton_Raphson(f, J, U0_L, max_iter, epsilon)
#Overflow error after 53 or 54 loops

U0_L = np.array([ 0.15 ])

L1 = lambda r: np.array([m2/(r**2.0) + (m1*R/(m1+m2) - r)*(m1+m2)/(R**3.0) - m1/((R - r)**2.0)])
J_L1 = [ [lambda r: 2.0*m1/((r-R)**3.0) - 2.0*m2/(r**3.0) - (m1+m2)/(R**3.0)] ]
R1 = Newton_Raphson_bt(L1,J_L1,U0_L,max_iter,epsilon)

L2 = lambda r: np.array([(m1*R/(m1+m2) + r)*(m1+m2)/(R**3.0) - m2/(r**2.0) - m1/((R + r)**2.0)])
J_L2 = [ [lambda r: 2.0*m1/((r+R)**3.0) + 2.0*m2/(r**3.0) + (m1+m2)/(R**3.0)] ]
R2 = Newton_Raphson_bt(L2,J_L2,U0_L,max_iter,epsilon)

L3 = lambda r: np.array([(m2*R/(m1+m2) + R - r)*(m1+m2)/(R**3.0) - m1/(R - r)**2.0 - m1/(R - r)**2.0])
J_L3 = [ [lambda r: 4.0*m1/((r-R)**3.0)-(m1+m2)/(R**3.0)] ]
R3 = Newton_Raphson_bt(L3,J_L3,U0_L,max_iter,epsilon)

y123 = 0.0

x45 = (x2 - x1) / 2.0

y4 = math.sqrt( (x2 - x1)**2.0 - x45**2.0 )
y5 = -y4

X = np.array([R1,R2,R3,x45,x45,x1,x2])
Y = np.array([y123,y123,y123,y4,y5,y1,y2])

fig = plt.figure()
ax = fig.add_subplot(111,aspect='equal')
ax = plt.gca()
ax.cla()
ax.set_xlim((-R-0.3,R+0.3))
ax.set_ylim((-R-0.3,R+0.3))
ax.plot(X,Y,'ro')
circle=plt.Circle((0,0),1,fill=False)
fig.gca().add_artist(circle)
plt.title("Langrangian Points")
plt.annotate('L1', xy = (R1, y123), xytext = (R1+0.2, y123+0.2), textcoords = 'offset points')
plt.annotate('L2', xy = (R2, y123), xytext = (R2+0.2, y123+0.2), textcoords = 'offset points')
plt.annotate('L3', xy = (R3, y123), xytext = (R3+0.2, y123+0.2), textcoords = 'offset points')
plt.annotate('L4', xy = (x45, y4), xytext = (x45+0.2, y4+0.2), textcoords = 'offset points')
plt.annotate('L5', xy = (x45, y5), xytext = (x45+0.2, y5+0.2), textcoords = 'offset points')
plt.annotate('Sun', xy = (x1, y1), xytext = (x1+0.2, y1+0.2), textcoords = 'offset points')
plt.annotate('Earth', xy = (x2, y2), xytext = (x2+0.2, y2+0.2), textcoords = 'offset points')
plt.savefig("Lagrangian_Points")
plt.show()
