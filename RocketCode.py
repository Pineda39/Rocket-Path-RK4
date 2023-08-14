# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 20:14:16 2022

@author: Francisco
"""

from pylab import *
from numpy import *

""" Exercise 1 """

print('----- Exercise 1 ---------\n')

""" Firstly, we remind the constants of the problem"""

g = 9.81 # Gravitational constant
c = 0.2 # Dimensionless drag coefficient
p = 1.29 # Air density
s = 0.25 # Mean cross-sectional area of ​​the rocket
M = 7.5 # Mass of the rocket (not counting the mass of the fuel)
k = 0.01 # Constant of proportionality fuel consumption/thrust
T0 = 50 # Rocket thrust force


v0 = 50 # Initial speed
o0 = pi/4 # Initial angle
m0 = 7.5 # Initial mass of the fuel

def T(m,T0): # Function representing thrust force of the rocket
    if m > 0:
        return T0
    else:
        return 0

""" Function that solves, applying the RK4 method, the Cauchy problem of the evolution of a rocket with an engine.
 It takes as input a series of parameters corresponding to the conditions of the problem, making it possible to
 general implementation.

 Faced with the challenge of not knowing the upper end of the interval in which we want to solve the problem, we design an algorithm
 which iterates the Runge-Kutta 4 method until a stop condition is met"""

def rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T): 
    def f(t,z): # Function that defines the differential system
        (x,y,v,o,m) = z
        f1 = v*cos(o)
        f2 = v * sin(o)
        f3 = 1/(M+m)*(T(m,T0)-0.5*c*p*s*v**2)-g*sin(o)+1/(M+m)*k*T(m,T0)
        f4 = -g/v*cos(o)
        f5 = -k*T(m,T0)
        return array([f1,f2,f3,f4,f5])
    
    a = 0 # Lower extreme of the interval (corresponding to time)
    h = 0.01 # Mesh pitch (initially fixed)
    z0 = array([0,0,v0,o0,m0]) # initial values of the problem
    
    t = zeros(1) # zero array that will start the algorithm
    z = zeros((5,2)) 
    z[:,0]= z0
    i = 1 # Iteration counter
    
    """ The reason for creating z as a two-column array, instead of a single column made up of the initial values,
    is the convenience that it gives us when applying the Runge-Kutta 4 method in a general way, without an initial case
    to work with a single-column array with unique index. The existence of the extra column
    it will be corrected later in a very simple way"""

    while z[1,i-1] > 0 or i==1 : # Loop that will repeat as long as the y coordinate is positive (we also add the condition i==1 so that the first iteration occurs)
    # The coefficients of the rk4 method are calculated
        k1 = f(t[i-1],z[:,i-1]) 
        k2 = f(t[i-1]+h/2,z[:,i-1]+k1*h/2)
        k3 = f(t[i-1]+h/2,z[:,i-1]+k2*h/2)
        k4 = f(t[i-1]+h,z[:,i-1]+h*k3)
        
        # Create new empty helper arrays to extend t and z with the values ​​computed in the iteration of the rk4 method
        zi = zeros((5,i+1))
        ti= zeros(i+1)
        ti[:i]= t
        ti[i] = t[i-1]+h
        
        # If we are not in the first iteration (i>1), zi is filled in intuitively, adding the value obtained in the iteration of the numerical method
        if i > 1:
            zi[:,:i] = z
            
            zi[:,i] = zi[:,i-1]+ (h/6) * (k1+2*k2+2*k3+k4)
           
        # In the first iteration, we get rid of the extra column of z mentioned above
        else:
            zi[:,0] = z[:,0]
            zi[:,i] = zi[:,i-1]+ (h/6) * (k1+2*k2+2*k3+k4)
        
        """ To give the algorithm more precision regarding the point of impact of the rocket, we add some extra conditions:
        If in the current iteration a negative value of the "y" coordinate is obtained, and in the previous iteration
        a value was obtained that differs from the ground by more than 1mm, the current iteration is not valid, that is, the
        solution and it will continue iterating with a mesh step corresponding to half of the previous
        
        Thus, we make sure that the last point of our solution is, quite approximately, a point on the ground"""
        
        if zi[1,i] < 0 and zi[1,i-1]>0.001:
            h = h/2
        
        # If the previous conditions are not verified, the solution is modified with the value obtained in this iteration
        else:
            z = zi
            t = ti
            i += 1
        
    return (t[:i-1], z[:,:i-1]) # The function returns an array consisting of t and z, where z is an array corresponding to (x,y,v,o,m)
    
# Note that the solution does not include the last column, since it does not check the conditions that were requested in the While loop (In particular, the "y" coordinate is negative)

ymax = max(rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)[1][1,:]) # Maximum height
alcance = max(rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)[1][0,:]) # Rocket scope
print('Maximum height of the rocket is:',round(ymax,3), 'metres')
print('The scope of the rocket is:',round(alcance,3), 'metr')

(t,z) = rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)

# We draw the graphs of the speed and mass of fuel with respect to time
fig1=figure() 
plot (t , z [2,:]) 
grid(True)
xlabel('t (s)')
ylabel('v (m/s)')
legend(['Speed'])
title("Evolution of the speed")

fig11 = figure()
plot(t , z [4 ,:])
grid(True)
xlabel('t (s)')
ylabel('m (kg)')
legend(['Mass of the fuel'])
title("Evolution of the mass of the fuel")

fig2=figure()
plot ( z [0 ,:] , z [1 ,:]) # We draw the trajectory of the rocket with engine
grid(True)
xlabel('x')
ylabel('y')
title("Exercise 1: Trajectory")
legend(['Trajectory'])

""" Exercise 2 """

print('\n----- Exercise 2 ---------\n')

# Item a)
T0 = 0 
c= 0
(t1,z1) = rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('Item a)\n')
ymax1 = max(z1[1,:])
alcance1 = max(z1[0,:])
print('Maximum height of the rocket is:',round(ymax1,3), 'metres')
print('The scope of the rocket is:',round(alcance1,3), 'metres')

# Item b)
T0 = 0
c= 0.2
(t2,z2) = rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\n Item b)\n')
ymax2 = max(z2[1,:])
alcance2 = max(z2[0,:])
print('Maximum height of the rocket is:',round(ymax2,3), 'metres')
print('The scope of the rocket is:',round(alcance2,3), 'metres')

# Item c)
T0 = 10
c= 0.2
(t3,z3) = rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\n Item c)\n')
ymax3 = max(z3[1,:])
alcance3 = max(z3[0,:])
print('Maximum height of the rocket is:',round(ymax3,3), 'metres')
print('The scope of the rocket is:',round(alcance3,3), 'metres')

# Item d)
T0 = 50
c= 0
(t4,z4) = rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\n Item d)\n')
ymax4 = max(z4[1,:])
alcance4 = max(z4[0,:])
print('Maximum height of the rocket is:',round(ymax4,3), 'metres')
print('The scope of the rocket is:',round(alcance4,3), 'metres')

# Item e)
# The solution of this section coincides with that of exercise 1
print('\n Item e)\n')
ymax = max(z[1,:]) 
alcance = max(z[0,:]) 
print('Maximum height of the rocket is:',round(ymax,3), 'metres')
print('The scope of the rocket is:',round(alcance,3), 'metres')

fig3=figure() # We compare the trajectories in the same graph

plot (z1[0,:] , z1 [1 ,:], z2[0,:] , z2 [1 ,:], z3[0,:] , z3 [1 ,:],z4[0,:], z4 [1 ,:], z[0,:] , z[1 ,:])
grid(True)
xlabel('x')
ylabel('y')
title("Exercise 2: Trajectories")
legend(['Trajectory A', 'Trajectory B','Trajectory C','Trajectory D','Trajectory E'])
    

""" Exercise 3 """

print('\n----- Exercise 3 ---------\n')

o0 = pi/2.1

# Item a)
T0 = 0
c= 0
(t1,z1) = rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('Item a)\n')
ymax1 = max(z1[1,:])
alcance1 = max(z1[0,:])
print('Maximum height of the rocket is:',round(ymax1,3), 'metres')
print('The scope of the rocket is:',round(alcance1,3), 'metres')

# Item b)
T0 = 0
c= 0.2
(t2,z2) = rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\nItem b)\n')
ymax2 = max(z2[1,:])
alcance2 = max(z2[0,:])
print('Maximum height of the rocket is:',round(ymax2,3), 'metres')
print('The scope of the rocket is:',round(alcance2,3), 'metres')

# Item c)
T0 = 10
c= 0.2
(t3,z3) = rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\nItem c)\n')
ymax3 = max(z3[1,:])
alcance3 = max(z3[0,:])
print('Maximum height of the rocket is:',round(ymax3,3), 'metres')
print('The scope of the rocket is:',round(alcance3,3), 'metres')

# Item d)
T0 = 50
c= 0
(t4,z4) = rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\nItem d)\n')
ymax4 = max(z4[1,:])
alcance4 = max(z4[0,:])
print('Maximum height of the rocket is:',round(ymax4,3), 'metres')
print('The scope of the rocket is:',round(alcance4,3), 'metres')

# Item e)
T0 = 50
c= 0.2
(t5,z5) = rocket_model(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\nItem e)\n')
ymax5 = max(z5[1,:])
alcance5 = max(z5[0,:])
print('Maximum height of the rocket is:',round(ymax5,3), 'metres')
print('The scope of the rocket is:',round(alcance5,3), 'metres')


fig4=figure() # We compare the trajectories in the same graph
plot (z1[0,:] , z1 [1 ,:], z2[0,:] , z2 [1 ,:], z3[0,:] , z3 [1 ,:],z4[0,:], z4 [1 ,:],z5[0,:], z5[1 ,:])
xlabel('x')
ylabel('y')
title("Exercise 3: Trayectories")
grid(True)
legend(['Trayectory A', 'Trayectory B','Trayectory C','Trayectory D','Trayectory E' ])




