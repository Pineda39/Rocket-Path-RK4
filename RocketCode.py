# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 20:14:16 2022

@author: Francisco
"""

from pylab import *
from numpy import *

""" Ejercicio 1 """

print('----- Ejercicio 1 ---------\n')

""" En primer lugar recogemos las distintas constantes que aparecen en el problema"""

g = 9.81 # Constante gravitatoria
c = 0.2 # Coeficiente adimensional de arrastre 
p = 1.29 # Densidad del aire
s = 0.25 # Área media de la sección tranversal del cohete
M = 7.5 # Masa del cohete (sin contar la masa del combustible)
k = 0.01 # Constante de proporcionalidad consumo de combustible/empuje
T0 = 50 # Fuerza de empuje del cohete


v0 = 50 # Velocidad inicial
o0 = pi/4 # Ángulo inicial
m0 = 7.5 # Masa inicial de fuel

def T(m,T0): # Función que representa la fuerza que ejerce el motor del cohete
    if m > 0:
        return T0
    else:
        return 0

""" Función que resuelve, aplicando el método RK4, el problema de Cauchy de la evolución de un cohete con motor.
 Toma como entrada una serie de parámetros correspondientes a las condiciones del problema, haciendo posible su 
 implementación de forma general.

 Ante el reto de no saber el extremo superior del intervalo en el que se desea resolver el problema, diseñamos un algoritmo
 que realiza iteraciones del método de Runge-Kutta 4 hasta que se cumple una condición de parada"""


def modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T): 
    def f(t,z): # Función que define el sistema diferencial
        (x,y,v,o,m) = z
        f1 = v*cos(o)
        f2 = v * sin(o)
        f3 = 1/(M+m)*(T(m,T0)-0.5*c*p*s*v**2)-g*sin(o)+1/(M+m)*k*T(m,T0)
        f4 = -g/v*cos(o)
        f5 = -k*T(m,T0)
        return array([f1,f2,f3,f4,f5])
    
    a = 0 # Extremo inferior del intervalo (correspondiente al tiempo)
    h = 0.01 # Paso de malla (fijo inicialmente)
    z0 = array([0,0,v0,o0,m0]) # valores iniciales del problema
    
    t = zeros(1) # Creamos arrays de zeros que iniciarán el algoritmo
    z = zeros((5,2)) 
    z[:,0]= z0
    i = 1 # Contador de las iteraciones
    
    """ La razón de crear z como un array de dos columnas, en vez de una única columna formada por los valores iniciales,
    es la comodidad que nos brinda a la hora de aplicar el método de Runge-Kutta 4 de forma general, sin un caso inicial
    para trabajar con un array de una sola columna con índice único. La existencia de la columna extra
    se subsanará más tarde de forma muy simple 
    """

    while z[1,i-1] > 0 or i==1: # Bucle que se repetirá mientras la coordenada y sea positiva(Añadimos además la condición i==1 para que se produzca la primera iteración)
        # Se calculan los coeficientes del método rk4
        k1 = f(t[i-1],z[:,i-1]) 
        k2 = f(t[i-1]+h/2,z[:,i-1]+k1*h/2)
        k3 = f(t[i-1]+h/2,z[:,i-1]+k2*h/2)
        k4 = f(t[i-1]+h,z[:,i-1]+h*k3)
        
        # Se crean nuevos arrays auxiliares vacíos para ampliar t y z con los valores calculados en la iteración del método rk4
        zi = zeros((5,i+1))
        ti= zeros(i+1)
        ti[:i]= t
        ti[i] = t[i-1]+h
        
        # Si no estamos en la primera iteración (i>1), se rellena zi de forma intuitiva, añadiendo el valor obtenido en la iteración del método numérico
        if i > 1:
            zi[:,:i] = z
            
            zi[:,i] = zi[:,i-1]+ (h/6) * (k1+2*k2+2*k3+k4)
           
        # En la primera iteración, nos deshacemos de la columna extra de z mencionada anteriormente
        else:
            zi[:,0] = z[:,0]
            zi[:,i] = zi[:,i-1]+ (h/6) * (k1+2*k2+2*k3+k4)
        
        """ Para otorgar al algoritmo de mayor precisión respecto al punto de impacto del cohete, añadimos unas condiciones extras:
        Si en la iteración actual se obtiene un valor negativo de la coordenada "y", y en la iteración anterior 
        se obtuvo un valor que dinsta del suelo más de un 1mm, la iteración actual no vale, es decir, no se modificará la 
        solución y se seguirá iterando con un paso de malla correspondiente a la mitad del anterior
        
        Así, nos aseguramos que el último punto de nuestra solución es, de forma bastante aproximada, un punto del suelo"""
        
        if zi[1,i] < 0 and zi[1,i-1]>0.001:
            h = h/2
        
        # Si no se verifican las condiciones anteriores, se modifica la solución con el valor obtenido en esta iteración
        else:
            z = zi
            t = ti
            i += 1
        
    return (t[:i-1], z[:,:i-1])  # La función devuelve un array formado por t y z, donde z es un array correspondiente a (x,y,v,o,m)
    
# Nótese que la solución no incluye la última columna, pues esta no verifica las condiciones que se pedían en el bucle While (En particular, la coordenada "y" es negativa)

ymax = max(modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)[1][1,:]) # Altura máxima
alcance = max(modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)[1][0,:]) # Alcance del cohete
print('La altura máxima del cohete es:',round(ymax,3), 'metros')
print('El alcance del cohete es:',round(alcance,3), 'metros')

(t,z) = modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)

# Dibujamos las gráficas de la velocidad y masa de fuel respecto del tiempo
fig1=figure() 
plot (t , z [2,:]) 
grid(True)
xlabel('t (s)')
ylabel('v (m/s)')
legend(['Velocidad'])
title("Evolución de la velocidad")

fig11 = figure()
plot(t , z [4 ,:])
grid(True)
xlabel('t (s)')
ylabel('m (kg)')
legend(['Masa de fuel'])
title("Evolución de la masa de fuel")

fig2=figure()
plot ( z [0 ,:] , z [1 ,:]) # Dibujamos la trayectoria del cohete con motor
grid(True)
xlabel('x')
ylabel('y')
title("Ejercicio 1: Trayectoria")
legend([' trayectoria '])

""" Ejercicio 2 """

print('\n----- Ejercicio 2 ---------\n')

# Apartado a)
T0 = 0 
c= 0
(t1,z1) = modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('Apartado a)\n')
ymax1 = max(z1[1,:])
alcance1 = max(z1[0,:])
print('La altura máxima del cohete es:',round(ymax1,3), 'metros')
print('El alcance del cohete es:',round(alcance1,3), 'metros')

# Apartado b)
T0 = 0
c= 0.2
(t2,z2) = modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\nApartado b)\n')
ymax2 = max(z2[1,:])
alcance2 = max(z2[0,:])
print('La altura máxima del cohete es:',round(ymax2,3), 'metros')
print('El alcance del cohete es:',round(alcance2,3), 'metros')

# Apartado c)
T0 = 10
c= 0.2
(t3,z3) = modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\nApartado c)\n')
ymax3 = max(z3[1,:])
alcance3 = max(z3[0,:])
print('La altura máxima del cohete es:',round(ymax3,3), 'metros')
print('El alcance del cohete es:',round(alcance3,3), 'metros')

# Apartado d)
T0 = 50
c= 0
(t4,z4) = modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\nApartado d)\n')
ymax4 = max(z4[1,:])
alcance4 = max(z4[0,:])
print('La altura máxima del cohete es:',round(ymax4,3), 'metros')
print('El alcance del cohete es:',round(alcance4,3), 'metros')

# Apartado e)
# La solución de este apartado coincide con la del ejercicio 1
print('\nApartado e)\n')
ymax = max(z[1,:]) 
alcance = max(z[0,:]) 
print('La altura máxima del cohete es:',round(ymax,3), 'metros')
print('El alcance del cohete es:',round(alcance,3), 'metros')

fig3=figure() # Comparamos las trayectorias en una misma gráfica

plot (z1[0,:] , z1 [1 ,:], z2[0,:] , z2 [1 ,:], z3[0,:] , z3 [1 ,:],z4[0,:], z4 [1 ,:], z[0,:] , z[1 ,:])
grid(True)
xlabel('x')
ylabel('y')
title("Ejercicio 2: Trayectorias")
legend(['Trayectoria A', 'Trayectoria B','Trayectoria C','Trayectoria D','Trayectoria E'])
    

""" Ejercicio 3 """

print('\n----- Ejercicio 3 ---------\n')

o0 = pi/2.1

# Apartado a)
T0 = 0
c= 0
(t1,z1) = modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('Apartado a)\n')
ymax1 = max(z1[1,:])
alcance1 = max(z1[0,:])
print('La altura máxima del cohete es:',round(ymax1,3), 'metros')
print('El alcance del cohete es:',round(alcance1,3), 'metros')

# Apartado b)
T0 = 0
c= 0.2
(t2,z2) = modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\nApartado b)\n')
ymax2 = max(z2[1,:])
alcance2 = max(z2[0,:])
print('La altura máxima del cohete es:',round(ymax2,3), 'metros')
print('El alcance del cohete es:',round(alcance2,3), 'metros')

# Apartado c)
T0 = 10
c= 0.2
(t3,z3) = modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\nApartado c)\n')
ymax3 = max(z3[1,:])
alcance3 = max(z3[0,:])
print('La altura máxima del cohete es:',round(ymax3,3), 'metros')
print('El alcance del cohete es:',round(alcance3,3), 'metros')

# Apartado d)
T0 = 50
c= 0
(t4,z4) = modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\nApartado d)\n')
ymax4 = max(z4[1,:])
alcance4 = max(z4[0,:])
print('La altura máxima del cohete es:',round(ymax4,3), 'metros')
print('El alcance del cohete es:',round(alcance4,3), 'metros')

# Apartado e)
T0 = 50
c= 0.2
(t5,z5) = modelo_cohete(v0, o0, m0, M, c, p, k, s, g, T0, T)

print('\nApartado e)\n')
ymax5 = max(z5[1,:])
alcance5 = max(z5[0,:])
print('La altura máxima del cohete es:',round(ymax5,3), 'metros')
print('El alcance del cohete es:',round(alcance5,3), 'metros')


fig4=figure() # Comparamos las trayectorias en una misma gráfica
plot (z1[0,:] , z1 [1 ,:], z2[0,:] , z2 [1 ,:], z3[0,:] , z3 [1 ,:],z4[0,:], z4 [1 ,:],z5[0,:], z5[1 ,:])
xlabel('x')
ylabel('y')
title("Ejercicio 3: Trayectorias")
grid(True)
legend(['Trayectoria A', 'Trayectoria B','Trayectoria C','Trayectoria D','Trayectoria E' ])




