import math
import matplotlib.pyplot as plt
import numpy as np

#funcion define la funcion dada
def f(x):
    return np.float128(x**3 - 2*x**2 + (4/3)*x - 8/27)


puntos = []
itera = []
x = np.arange(-1,2, 0.1)

#funcion posfalsa
def posfalsa(x0,x1,e):
    paso = 1
    condicion = True
    while condicion:
        #formula ecuacion de la recta
        x2 = np.float128(x0 - (x1-x0) * f(x0)/( f(x1) - f(x0) ))
        print('Iteracion-%d, x2 = %0.32f and f(x2) = %0.33f' % (paso, x2, f(x2)))
        #valida si son positivos o negativos y se iguala a x2 si es el caso
        if f(x0) * f(x2) < 0:
            x1 = x2
            #si es positivo x0=x2
        else:
            x0 = x2

        paso = paso + 1
        condicion = abs(f(x2)) > e
        puntos.append(x2)
        itera.append(paso)

    print('\nRaiz requerida es: %0.8f' % x2)



#x0 = input('Primer intervalo ')
#x1 = input('Segundo Intervalo ')
#e = input('Error torelable: ')

x0 = -1
x1 = 1
e = 10**-2

# casteo a float
x0 = float(x0)
x1 = float(x1)
e = float(e)



    

puntosBise = []
iteraBise = []

def biseccion(f, x0, x1, tol):
    cont = 0
    while x1-x0>=tol:
        cont = cont + 1
        x2=(x0+x1)/2
        if f(x2)==0:
            return (x2, cont)
        else:
            if f(x0)*f(x2)>0: 
                x0=x2
            else:
                x1=x2
        puntosBise.append(x2)
        iteraBise.append(cont)
    return (x2, cont)


funcion = []
tolerancias = [10**-4, 10**-6, 10**-8]
for i in x:
    funcion.append(f(i))
    

for k in tolerancias: 
    
    puntos.clear()
    itera.clear()
    puntosBise.clear()
    iteraBise.clear()
    if f(-1) * f(1) > 0.0:
        print('Error.')
        print('Intente denuevo con otros valores.')
    else:
        posfalsa(-1,1,k)
    
    plt.plot(itera, puntos, 'r');
    plt.title("Posición falsa " + str(k))
    plt.xlabel("Iteraciones")
    plt.ylabel("Resultados")
    plt.xlim(-200, 3000)
    plt.show()
    
    
    biseccion(f, -1, 1, k)
    plt.plot(iteraBise, puntosBise, 'm')
    plt.xlabel("Iteraciones")
    plt.ylabel("Resultados")
    plt.title("Bisección " + str(k))
    plt.show()


