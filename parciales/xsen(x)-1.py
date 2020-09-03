import numpy as np
import matplotlib.pyplot as plt

#funcion define la funcion dada
def f(x):
    return x * np.sin(x) - 1


puntos = []
itera = []

#funcion posfalsa
def posfalsa(x0,x1,e):
    paso = 1
    condicion = True
    while condicion:
        #formula ecuacion de la recta
        x2 = x0 - (x1-x0) * f(x0)/( f(x1) - f(x0) )
        print('Iteracion0-%d, x2 = %0.6f and f(x2) = %0.6f' % (paso, x2, f(x2)))
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

# casteo a float
#x0 = float(x0)
#x1 = float(x1)
#e = float(e)


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

tolerancias = [10**-6, 10**-8, 10**-15]   

for k in tolerancias: 
    
    puntos.clear()
    itera.clear()
    puntosBise.clear()
    iteraBise.clear()
    cont = 0 
    if f(-1) * f(2) > 0.0:
        print('Error.')
        print('Intente denuevo con otros valores.')
    else:
        posfalsa(-1,2,k)
    
    plt.plot(itera, puntos, 'r');
    plt.title("Posición falsa " + str(k))
    plt.xlabel("Iteraciones")
    plt.ylabel("Resultados")
    plt.show()

    biseccion(f, -1, 2, k)
    plt.plot(iteraBise, puntosBise, 'm')
    plt.xlabel("Iteraciones")
    plt.ylabel("Resultados")
    plt.title("Bisección " + str(k))
    plt.show()