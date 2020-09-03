import math
#funcion define la funcion dada
def f(x):
    return math.cos(2*x)*2 - x*2

#funcion posfalsa
def posfalsa(x0,x1,e):
    paso = 1
    condicion = True
    while condicion:
        #formula ecuacion de la recta
        x2 = x0 - (x1-x0) * f(x0)/( f(x1) - f(x0) )
        print('Iteracion-%d, x2 = %0.32f and f(x2) = %0.33f' % (paso, x2, f(x2)))
#valida si son positivos o negativos y se iguala a x2 si es el caso
        if f(x0) * f(x2) < 0:
            x1 = x2
            #si es positivo x0=x2
        else:
            x0 = x2

        paso = paso + 1
        condicion = abs(f(x2)) > e

    print('\nRaiz requerida es: %0.33f' % x2)



x0 = input('Primer intervalo ')
x1 = input('Segundo Intervalo ')
e = input('Error torelable: ')

# casteo a float
x0 = float(x0)
x1 = float(x1)
e = float(e)


if f(x0) * f(x1) > 0.0:
    print('Error.')
    print('Intente denuevo con otros valores.')
else:
    posfalsa(x0,x1,e)

y = lambda x : math.cos(2*x)*2 - x*2
print('%0.38f' % y(0.514933266590332561207787875900976))

    #Error 10e8: 0.00000001
    #Error 10e16: 0.0000000000000001
    #Error 10e32: 0.00000000000000000000000000000001