# Polinomio: 2X⁴-3X²+3X-4 -> X = -2

def horner(x):
    coef = [1, 0, -3, 2, -4]
    res = 0

    for i in range(len(coef)):
        res = res * x + coef[i]

    print("El resultado del polinomio evaluado en X= ", x, " es igual a: ", res)


if __name__ == "__main__":
    horner(-2)