# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 22:35:54 2020

@author: USA
"""

import cmath

def poly(x):
    return x**4 -5*x**3 - 9*x**2 + 155*x -250
def g(x):
    return 4*x*3 - 15*x*2 - 18*x + 155
def h(x):
    return 12*x**2 - 30*x - 18

def laguerre(poly, g, h, start, epsilon):

    p_d, p_d_, n = g, h, 4
    start = complex(start)
    while True:
        px = poly(start)
        if not px:
            return start
        g = p_d(start) / px
        h = g ** 2 - p_d_(start) / px
        dp = cmath.sqrt((n - 1) * (n * h - g**2))
        d1 = g + dp
        d2 = g - dp
        if abs(d2) > abs(d1):
            d = d2
        else:
            d = d1
        a = n / d
        x_n = start - a
        if str(start) == str(x_n) or abs(start - x_n) < epsilon:
            return start
        start = x_n


print(laguerre(poly, g, h, 6+8j , 10**-8 ))
