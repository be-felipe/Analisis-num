################
#Horner reto 1.1
################

#Ecuacion: x^4-5*x^3-9*x^2+155*x-250

#Variables 
n=5     #Numero de grados

coeficientes=c(1,-5,-9,155,-250)      #Coeficientes del polinomio     
coeficientesDerivada=c()      #Coeficientes priemera derivada del polinomio
coeficientesDerivada2=c()     #Coeficientes segunda derivada del polinomio
resultado1=c()
resultado2=c()
res=0
j=n-1
x=3     #Valor x=3




##1°Derivada##
i=1;
while (i<n){
  coeficientesDerivada[i]=coeficientes[i]*(j)
  i=i+1
  j=j-1
}

##Horner Primera Derivada##
j=1
while(j<n){
  res=res*x+coeficientesDerivada[j]
  resultado1[j]= res;
  j=j+1;
}

res=0

##2°Derivada##
j=n-2
i=1;
while (i<n){
  coeficientesDerivada2[i]=coeficientesDerivada[i]*(j)
  i=i+1
  j=j-1
}


##Horner Segunda Derivada##
j=1
while(j<n-1){
  res=res*x+coeficientesDerivada2[j]
  resultado2[j]= res;
  j=j+1;
}


TablaHornerPrimeraDerivada = data.frame(resultado1) 
TablaHornerSegundaDerivada = data.frame(resultado2) 

TablaHornerPrimeraDerivada
TablaHornerSegundaDerivada



expresion <- expression (coeficientes[1]x^4 + coeficientes[2]*x^3 + coeficientes[3]*x^2 + coeficientes[4]*x + coeficientes[5]) # escribimos el polinomio
derivada <- D(expresion, "x") # Derivada del polinomio
valores=c(5, -5, 2+6i, 2-6i) 

for (i in valores) {
  x <- 10**-8  # Cualquier valor diferente de aprox
  aprox <- i # valor puntoinicial
  
  while ( x != aprox) {
    
    x <- aprox # Se le asigna el valor aproximado a x.
    
    reemplazoexpresion <- eval(expresion) #Reemplaza el valor de x en "expresión"
    
    reemplazoderiv <- eval(derivada) #Reemplaza el valor de x en "derivada"
    
    #newton
    
    aprox <- x - (reemplazoexpresion/reemplazoderiv) #Ecuación método de Newton
  }
  print("Raiz")
  print(x)
}



