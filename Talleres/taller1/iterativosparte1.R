f = function(x){log(x+2)}
g = function(x){sin(x)}
Tolerancia = 0.0000000000000001
q = function(x){log(x+2)-sin(x)}
Teorico = uniroot(q, c(-1.8,-1), tol = 1e-8) 

plot(f, xlim = c(-4,4), ylim=c(-4,4), ylab = "f(x)", col = "green")
par(new = TRUE)


plot(g, xlim = c(-4,4), ylim=c(-4,4), ylab = "g(x)", col = "blue")
plot(q, xlim = c(-4,4), ylim=c(-4,4), ylab = "f(x) - g(x)", col = "brown")
par(new = TRUE)

abline(0,0, col = "red")

#Variables
a = -1.9
b = -1
arr = c(0)
arr[1] = a
arr[2] = b
i=3
Erroresm1 = c() 
Erroresm1_1 = c()

#Inicio del while
while(abs(Teorico$root-arr[i-2]) > Tolerancia){
  Erroresm1[i-2]=abs(Teorico$root-arr[i-2])
  arr[i]= arr[i-1]-((q(arr[i-1])*(arr[i-1]-arr[i-2]))/(q(arr[i-1])-q(arr[i-2])))
  i = i+1
}

i=1
while(i<length(Erroresm1)+1){
  Erroresm1_1[i]=Erroresm1[i+1]
  i = i+1
}
plot(Erroresm1,Erroresm1_1, type="l")

numero_iteraciones =length(arr)

Metodo2 = c(0)
Erroresm2 = c(0)
Erroresm2_1 = c(0)
Metodo2[1] = a
Metodo2[2] = b
i=3

while(abs(Teorico$root-Metodo2[i-2]) > Tolerancia){
  Erroresm2[i-2]=abs(Teorico$root-arr[i-2])
  Metodo2[i]= Metodo2[i-1]-((q(Metodo2[i-1]))*(Metodo2[i-1]-Metodo2[i-2])/(q(Metodo2[i-1])-q(Metodo2[i-2])))
  i = i+1
}

i=1
while(i<length(Erroresm2)+1){
  Erroresm2_1[i]=Erroresm2[i+1]
  i = i+1
}

plot(Erroresm2,Erroresm2_1, type ="l",col ="blue")

tabla = data.frame(Erroresm1,Erroresm1_1,Erroresm2,Erroresm2_1)
tabla