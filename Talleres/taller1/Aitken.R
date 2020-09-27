#################################
#Problema Convergencia Acelerada
################################
#install.packages("rootSolve")
require("rootSolve")
###############
#Variables
###############
#Método de biseccion
a  = 0.1;
b = 1.6;
tolerancia = 0.00000001;
contador = 1;
aiken = c(0)
i=1
resultados = c (0)
valor = (a + b)/2
error = c(0)
errormas1= c(0)
iterador=1


#################
#Funcion
#################

f = function(contador){return(cos(1/contador))}

#Evaluamos la funcion en 1

cat("La funcion evaluada en el punto 1 da : ",f(1));
plot(f, xlim = c(-1,1), ylim=c(-1,1), ylab = "f(x)", col = "green")
abline(1,0, col = "red")

############
#funciones 2
############

#for
Vectorg =c (0)
vectorh = c (0)
g = function(x)3*(sin(x))^3-1
h = function(x)4*sin(x)*cos(x)
vectorh = (uniroot.all(h,c(1,20)))
vectorg =(uniroot.all(g,c(1,20)))
vectorg
vectorh
i = 1
j = 1
n = length(vectorg)
k = length(vectorh)
vectorsolu =c (0)
contador=1

for (i in n) {
  
  for (j in k) {
    
    if (vectorg[i]== vectorh[j])
    {
      vectorsolu[contador] =  vectorg[i];
      contador = contador +1;
    }
    j = j +1;
  }
  i = i +1
}

vectorg
###############
vectorh
vectorsolu

############
#proceso
############

while((b - a) > tolerancia){
  valor = (a + b)/2;
  resultados[contador]= valor
  error[contador] = (abs(b-a))/2
  if (valor != 0 & a != 0){
    
    
    if((f(valor)*f(a)) <0){
      
      b = valor;
    }else{
      a = valor;
    }
    
  }
  contador = contador + 1;
}

while(iterador<contador){
  errormas1[iterador]=error[iterador+1]
  iterador = iterador+1
}

#Aiken
while(i<contador){
  aiken[i]= resultados[i+2]-(((resultados[i+2]-resultados[i+1])^2)/(resultados[i+2]-(2*resultados[i+1])+resultados[i]))
  i = i+1;
}

###############
#Creacion Tabla
###############
Datos_Convergencia = data.frame(error, errormas1,resultados,aiken )
Datos_Convergencia

plot(error,errormas1,type ="l", col ="green")