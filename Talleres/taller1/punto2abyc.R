
binarioEnt<-function(numero) {
  sum<-0
  exp<-1
  while (numero > 0) {
    digito<-numero %% 2
    numero<-floor(numero / 2)
    sum<-sum + digito * exp
    exp<-exp * 10
  }
  return(sum)
}
numero<-readline("Numero decimal ?: ")
numero<-as.numeric(numero)


binarioDec <- function (decimal) {
  
  bits <- 15
  round (decimal , 1)
  binario <- c()
  
  for (i in 1 : bits) {
    
    decimal = decimal * 2;
    
    
    if(decimal >= 1) {
      
      decimal = decimal - 1
      binario[i] = 1
      
    }
    
    else
    {
      binario[i] = 0
    }
  }
  return(binario)
}

#Impresiones 

cat("La represenacion de pi ","En binario es :",binarioEnt(3),".",binarioDec(0.141592))
cat ("En decimal da :", 3.141592)
cat("El  punto c da ","En binario :",binarioEnt(11),".",binarioDec(0.25))
cat ("En decimal da :", 11.25)
cat("El  punto c da ","En binario :",binarioEnt(0),".",binarioDec(0.6666666))
cat ( "En decimal da :",2/3)
cat("El  punto c","En binario da :",binarioEnt(30),".",binarioDec(0.6))
cat ( "En decimal da :",30.6)
cat("El punto c  ","En binario da:",binarioEnt(99),".",binarioDec(0.9))
cat ("En decimal da :",99.9)
