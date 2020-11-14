library(pracma)


temperaturas1 = read.csv(file = "Itatira.csv",header = TRUE,sep = ";")

temperaturas2 = read.csv(file = "Santa Quiteira.csv",header = TRUE, sep = ";")


x = seq(from = 1, to = 720, by = 1)
y = temperaturas1$Temp

diasIdeales = temperaturas1$Dia

horasIdeales = temperaturas1$Hora

indicesIdeales = x


ones = rep(1, 720)
eliminate = sample.int(720,720*0.2)
for (e in eliminate) {
  ones[e] = 0
}


newX = c()
newY = c()
i = 1
j = 1


for (o in ones)
{
  if(o == 1)
  {
    newX[i] = x[j]
    newY[i] = y[j]
    i = i + 1
  }
  j = j +1
}



plot(x,y,type='l', ylab = "Temperaturas", xlab = "Indices Ideales")

lines(spline(newX,newY,n=200),col=2)

interpolado = splinefun(newX,newY)

arregloInterpolados = c()
k = 1

error = c()

for(var in x)
{
  errorY = interpolado(var)
  
  auxi = abs((y[k] - errorY)/y[k])
  error = c(error, abs((y[k] - errorY)/y[k]))
  k = k + 1
}

max(error)
min(error)
mean(error)



arregloDeCalculos = c()

for (i in 1:length(temperaturas2$Dia)) 
{
  
  auxDia = temperaturas2$Dia[i]
  auxHora = temperaturas2$Hora[i]
  
  for(j in 1:720)
  {
    
    if((diasIdeales[j] == auxDia) && (horasIdeales[j] == auxHora))
    {
      arregloDeCalculos = c(arregloDeCalculos,indicesIdeales[j])
    }
  }
}

nuevosY = c()

errorNuevaEstacion = c()

z = 1
for (xd in arregloDeCalculos) 
{
  nuevosY = c(nuevosY, interpolado(xd))
  errorNuevaEstacion = c(errorNuevaEstacion, abs((temperaturas2$Temp[z] - nuevosY[z])/temperaturas2$Temp[z]))
  z = z + 1
}

plot(arregloDeCalculos,temperaturas2$Temp, ylab = "Temperaturas", xlab = "Indices Calculados", type = 'l')

lines(arregloDeCalculos, nuevosY, col = 3)



maximo = 0

media = 0

for (error in errorNuevaEstacion) {
  
  if(error > maximo)
    maximo = error
  
  media = media + error
  
}

minimo = 200

for (error in errorNuevaEstacion) {
  
  if(error < minimo)
    minimo = error
  
}

print(minimo)


qqnorm(errorNuevaEstacion)
qqline(errorNuevaEstacion)

print(maximo)