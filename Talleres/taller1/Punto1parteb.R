resultadoIt = 0;
i = 0;
x = 1.0001;

repeat{
  
  resultadoIt=resultadoIt+x^i;
  i=i+1;
  if(i>50) break;
}


q=function(b) ((b^51)-1)/(b-1)
resultadoDirecto=q(x);


errorA=abs(resultadoDirecto-resultadoIt);


errorR=((resultadoIt/(resultadoIt+resultadoDirecto))*(errorA/resultadoIt))+((resultadoDirecto/(resultadoIt+resultadoDirecto))*(errorA/resultadoDirecto));


Datos = data.frame(i,errorA,errorR,resultadoIt,resultadoDirecto)
print(Datos)

