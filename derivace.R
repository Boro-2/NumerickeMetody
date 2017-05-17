fce<-function(x){return(x*x+1)} # Integral from 1 to 10 is 342
fce2<-function(x){return(2*x+5)} # Integral from 0 to 5 is 50
fce3<-function(x){return(x*x*x-3*x*x-3*x+20)} # Integral from -2 to 5 is 127,75
fce4<-function(x){return(5^(2*x))} # Integral from -2 to 5 is 3 033 862
fce5<-function(x){return(abs(2*cos(x^2)*x^2)+sin(2*x)+exp(x)/10+10)} ## romberg to pocita presne
a<--10.1
b<--9.9
x<-c(sort(runif(10000,a,b)))
y<-c()
for(i in 1:length(x)){y[i]<-fce5(x[i])}



plot(x,y,xlim = c(min(x)-1,max(x)+1),ylim = c(min(y)-2,max(y)+2))

derivace<-function(f,x0,h,derivace){#ze dvou stran
  if(derivace==1)
  {v<-(f(x0+h)-f(x0-h))/(2*h)}
  else if(derivace==2){v<-(f(x0+h)-2*f(x0)+f(x0-h))/(h*h)}
  print(v)
  return(v)
}

richardson<-function(f,x0,h0,iter,derivace){
  mat<-matrix(nrow = iter,ncol = iter)
  initH=h0
  h<-initH
  parts<-1
  for(i in 1:iter){
    mat[i,1]<-derivace(f,x0,h,derivace)
    h<-h/2
  }
  print(mat[,1])
  for(j in 2:iter){
    for(i in 1:(iter+1-j)){
      mat[i,j]<-(4^(j-1)*mat[i+1,j-1]-mat[i,j-1])/(4^(j-1)-1)
    }
  }
  print(mat[,1])
  print(mat[1,])
  return(mat[1,iter])
}

derivace(fce,0,0.01,1)
h=2
#richardson(fce3,-10,h,10,2)

v<- c(1,2,3)
k<- c(1,2,10)
l<- v %*% k
print(l)