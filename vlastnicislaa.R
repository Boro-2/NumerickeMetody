ml<-matrix(c(15,4,8,3,2,3,5,3,-1),nrow=3, byrow = TRUE)
mhl<-matrix(c(2.2,1,0.5,2,1,1.3,2,1,0.5,2,0.5,1.6,2,1,1.6,2),4, byrow = TRUE)
print(mhl)
library("plotrix")
gersgerin<-function(m){
  s<-nrow(m)
  nums<-c()
  radiuses<-c()
  
  for(i in 1:s){
    sum<-0
    for(j in 1:s){
      if(j==i) { nums[i]<-m[i,i] }
      else{ sum<-sum+abs(m[i,j])}
    }
    radiuses[i]<-sum
  }
  print(nums)
  print(radiuses)
  x<-c(min(nums),max(nums))
  y<-c(-max(radiuses),max(radiuses))
  plot(x,y, asp = 1, xlim = c(min(nums)-max(radiuses),max(nums)+max(radiuses)))
  abline(h=0, lwd=1,col='black')
  for(i in 1:length(nums)){
    points(nums[i],0,type='p',col='red', bg='red',pch=21, cex=1)
    draw.circle(nums[i], 0, radiuses[i], nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)
  }
  
}

nejmensiC<-function(m){
  r<-nrow(m)
  n<-ncol(m)
  res<-c()
  for(i in 1:r){
    sum<-0
    for(j in 1:n){
      sum<-sum+m[i,j]#*m[i,j]
    }
    res[i]<-sqrt(sum)
  }
  sum<-1
  for(k in 1:length(res)){sum<-sum*res[k]}
  print(sum^(1/r))
}

nasobmv<-function(m,v){
  b<-c()
  n<-nrow(m)
  for(i in 1:n){
    sum<-0
    for(j in 1:n){
      sum<-sum+m[i,j]*v[j]
    }
    b[i]<-sum
  }
  return(t(b))
}

krylovaM<-function(A){
  y<- t(c(1,0,0,0))
  s<-matrix(y)
  n<-nrow(A)
  for(i in 1:(n-1)){
    s1<-nasobmv(A,s[,i])##nekde tu bude chyba
    s<-cbind(s,t(s1))
  }
  print(s)
  p<-t(nasobmv(-A,s[,n]))
  #print(p)
  re<-solve(s,p)##tohle jsou koeficienty polynomu x^4+re[1]*x^3+re[2]*x^2...
  print(re) # udelam funkci a jeji derivaci a pocitam koreny
  return(re)
  
}


fcehl<-function(x){
  return(x^4-6*x*x*x-0.2*x*x+12.7350*x-2.7616)
}
fceder<-function(x){
  return(4*x*x*x-18*x*x+12.7350)
}

skalar<-function(x,y){
  sum<-0
  if(length(x)==length(y)){
    for(i in 1:length(x))
    {
      sum<-sum+x[i]*y[i]
    }
    return(sum)
  }
  else print("error")
  return(0)
  
}

norma <- function(x){
  n<-length(x)
  sum<-0
  for(i in 1:n){
    sum<-sum+ x[i]*x[i]
  }
  return(sum)
}

#ROzpocet dominnatniho vlastniho cisla
rozpocetDom<-function(m){
  v0<-t(c(1,0,0,0))
  y<-t(c(1,0,0,0))
  a<- vector()
  e<- 0.0001
  v<-nasobmv(m,v0)
  k<-1
  a1<- skalar(y,v)/skalar(y,v0)
  a[k]<- a1
  while(norma(nasobmv(m,v)-a[k]*v)>=e){
    v0<-v
    v<-nasobmv(m,v)
    k<-k+1
    a[k]=skalar(y,v)/skalar(y,v0)
  }
  amax<-a[length(a)]
  print(amax)
  return(amax)
  plot(a, pch=20,type = 'b')
}

posunutiSpektra<-function(m){
  jed<-matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,byrow = TRUE)#jednotkova mat
  print(jed)
  amax<-rozpocetDom(m)
  A1<-m-(amax*I)
  print(A1)
  amin<-amax+rozpocetDom(A1)
  print(amin)
}

#gersgerin(mhl)
#nejmensiC(mhl)
#krylovaM(mhl)
x<--2
for(i in 1:100){
  print(NewtonRoot(x,fcehl,fceder)) #nahodne volim x abych zjistil vsechny vlastni cisla
  x<-x+0.1
}
#posunutiSpektra(mhl)

