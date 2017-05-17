readkey <- function()
{
  cat("[press [enter] to continue]")
  number <- scan(n=1)
}
#kvadraticka rovnice najjit oba koreny pomoci 3 metod
#plocha pod osou x napr 1NewtonCoatsuv, 1 z Gaussovy kvadratury
#prolozit navzorkovanou fci kubickymi spliny
x<-sort(runif(400,-200,10))
y<-c()
for(i in 1:length(x))
{ y[i]<-f(x[i])}
#print(y)

plot(x,y,type='b',xlim = c(-300,20),ylim = c(-400,500))#
abline(h=0, col='gray2', lwd=2)



NewtonRoot<-function(x,origf,derf){
  
  x1<-x
  y<-origf(x1)
  yd1<-derf(x1)
  x2<-x1-y/yd1
  
  while(abs(x2-x1)>0.000001){
    x1<-x2
    y<-origf(x1)
    yd1<-derf(x1)
    x2<-x1-y/yd1
    #print(x2)
    
    c<-y-yd1*x1
    prusecikX<-(-c)/yd1
    
   # abline(a=c,b=yd1, col='blue3', lwd=1.5) # y= a + bx
    #abline(v=prusecikX, col='maroon', lwd=1.5) 
    #points(x1,y,type='p',col='red', bg='red',pch=21, cex=1) #bg barva pozadi, pch tvar bodu, cex velikost bodu
    #readkey() #zmackni cislo, jedu dal, jinak koncim
  }
  #print("hotovo")
  #print(x2)
  return(x2)
}

goldenSection <- function(f, lower.bound, upper.bound){
  golden.ratio = 2/(sqrt(5) + 1)
  f1 = f(lower.bound)
  f2 = f(upper.bound)
  
  iteration = 0
  while (abs(upper.bound - lower.bound) > 0.0001)
  {
    iteration = iteration + 1
    if (f2 < f1)
    {
        x2<-lower.bound + golden.ratio*(upper.bound - lower.bound)
        f2<-f(x2)
        upper.bound<-x2
    } 
    else 
    {
        x1<-upper.bound - golden.ratio*(upper.bound - lower.bound)
        f1<-f(x1)
        lower.bound<-x1
    }
  }
  
  maxim = (lower.bound + upper.bound)/2
  cat('max =', maxim, '\n')
}

bisection<-function(f, x1, x2){
  citac<-0
  plot(f, x1, x2, col=2)
  abline(v=x1, col=3)
  abline(v=x2, col=3)
  abline(h=0, lwd=2)
  if(f(x1)*f(x2)<0)
  {
    repeat{
      citac<-citac+1
      c<-(x1+x2)/2
      abline(v=c, lwd=1, col=4)
      if(c==x1 || c==x2) 
        break
      
      if(f(x1)*f(c)<0) 
        x2<-c
      else 
        x1<-c
    }
    return(c)
  }
  else print("error")
}

secny<-function(f,x1,x2){
  
  repeat{
    x3<-x1-(f(x1)*(x2-x1)/(f(x2)-f(x1)))
    if(abs(x3-x2)<0.001)
      break
    x1<-x2
    x2<-x3
  }
  return(x3)
}

regulafalsi<-function(f, x1, x2){
  citac<-0
  plot(f, x1, x2, col=2)
  abline(v=x1, col=3)
  abline(v=x2, col=3)
  abline(h=0, lwd=2)
  lines(c(x1,x2),c(f(x1),f(x2)), col="blue")
  if(f(x1)*f(x2)<0)
  {
    
    repeat{
      citac<-citac+1
      c<-x1-(f(x1)*(x2-x1)/(f(x2)-f(x1)))
      lines(c(x1,x2),c(f(x1),f(x2)), col="blue")
      if(c==x1 || c==x2) 
        break
      
      if(f(x1)*f(c)<0) 
        x2<-c
      else 
        x1<-c
    }
    return(c)
  }
  else print("error")
}
jednoduchaiterace<-function(f, x){
  plot(f, col=2, xlim=c(-5,5), ylim=c(-5,5))
  abline(h=0,col="black")
  abline(a=0,b=1,col="3")
  repeat{
    x1<-f(x)
    if(abs(x1-x)<0.0001)
      break
    x<-x1
  }
  return(x1)
}


#

f<-function(x){return((x*x)-3*x-4)}#puvodni fce koreny -1,4
fd<-function(x){return(2*x-3)}#jeji derivace
x1<-NewtonRoot(50,f,fd)
x2<-NewtonRoot(-7,f,fd)
print(x1)
print(x2)
x1<-bisection(f,0,10)
x2<-bisection(f,0,-5)
print(x1)
print(x2)
x1<-secny(f,0,6)
x2<-secny(f,0,-3)
print(x1)
print(x2)
x<-c(seq(-1,4,(5/10)))
y<-c(f(x)) 
res<-simpsonRule(x,y,-1.0001,4.0001)
print(res)
res2<-GaussIntegral3(x,y,-1,4)
res3<-GaussIntegral3Orig(f,-1,4)
print(res3)

plot(x,y, type = "l",xlim = c(min(x)-3,max(x)+3),ylim = c(min(y)-1,max(y)+1))
x0s <- c(seq(-3,7,(10/1000)))
#print(x0s)
y0s<-kubSplVectorx0(x,y,x0s)
#print(y0s)
for(i in 1:length(x0s)){
  points(x0s[i], y0s[i],type='p',col='red',bg='red',pch=21,cex=0.3)
}
