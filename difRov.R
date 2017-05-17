#Rung-kuttovy metody

difder<-function(x,y){
  return(-2*(y-x)+1)
}
fceex<-function(x){return(exp(-2*x)+x)}
plot(0,1,xlim = c(0,2),ylim=c(0,1))

analyticky<-function(x0,h,iter){
  x<-c()
  y<-c()
  x[1]<-x0
  y[1]<-fceex(x[1])
  for(i in 1:iter){
    x[i+1]<-x[i]+h
    y[i+1]<-fceex(x[i+1])
    points(x[i+1], y[i+1], type='p',col='green', bg='green',pch=21, cex=0.2)
  }
}

euler1r<-function(dd,x0,y0,h,iter){
  y<-c()
  x<-c()
  x[1]<-x0
  y[1]<-y0
  for(i in 1:iter){
    x[i+1]<-x[i]+h
    y[i+1]<-y[i]+h*dd(x[i],y[i])
    points(x[i+1], y[i+1], type = "p",col='blue')
  }
}

euler1rprumer<-function(dd,x0,y0,h,iter){
  y<-c()
  x<-c()
  x[1]<-x0
  y[1]<-y0
  for(i in 1:iter){
    x[i+1]<-x[i]+h
    fxy<-dd(x[i],y[i])
    y[i+1]<-y[i]+h*(fxy+dd(x[i+1],y[i]+h*fxy))
    points(x[i+1], y[i+1], type = "p",col='red3')
  }
}

euler1r4interval<-function(dd,x0,y0,h,iter){
  y<-c()
  x<-c()
  x[1]<-x0
  y[1]<-y0
  for(i in 1:iter){
    x[i+1]<-x[i]+h
    k1<-dd(x[i],y[i])
    k2<-dd(x[i]+h/2,y[i]+h*k1/2)
    k3<-dd(x[i]+h/2,y[i]+h*k2/2)
    k4<-dd(x[i]+h,y[i]+h*k3)
    y[i+1]<-y[i]+h*(k1/6+k2/3+k3/3+k4/6)
    points(x[i+1], y[i+1], type = "p",col='gold')
  }
}

analyticky(0,0.01,200)
euler1r(difder,0,1,0.1,7)
euler1rprumer(difder,0,1,0.1,7)
euler1r4interval(difder,0,1,0.1,7)
