fce<-function(x){return(x*x+1)} # Integral from 1 to 10 is 342
fce2<-function(x){return(2*x+5)} # Integral from 0 to 5 is 50
fce3<-function(x){return(x*x*x-3*x*x-3*x+20)} # Integral from -2 to 5 is 127,75
fce4<-function(x){return(5^(2*x))} # Integral from -2 to 5 is 3 033 862
fce5<-function(x){return(abs(2*cos(x^2)*x^2)+sin(2*x)+exp(x)/10+10)} ## romberg to pocita presne

a<--2
b<-5
x<-c(sort(runif(400,a-1,b+1)))
y<-c()
for(i in 1:length(x)){y[i]<-fce5(x[i])}
min(y)

plot(x,y,xlim = c(min(x)-5,max(x)+5),ylim = c(min(y)-10,max(y)+10))
abline(v=a,lwd=2,col="gray")
abline(v=b,lwd=2,col="gray")


#NASTUDOVAT METODA ZLATEHO REZU PRO MIN A MAX FCE

#MC-Monte Carlo

intgAnal<-function(pol,a,b){
  
  intg<-c()
  down<-1
  for(i in (length(pol)):1){
    intg[i]<-pol[i]/down
    down<-down+1
  }
  intg[length(pol)+1]<-0
  
  suma<-0
  sumb<-0
  for(j in 1:length(intg)){
    suma<-(suma)*a+intg[j]
    sumb<-(sumb)*b+ intg[j]
  }
  
  print(sumb-suma)
  return(sumb-suma)
  
}


MCGeo<-function(x,y,a,b,co,n){#vstupni sour x,y; interval a,b, ; uplimimit max(f(x))*co >f(x) should be 1< ; n-pocet nahodnych nastrelu
  xRan<-sort(c(runif(n,a,b)))
  yFun<-kubSplVectorx0(x,y,xRan)
  uplim<-co*max(yFun)
  
  yRan<-c(runif(n,0,(uplim-1)))
  S<-0
  
  for(i in 1:length(xRan)){if(yRan[i]<yFun[i]){S<-S+1}}
  P<-S/n #percentage of points below curve
  print((b-a)*uplim*P)
  return((b-a)*uplim*P)
}

MCMeanY<-function(x,y,a,b,n){#points x,y; interval a,b, ; n number of randomly choosen points from <a,b>
  return(mean(kubSplVectorx0(x,y,c(runif(n,a,b))))*(b-a))
}

simpsonRule<-function(x,y,a,b){
  val<-kubSplVectorx0(x,y,c(a,(a+b)/2,b))
  print(val)
  res<-(b-a)/6*(val[1]+4*val[2]+val[3])
  return(res)
}

f<-function(x){return(sin(x))}

SimpsonRulefce<-function(f,a,b){ #melo bz tam vstoupit n, jako pocet parametru
  return((f(a)+4*f((a+b)/2)+f(b))*(b-a)/6)
}

GaussIntegral3Orig<-function(f,a,b){
  h<-abs(b-a)/2
  fmiddle<-f((a+b)/2)
  fleft<-f((a+b)/2-h*sqrt(3/5))
  fright<-f((a+b)/2+h*sqrt(3/5))
  integral<- h/9*(8*fmiddle + 5*fleft +5*fright)
  return(integral)
}
GaussIntegral3<-function(x,y,a,b){
  h<-abs(b-a)/2
  y123<-kubSplVectorx0(x,y,c((a+b)/2-h*sqrt(3/5),(a+b)/2,(a+b)/2+h*sqrt(3/5)))
  integral<- h/9*(8*y123[2] + 5*y123[1] +5*y123[3])
  #print(integral)
}

RombergovaKvadratura<-function(f,a,b,iter){
  mat<-matrix(nrow = iter,ncol = iter)
  initH=abs(b-a)
  h<-initH/100
  parts<-100
  for(i in 1:iter){
    sum<-0
    begin<-a
    
    for(k in 1:parts){
      sum<-sum+GaussIntegral3Orig(f,begin,(begin+h))
      begin<-begin+h
    }
    mat[i,1]<-sum#naplneni 1. sloupce obsahy s ruznymi kroky
    parts<-parts*2
    h<-h/2
    print(i)
  }
  
  for(j in 2:iter){#pro kazdy dalsi sloupec (1. mam, tazke 2 az pocet iteraci)  
    for(i in 1:(iter+1-j)){
      mat[i,j]<-(4^(j-1)*mat[i+1,j-1]-mat[i,j-1])/(4^(j-1)-1) #pomoci dvou predchozich obsahu spocitam dalsi, s tim ze obsah spocitany pomoci
    } # mensiho kroku ma mnohem vetsi vahu 
    print(j)
  }
  
  j<-c(seq(1,iter,1))
  print(mat-141.3806)
  #print(mat[1,])
  plot(j,print(mat[,1]),type = "b")
   points(j,print(mat[1,]),type = "b",col='red')
  
}
#simpsonRule()
#RombergovaKvadratura(fce5,a,b,12)
#GaussIntegral3Orig(fce5,a,b)
#GaussIntegral3(x,y,a,b)

#print(SimpsonRuleSkvor(f,0,pi))
# vyslGeo<-c()
# vyslMean<-c()
# for(p in 1:4){
  #print((vyslGeo[p]<-MCGeo(x,y,a,b,1.1,800)))
  #print((vyslMean[p]<-MCMeanY(x,y,a,b,8000)))
# }
# meanG<-mean(vyslGeo)
# sdG<-sd(vyslGeo)
# meanM<-mean(vyslMean)
# sdM<-var(vyslMean)
# print(meanG)
# print(sdG)
# print(meanM)
# print(sdM)



#   for(i in 1:10000)
# {
#   xx<-runif(1,0,100)
#   points(xx,f(xx),type='p',col='gold',bg='gold',pch=21,cex=0.5)
#   
# }


