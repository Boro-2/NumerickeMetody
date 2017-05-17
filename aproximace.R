#1) po castech linearni
# 2.3.4) POLYNOM vandermont, lagrangeuv tvar, newtonuw dodelat newtona
myPlot<-function(x,y,consX,consY){
  rozsiritX<-consX
  rozsiritY<-consY
  intervalX<-((x[length(x)])-x[1])*(rozsiritX) #rozsiri graf horizontalne
  intervalY<-(max(y)-min(y))*rozsiritY #rozsiri graf vertikalne 
  plot(c(x[1]-intervalX,x,x[length(x)]+intervalX),c(min(y)-intervalY,y,max(y)+intervalY),type="b")
}
funcVector<-function(iks,po){#c je polynom stupne
  fcVct<-matrix(nrow = po+1,ncol = length(iks))
  n=length(iks)
  fcVct[1,]=rep(1,n);
  for(i in 1:po){
    fcVct[i+1,]<-c(fcVct[i,]*iks)
  }
  return(fcVct)}

x<-c(1,2,3,4,5)

y<-c(0,-1,-3,-2,1)

# 
# x<-c(1,2,3,4)
# y<-c(3,7,12,25)

# y2<-c(10,5,-5,30,-6,-15,2,-22)
# y3<-c(6,1,-14,30,2,16,13,-22)
# yc<-c(-6,1,4,7,9,11,12,17,10,5,-5,30,-6,-15,2,-22,6,1,-14,30,2,16,13,-22)
# ym<-matrix(c(y,y2,y3),ncol=length(y))

xs<--pi
delta<-pi/(1.9)
xsin<-vector(mode = "numeric")
ysin<-c()
ind<-0
while(xs<(3*pi)){
  ys<-sin(xs)
  #print(ys)
    xsin[ind]<-round(xs,6)
    ysin[ind]<-round(ys,6)
    xs<-xs+delta
    ind<-ind+1
}
#print(xsin)
#print(ysin)


# x<-c(1,3,4)
# y<-c(2,5,6)


#myPlot(xsin,ysin,0.3,6)#vektory x, y a konstanty rozsireni grafu pro osu x, y

plot(x,y, type = "b",xlim = c(min(x)-3,max(x)+3),ylim = c(min(y)-30,max(y)+30))
# abline(v=xsin[1], col='maroon', lwd=1.5)
# abline(v=xsin[length(xsin)], col='maroon', lwd=1.5)

linearPol <- function(x,y,inp){
  pos<-0
  outp<-0
  if(inp<x[1]){ pos<-1}
  else if(inp>x[length(x)]){ pos<-length(x)-1 }
  else{
  for(i in 1:length(x)){
    if(x[i]>inp){ pos<-(i-1)
    break}
  }}
  outp<-(y[pos]*(x[pos+1]-inp)/(x[pos+1]-x[pos]))+(y[pos+1]*(inp-x[pos])/(x[pos+1]-x[pos]))
  
  points(inp,outp,type='p',col='aquamarine4')
  return(outp)
} # works fine

vandermondeM<- function(x,y,x0){
  b<-y
  matic<- matrix(nrow=length(x),ncol=length(x),byrow=TRUE)
  for(i in 1:length(x)){
    for(j in 1:length(x)){
      matic[i,j]<-x[i]^(j-1)
    } }
  res<-solve(matic,y)
  y0<-0
  for(i in length(res):1){
    y0<-y0*x0+res[i]
  }
  points(x0,y0, type='p',col='sienna')
} #funguje

lagrange<-function(x,y,x0){
  L<-0
  for(i in 1:length(x)){
    l<-1
    for(j in 1:length(x)){
      if(j!=i){
        l<-l*(x0-x[j])/(x[i]-x[j])
      }
    }
    L<-L+y[i]*l
  }
  points(x0,L,type='p',col='steelblue2')
  return(L)
  
  
} # funguje 

kubSpl0<-function(x,y,x0){ #nekde je chyba ve vkladani derivaci do matice, reseni avraci na posled 2 mistech 0-ly
  s<-(length(x)-1) #pocet intervalu
  mat<-matrix(nrow=4*s,ncol = 4*s)
  rightS<-c()
  
  for(i in 1:nrow(mat)){
    for(j in 1:nrow(mat))
      mat[i,j]<-0
  }
  
  ind<-1
  for(a in seq(1,2*s,2)){
    for(b in 1:4){
      mat[a,(ind-1)*4+b]<-x[ind]^(b-1)
      mat[a+1,(ind-1)*4+b]<-x[ind+1]^(b-1)
      
    }
    ind<-ind+1
  }
  
  for(a in 1:(s-1)){#navaznost 1 derivaci
    for(b in 2:4){
      mat[a+2*s,(a-1)*4+b]<-(b-1)*x[a+1]^(b-2)
      mat[a+2*s, a*4+b]<-(-1)*(b-1)*x[a+1]^(b-2) #bylo tu x[a+1]
      
    }
  }
  ind<-2 
    for(a in (3*s):(4*s-2)){#navaznost druhych derivaci
      
        mat[a,(ind-2)*4+3]<-2
        mat[a, (ind-1)*4+3]<-(-2) # minus ano ci ne
        
        mat[a,(ind-2)*4+4]<-6*x[ind]
        mat[a, (ind-1)*4+4]<-(-6)*x[ind+1] # minus ano ci ne? bylo x[ind+1]
        ind<-ind+1
    }
  
  mat[4*s-1,3]<-2
  mat[4*s-1,4]<-6*x[1]
  mat[4*s,(s-1)*4+3]<-2
  mat[4*s,(s-1)*4+4]<-6*x[length(x)]

  
  
  ind<-2
  
  for(i in 2:(4*s)){
    if(i==(2*s)) {next}
    if(i<(2*s)){#bylo jen 2s
      rightS[i]<-y[ind]
      ind<-ind+0.5 #po dvou cyklech zvedne ind o 1
    }
    else {rightS[i]<-0
    }
  }
  rightS[1]<-y[1]
  rightS[2*s]<-y[length(y)]

  xx<-solve(mat,rightS)
  
  print(xx)
  i<-(-10)
  if(x0<x[1]) {i<-0} #na ktere casti kub splinu budu fcni hodnotu pro x0 hledat
  else if(x0>x[length(x)]) {i<-length(x)-2}
  else {
    for(k in 2:length(x)){ 
      if(x[k]>x0) {
        i<-k-2
        break}}} 
 
  y0<-xx[i*4+1]+xx[i*4+2]*x0+xx[i*4+3]*x0*x0+xx[i*4+4]*x0*x0*x0
 
  points(x0, y0, type = "p",col='darkgreen')
  print(y0)
  return(y0)
}


kubSpl<-function(x,y,x0){
  ma<-matrix(nrow=(length(x)-2),ncol=(length(x)-2))
  g<-c()
  
  ma[1,1]= (x[3]-x[1])/3
  g[1]= (y[3]-y[2])/(x[3]-x[2])-(y[2]-y[1])/(x[2]-x[1])
  
  for(o in 2:nrow(ma)){#vytvorim 3diagonalni matici
    dx1<-(x[o+1]-x[o])
    dx2<-(x[o+2]-x[o+1])
    ma[o,o]<-(dx1+dx2)/3
    ma[o-1,o]<-dx1/6
    ma[o,o-1]<-ma[o-1,o]
    g[o]<-(y[o+2]-y[o+1])/(x[o+2]-x[o+1])-(y[o+1]-y[o])/(x[o+1]-x[o])
  }
  
  mi<-c(0,tridiagonal(ma,g),0) #vektor s resenim 3diagonal, prvni a posledni druha derivace =0
  ai<-c()
  bi<-c()
  for(u in 1:(length(x)-1)){
    bi[u]<-y[u]-mi[u]*(x[u+1]-x[u])*(x[u+1]-x[u])/6 
    ai[u]<-(y[u+1]-y[u])/(x[u+1]-x[u])+mi[u]*(x[u+1]-x[u])/6-mi[u+1]*(x[u+1]-x[u])/6
  }
  
  
  i<-1
  y0<-0
  
  if(x0<x[1]) {i<-1} #na ktere casti kub splinu budu fcni hodnotu pro x0 hledat
  else if(x0>x[length(x)]) {i<-length(x)-1}
  else {
    for(k in 1:length(x)){ 
      if(x[k]>x0) {
        i<-(k-1)
        break}      }} 
  
  y0<-(x0-x[i])^3/(6*(x[i+1]-x[i]))*mi[i+1]+(x[i+1]-x0)^3/(6*(x[i+1]-x[i]))*mi[i]+ai[i]*(x0-x[i])+bi[i]
  #print(y0)
  points(x0, y0, type = "p",col="red")
  return(y0)
  
  
} #funguje 

kubSplVectorx0<-function(x,y,x0){ #vstup x0 je vektor
  ma<-matrix(nrow=(length(x)-2),ncol=(length(x)-2))
  g<-c()
  
  ma[1,1]= (x[3]-x[1])/3
  g[1]= (y[3]-y[2])/(x[3]-x[2])-(y[2]-y[1])/(x[2]-x[1])
  
  for(o in 2:nrow(ma)){
    dx1<-(x[o+1]-x[o])
    dx2<-(x[o+2]-x[o+1])
    ma[o,o]<-(dx1+dx2)/3
    ma[o-1,o]<-dx1/6
    ma[o,o-1]<-ma[o-1,o]
    g[o]<-(y[o+2]-y[o+1])/(x[o+2]-x[o+1])-(y[o+1]-y[o])/(x[o+1]-x[o])
  }
  
  mi<-c(0,tridiagonal(ma,g),0) #vektor s resenim 3diagonal, prvni a posledni druha derivace =0
  ai<-c()
  bi<-c()
  for(u in 1:(length(x)-1)){
    bi[u]<-y[u]-mi[u]*(x[u+1]-x[u])*(x[u+1]-x[u])/6 
    ai[u]<-(y[u+1]-y[u])/(x[u+1]-x[u])+mi[u]*(x[u+1]-x[u])/6-mi[u+1]*(x[u+1]-x[u])/6
  }
  #print(ai)
  #print(bi)
  y0<-c()
  for(r in 1:length(x0)){
  i<-1
  y0[r]<-0
  
  if(x0[r]<x[1]) {i<-1} #na ktere casti kub splinu budu fcni hodnotu pro x0 hledat
  else if(x0[r]>x[length(x)]) {i<-length(x)-1}
  else {
    for(k in 1:length(x)){ 
      if(x[k]>x0[r]) {
        i<-(k-1)
        break}}} 
  
  y0[r]<-(x0[r]-x[i])^3/(6*(x[i+1]-x[i]))*mi[i+1]+(x[i+1]-x0[r])^3/(6*(x[i+1]-x[i]))*mi[i]+ai[i]*(x0[r]-x[i])+bi[i]
  #if(draw==1)
  points(x0[r], y0[r],type='p',col='red',bg='red',pch=21,cex=0.3)
  }
  
  return(y0)
  
} #funguje 
mNejCtv<-function(x,y,x0,pol){##jednoducha verze, uz funkcni
  mat<-matrix(ncol = (pol+1),nrow = (pol+1))
  s<-(pol+1)
  c<-length(x)
  vysl<-c()
  for(i in 1:s){#radky
    su<-0
    for(j in i:s){#sloupce
      sum<-0
      for(k in 1:c){#vsechny x
        sum<-sum+x[k]^(i+j-2)
        if(j==i)#prava strana 
        {
          su<-su+y[k]*x[k]^(i-1)
        }
      }
      mat[i,j]<-sum
      mat[j,i]<-sum
    }
    vysl[i]<-su #prava strana
  }
 # print(mat)
 # print(vysl)
  xx<-solve(mat,vysl)
  print(xx)
  y0<- xx %*% funcVector(x0,pol)
  points(x0, y0, type = "l",col='red')
}

mNejCtvMultiY<-function(x,y,x0,pol,testing){#y je matice kde radek je mnozina y, NA 95% NEFUNGUJE, NUTNO OPRAVIT NA ZAKLADE mNejCtv
  #vytvari graf pro kazdy set y, ale pomoci lUDekomozice resi vse rychle
  mat<-matrix(ncol = (pol+1),nrow = (pol+1))
  s<-(pol+1)
  c<-length(x)
  vysl<-matrix(ncol = s,nrow = nrow(y))
  for(v in 1:nrow(y)){
    for(i in 1:s){
      su<-0
      for(j in i:s){
        sum<-0
        for(k in 1:c){
          sum<-sum+x[k]^(i+j-2)
        }
        mat[i,j]<-sum
        mat[j,i]<-sum
        
        su<-su+y[v,j]*x[j]^(i-1)
      }
      vysl[v,i]<-su
    }
  }
  print(vysl)
  xx<-lUD(mat,y)
  if(testing==TRUE){
    for(j in 1:nrow(xx))
    {
      myPlot(x,y[j,],0.5,2)
      for(p in 1:1000){
        y0<- xx[j,] %*% funcVector(x0,pol)
        points(x0, y0, type = "p",col='darkblue')
        x0<-runif(1,min(x)-2,max(x)+2)
      }
    }
  }
  #y0<-dejHod(xx[1,],x0)
  # print(y0)
  #points(x0, y0, type = "p",col='darkblue')
}


newtonInterpol<-function(x,y,x0){
  mat<-matrix(nrow = length(x),ncol = length(x)+1)
  mat[,1]<-x
  mat[,2]<-y
  k<-nrow(mat)
  for(j in 3:ncol(mat)){
    k<-k-1
    for(i in 1:k){
      mat[i,j]<-(mat[i+1,j-1]-mat[i,j-1])/(mat[i+j-2,1]-mat[i,1])
    }
  }
  #print(mat[1,2:ncol(mat)])
  y0<-mat[1,2]
  sumx<-1
  for(h in 3:ncol(mat))
  {
    sumx<-sumx*(x0-x[h-2])
    y0<-y0+mat[1,h]*sumx
  }
  points(x0, y0, type = "p",col='red3')
}

newtonInterpolFast<-function(x,y,x0,a){
  #http://mi21.vsb.cz/sites/mi21.vsb.cz/files/unit/numericke_metody_interaktivne.pdf (algoritmus 20)
  pocetX<-length(x)
  a[1]<-y[1]
  if(length(a)!=pocetX)
  for(i in (length(a)+1):pocetX){#resim pouze indexy a, ktere nejsou vypocitane
    S<-0
    P<-1
    for(j in 1:(i-1)){
      S<-S+a[j]*P
      P<-P*(x[i]-x[j])
    }
    a[i]<-(y[i]-S)/P
  }
  
  S<-0
  P<-1
  for(k in 1:(pocetX)){
    S<-S+a[k]*P
    P<-P*(x0-x[k])#tady sice pouziju posledni x, ale k vypoctu neslouzi
  }
  points(x0, S, type = "p",col='gold')
  return(a)
}

a<-c()

# print(system.time(for(p in 1:1000)a<- newtonInterpolFast(xsin,ysin,runif(1,-6,2*pi+10),a)))
#print(system.time(for(p in 1:1000)kubSpl(xsin,ysin,runif(1,-6,2*pi+10))))
#print(system.time(for(p in 1:5000)a<- mNejCtv(xsin,ysin,runif(1,-6,2*pi+10),8)))# minn<-min(x)-2

minn<-min(x)-2
 maxx<-max(x)+2
 #print(system.time(for(p in 1:5000)a<- newtonInterpolFast(x,y,runif(1,minn,maxx),a)))
# print(system.time(for(p in 1:5000) newtonInterpol(x,y,runif(1,minn,maxx))))
# print(system.time(for(p in 1:5000) vandermondeM(x,y,runif(1,minn,maxx)) ))
# print(system.time(for(p in 1:5000) lagrange(x,y,runif(1,minn,maxx)) ))

#mNejCtvMultiY(x,ym,runif(1,minn,maxx),5,TRUE) #nedavat do cyklu!!! ma v sobe vlastni
xss<-seq(0,100,2)
 #yss<-kubSplVectorx0(x,y,xss)
 #plot(xss,yss)
for(p in 1:1){
 #a<- newtonInterpolFast(x,y,2,a)
#newtonInterpol(x,y,runif(1,minn,maxx))
mNejCtv(x,y,sort(runif(100,-2,8)),2)
#kubSpl(x,y,runif(1,minn,maxx))
#linearPol(x,y,runif(1,minn,maxx))
#vandermondeM(x,y,runif(1,minn,maxx))
#lagrange(x,y,runif(1,minn,maxx))
 }

