mat1 <-function(r,c){#funkce vytvori matici {rxc} a naplni nahodnymi cisli (v tomhle pripade 1-7)
  matt<-matrix(rep(round(runif(r*c,-7,7),0)),r,byrow=T)
  return(matt)
}
op<-4
po<-4
m1<-matrix(c(1,2,0,0,2,5,4,0,0,4,3,1,0,0,3,2),nrow=4,byrow = TRUE)

b1<-c(2,1,4,7)

m2<-mat1(4,4)
b2<-c(5,4,1,2)

m3<-matrix(c(0,4,4,0,5,1,2,7,5),byrow=TRUE,nrow=3)
b3<-c(3,9,-6)

m4<-matrix(c(2,5,1,3),2)
b4<-c(-1,-1)

ml<-matrix(c(2,4,8,3,2,3,5,3,-1,-1,-3,-2,1,2,4,2),nrow=4)
bl<-c(4,6,12,6)

blx3<-matrix(c(4,6,12,6 ,8,-6,1,17,-5,3,4,6),ncol=ncol(ml),byrow=TRUE)

ThomasMatice<-matrix(c(-1,2,-1,0,1,1,1,-2,0),nrow = 3,byrow=TRUE) # byrow zn. vytvori se po radcich, kdyby tam nebzlo byrow = true, tak by si to musel zapsat vektorove to yn jeden sloupec druhy a treti
vektorPravychStran<-matrix(c(0,2,1), ncol=ncol(ThomasMatice),byrow=TRUE) #vektor pravych stran

#pripojeni reseni k matici
mplusright<-cbind(m3,b3)

msim<-matrix(c(2,2,2,2,8/9,2/9,2,32/81,2/81),ncol = 3,byrow = TRUE)
rsim<-c(2,2/3,2/5)

switchrows<-function(m,index){
  swap<-FALSE
  for(i in 1:nrow(m)){
    if(i==index)next
    if(m[i,index]!=0){
      temp<-m[index,]
      m[index,]<-m[i,]#prohodim radky
      m[i,]<-temp
      swap<-TRUE
    }
  }
  if(swap==FALSE){print("not able to swap rows")}
  return(m)
}
gaEl<-function(m){
  r<-nrow(m)
  s<-ncol(m)
  for(k in 1:(r-1)){
    #print(m)
    if(m[k,k]==0){m<-switchrows(m,k)}
    #print(m)
    for(i in (k+1):r){
      if(m[i,k]==0){next}
      for(j in (k+1):s){
        m[i,j]<-m[i,j]-m[i,k]*m[k,j]/m[k,k]
      }
    }
  }
  res<-c()#vektor pro zapis reseni pro x1:x(n)
  if(nrow(m)+1==ncol(m))#jen kontroluju jestli je to matice { n X (n+1)} posledni sloupec jsou prave strany 
  {
    sl<-ncol(m)#index posledniho sloupce/prave strany
    
    for(k in nrow(m):1){
      sm<-0
      
      if((k)<nrow(m))
      for(i in (k+1):nrow(m)){
        sm<-sm+res[i]*m[k,i]
      }
      
      res[k]<-(m[k,sl]-sm)/m[k,k]
    }
  }
  print(res)
  return(res)
  
}

lUD<-function(m, b){ # m i b jsou matice, 
  #kazdy radek matice b je prava strana linearni soustavy rovnic, reseni jednotlivych pravych stran budou take v radku
  
  siz<-ncol(m)
  mU<-matrix(nrow=siz,ncol=siz)
  mL<-matrix(nrow=siz,ncol=siz)
  for(r in 1:siz){#radek pro Upper, sloupc pro Lower
    
        mU[r,r]<-m[r,r] #najde prvkz na diagonale
        mL[r,r]<-1
        if(r>1){
          for(h in 1:(r-1)){
            mU[r,r]<-mU[r,r]-(mL[r,h]*mU[h,r])    }
        }
    
    
    for(i in r+1:siz){#sloupec pro Upper, radek pro Lower
      if(i==siz+1) break
      u<-0
      if(r>1){
        for(j in 1:(r-1)){
        u<-u+mL[r,j]*mU[j,i]  }
      }
      mU[r,i]<-m[r,i]-u
      mL[r,i]<-0
      
      u<-0
      if(r>1){
        for(j in 1:(r-1)){
        u<-u+mL[i,j]*mU[j,r]  }}
      mL[i,r]<-(m[i,r]-u)/mU[r,r]
      mU[i,r]<-0
    }
  }
  print(mU)
  print(mL)
  y<-matrix(nrow=nrow(b),ncol=ncol(b))
  x<-matrix(nrow=nrow(b),ncol=ncol(m))
  
  for(r in 1:nrow(b)){#pro kazdy radek v b matici
    for(i in 1:nrow(mL)){ # pro kazdou hodnotu b1-bn
      y[r,i]<-b[r,i]
      if(i>1){
        for(j in 1:(i-1)){
          y[r,i]<-y[r,i]-mL[i,j]*y[r,j]
        }
      }
    }
  }
  for(r in 1:nrow(b)){
    for(k in nrow(mU):1){
      x[r,k]<-y[r,k]
      if(k<nrow(mU)){
        for(j in (k+1):ncol(mU)){
          x[r,k]<-x[r,k]-mU[k,j]*x[r,j]
        }
      }
      x[r,k]<-x[r,k]/mU[k,k]
    }
  }
  #print(x)
  #print(xxx)
  return(x)
}



jacobi<-function(m,b){#funguje jen s nekterymi maticemi 
  siz<-nrow(m)
  citac<-0
  x0<-vector(mode = "numeric",siz)
  x1<-vector(mode = "numeric",siz)
  x1<-c(rep(0,siz))
  
  repeat{
    x0<-x1
    for(i in 1:siz){
      u<-0
      
      for(j in 1:siz) { 
        if(j==i){next}
        u<-u+m[i,j]*x0[j] }
      
      x1[i]<-(b[i]-u)/m[i,i]
      
    }
    cat(citac, "... ",x1,"\n")
    citac<-citac+1
    #if(citac>100){break}
    if(all(x0==x1)){break}
  }
}

gaussSeidl<-function(m,b){#funguje jen s nekterymi maticemi 
  siz<-nrow(m)
  citac<-0
  x0<-vector(mode = "numeric",siz)
  x1<-vector(mode = "numeric",siz)
  x1<-c(rep(0,siz))
  repeat{
    x0<-x1
    for(i in 1:siz){
      u<-0
      
      for(j in 1:siz) { 
        if(j==i)next
        if(j>i){u<-u+m[i,j]*x0[j]}
        else {u<-u+m[i,j]*x1[j]} 
      }
      
      x1[i]<-(b[i]-u)/m[i,i]
      
    }
    cat(citac, "... ",x1,"\n")
    citac<-citac+1
    if(all(x0==x1)){break}
  }
}
 
tridiagonal<-function(m,rs){#m is matrix, rs right side vector
  siz<-nrow(m)
  a<-c()
  b<-c()
  c<-c()
  
  for(i in 1:siz){
    a[i]<-m[i,i]
    if(i<siz){
      b[i]<-m[i+1,i]
      c[i]<-m[i,i+1]
    }
  }
  d<-rs
  
  u<-vector(mode="numeric",length=siz)
  y<-vector(mode="numeric",length=siz)
  x<-vector(mode="numeric",length=siz)
  u[1]<-a[1]
  y[1]<-d[1]
  
  for(i in 2:siz){
    k<-b[i-1]/u[i-1]#ve vicherovi je b[i] ale maji posunuty index
    u[i]<-a[i]-k*c[i-1]
    y[i]<-d[i]-k*y[i-1]
  }
  
  x[siz]<-y[siz]/u[siz]
  
  for(i in seq((siz-1),1,-1)){
    x[i]<-(y[i]-c[i]*x[i+1])/u[i]
      }
  return(x)
}


#gaEl(mplusright)
#gaEl(ml)
#gaussSeidl(ml,bl)
#jacobi(m4,b4)
aaa<-lUD(ml,blx3)
aaa<-lUD(ThomasMatice,vektorPravychStran)
#tridiagonal(m1, b)

print(aaa)