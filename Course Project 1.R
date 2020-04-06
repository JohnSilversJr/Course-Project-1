##########################################################################
##########################################################################
##### Problem 1 ##########################################################
##########################################################################
##########################################################################

##########################################################################
##### Supporting Functions: ##############################################
##########################################################################

##########################################################################
####  Norm - length of a vector       ####################################
##########################################################################
Norm = function(x){
  tot = 0
  for (i in 1:length(x)){
    tot=tot+x[i]^2
  }
  return(sqrt(tot))
}

##########################################################################
#### meshgrid  - used by vector field ####################################
##########################################################################
meshgrid = function(x,y){
  n = length(x)
  m = length(y)
  mat1 = matrix(0,nrow=m,ncol=n)
  mat2 = matrix(0,nrow=m,ncol=n)
  for (i in 1:m){
    mat1[i,]=x
  }
  for (i in 1:n){
    mat2[,i]=y
  }
  ans = list(X=mat1,Y=mat2)
  return (ans)
}

##########################################################################
#### Quiver  - used by vector field ######################################
##########################################################################

quiver <- function(x, y, u, v,
                   scale = 0.05, angle = 10, length = 0.1, ...) {
  stopifnot(is.numeric(x), is.numeric(y), is.numeric(u), is.numeric(v))
  
  arrows(x, y, x+scale*u, y+scale*v, angle=10, length=length, ...)
}

###########################################################################
##### Vector Field ####################
###########################################################################

VectorField = function(fun, xlim, ylim, n = 16,
                       scale = 0.05, col = "darkblue",xlab = "xlim", ylab="ylim",
                       main="",...) {
  stopifnot(is.numeric(xlim), length(xlim) == 2,
            is.numeric(ylim), length(ylim) == 2)
  
  xpts = seq(xlim[1],xlim[2],length.out=n)
  ypts = seq(ylim[1],ylim[2],length.out=n)
  
  M = meshgrid(xpts, ypts)
  
  x = M$X
  y = M$Y
  px=M$X
  py=M$Y
  for (i  in 1:n){
    for (j in 1:n){
      ans = fun(c(xpts[j],ypts[i]))
      px[i,j]=ans[1]
      py[i,j]=ans[2]
    }
  }
  
  
  plot(xlim, ylim, type="n",xlab=xlab,ylab=ylab,main=main); grid()
  quiver(x, y, px, py, scale = scale, col = col, ...)
  #return(list(px=px,py=py))
}

###########################################################################
##### Jacobian of a 2D vector function ####################################
###########################################################################
Jacobian2 = function(f,x0,h=1E-4){
  jax = matrix(0,nrow=2,ncol=2)
  xph = c(x0[1]+h,x0[2]);xmh=c(x0[1]-h,x0[2])
  yph = c(x0[1],x0[2]+h);ymh=c(x0[1],x0[2]-h)
  jax[,1]=(f(xph)-f(xmh))/(2*h)
  jax[,2]=(f(yph)-f(ymh))/(2*h)
  return(jax)
}

###########################################################################
##### Zeros - root finding for a vector function ##########################
###########################################################################
Zeros = function(f,x0,h=1E-4,tol=1E-4){
  i = 1
  p=c(1,1)
  while (Norm(p)>tol & i<100){
    p = solve(Jacobian2(f,x0),-f(x0)) # linear algebra step
    x0 = x0+p
    #print(i)
    #print(x0)
    i=i+1
  }
  return(x0)
}

###########################################################################
##### path - output is the points along a particular dynamic system########
###########################################################################
path = function(f,x0,deltat=0.01,N=1000,tol=1E-4){
  len = length(x0)
  points=matrix(0,ncol=len)
  points[1,] = x0
  n = 0
  p = c(1,1)
  while(Norm(p)>tol & n<N){
    n=n+1
    p = f(x0)*deltat
    x0=x0+p
    points = rbind(points,x0)
  }
  
  rownames(points)=0:n
  return(points)
}
###########################################################################
##### Find Root ###########################################################
###########################################################################
findRoot = function(f,xlim,ylim){
  eq = data.frame(x=1,y=1);idx=1
  x=seq(xlim[1],xlim[2],length.out=4)
  y=seq(ylim[1],ylim[2],length.out=4)
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      x0=c(x[i],y[j])
      eq[idx,]=round(Zeros(f,x0),4)
      idx=idx+1
    }
  }
  eq=unique(eq)
  return(eq)
}

##########################################################################
##########################################################################
##### Phase Portrait Function [Answer] ###################################
##########################################################################
##########################################################################

phasePortrait = function(f,xlim,ylim,scale){
  VectorField(f,xlim,ylim,scale=scale)
  eq=findRoot(f,xlim,ylim)
  print(eq)
  points(eq,bg="green",pch=21:25,col="green")
  XL = length(xlim)
  YL = length(ylim)
  x0 = c(-.75*XL,-.75*YL)
  N = 50
  for (i in 1:N){
    pts = path(f,x0,500,deltat=0.01)
    points(pts[,1],pts[,2],type="l")
  }
  x0 = c(-.75*XL,.5*YL)
  N = 50
  for (i in 1:N){
    pts = path(f,x0,500,deltat=0.01)
    points(pts[,1],pts[,2],type="l")
  }
  x0 = c(-.75*XL,.75*YL)
  N = 50
  for (i in 1:N){
    pts = path(f,x0,500,deltat=0.01)
    points(pts[,1],pts[,2],type="l")
  }
  x0 = c(.75*XL,-.75*YL)
  N = 50
  for (i in 1:N){
    pts = path(f,x0,500,deltat=0.01)
    points(pts[,1],pts[,2],type="l")
  }
  x0 = c(.75*XL,.5*YL)
  N = 50
  for (i in 1:N){
    pts = path(f,x0,500,deltat=0.01)
    points(pts[,1],pts[,2],type="l")
  }
  x0 = c(.75*XL,.75*YL)
  N = 50
  for (i in 1:N){
    pts = path(f,x0,500,deltat=0.01)
    points(pts[,1],pts[,2],type="l")
  }
}

f=function(x){c(-x[1]^3-4*x[1]-x[2],
                3*x[1])
}

phasePortrait(f,xlim=c(-3,3),ylim=c(-3,3),scale=0.01)

##########################################################################
##########################################################################
##### Problem 2 ##########################################################
##########################################################################
##########################################################################

###Individual recursion euqations###
#succeptible = x[1]-(alpha*x[1]*x[2])-((1-beta)*alpha*x[1]*x[2])
#infected = x[2]+(alpha*x[1]*x[2])-(lambda*x[2])-(beta*x[2])
#recovered = x[3]+(beta*x[2])+((1-beta)*alpha*x[1]*x[2])
#death = x[4]+(lambda*x[2])

### function to plot the dataframe as lines ###
PlotAll = function(df){
  n = ncol(df)
  colors = c("blue","red","orange","purple","green","yellow")
  plot(df[,1],df[,2],type = "l",col=colors[1],ylim=c(0,max(df)))
  for (i in 2:n){
    lines(df[,1],df[,i],col=colors[i])
  }
}

### recursion equations in a function ###
f = function(x){ c(((x[1]-(alpha*x[1]*x[2]))-((1-alpha)*alpha*x[1]*x[2])),
                   (x[2]+(alpha*x[1]*x[2])-(lambda*x[2])-(beta*x[2])),
                   (x[3]+(beta*x[2])+((1-alpha)*alpha*x[1]*x[2])),
                   (x[4]+(lambda*x[2])))
}
### rates ###
alpha = 0.00000001
beta = 0.8
lambda = 0.03

### start point ###
x0=c(7000000000,100,50,20)
print(f(x0))

### constraints to keep the function from pulling from a population that is not there ###
f2 = function(x){
  ans = 0
  #movement from succeptible to infected
  delta1 = (alpha*x[1]*x[2])+((1-alpha)*alpha*x[1]*x[2])
  if ((alpha*x[1]*x[2])+(1-alpha)*alpha*x[1]*x[2]<x[1]){delta1=(alpha*x[1]*x[2])+((1-beta)*alpha*x[1]*x[2])}else{delta1=(x[1])}
  #movement from succeptible to recovered
  delta2 = (alpha*x[1]*x[2])+((1-alpha)*alpha*x[1]*x[2])
  #movement from infected to recovered
  delta3 = (beta*x[2])+((1-alpha)*alpha*x[1]*x[2])
  #movement from infected to dead
  delta4 = (lambda*x[2])
  ans[1] = x[1]-delta1
  ans[2] = x[2]+delta1-delta3-delta4
  ans[3] = x[3]+delta2+delta3
  ans[4] = x[4]+delta4
  
  return(ans)
}

N=30
x=c(7000000000,100,50,20)
ans.x1=0
ans.x2=0
ans.x3=0
ans.x4=0
ans.sum=0
for (i in 1:N){
  ans.x1[i]=x[1]
  ans.x2[i]=x[2]
  ans.x3[i]=x[3]
  ans.x4[i]=x[4]
  ans.sum[i]=sum(x)
  x=f2(x)
}
result = data.frame(n=1:N,x1=ans.x1,x2=ans.x2,x3=ans.x3,x4=ans.x4,pop=ans.sum)
print(result)
tail(result)

PlotAll(result[1:5])