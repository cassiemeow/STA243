########## Question 4 ##########
RMM = function(A, B, r, seed.num=243) {
  set.seed(seed.num)
  n = dim(A)[2]

  prob <- sapply(1:n, FUN = function(x) norm(A[,x],"2") * norm(B[x,],"2"))
  prob = prob / sum(prob)

  S = 0
  for (i in 1:r) {
    l = sample(1:n, size=1, replace=T, prob=prob)
    S = S + outer(A[,l], B[l,]) / prob[l]
  }
  S = S/r
  
  return(S)
}

A = as.matrix(read.csv(file="STA243_homework_1_matrix_A.csv", header=F))
B = as.matrix(read.csv(file="STA243_homework_1_matrix_B.csv", header=F))
C = A %*% B

result_RMM = list()
result_error = c()
r_vec = c(20,50,100,200)

for (i in 1:length(r_vec)) {
  r = r_vec[i]
  M = RMM(A, B, r)
  error = norm(M-C,"F") / norm(A,"F") / norm(B,"F")
  result_RMM[[i]] = M
  result_error = c(result_error,error)
}

print(data.frame(r=r_vec,error=result_error))

par(mfrow=c(2,2))
for (i in 1:length(r_vec)) {
  r = r_vec[i]
  title.text = paste0("r=",r)
  image(result_RMM[[i]], main=title.text) 
}
par(mfrow=c(1,1))

image(C)



########## Question 5 ##########
power_iteration = function(A, v0, eps = 1e-6, maxiter=100, num.trunc=NA) {
  # Please implement the function power_iteration that takes in the matrix X and initial vector v0 and returns the eigenvector.
  c_old = v0
  step = 0
  while (step < maxiter) {
    c = A %*% c_old
    c = as.vector(c / sqrt(sum(c^2)))
    if (!is.na(num.trunc)) {
      c[rank(-c) > num.trunc] = 0
      c = as.vector(c / sqrt(sum(c^2)))
    }
    crit = sqrt(sum(c-c_old)^2)
    if (crit < eps) {
      break
    } else {
      c_old = c
      step = step + 1
    }
  }
  if (step == maxiter) {
    print("The optimal value does not converge")
    print(paste0("Max iteration = ",maxiter))
    print(paste0("Crit = ",crit))
  }
  return(c)
}

set.seed(5)
E = matrix(rnorm(100), 10, 10)
v = c(1, rep(0, 9))
lams = 1:10
prods = c()
for (lambda in lams) {
  X = lambda*outer(v, v) + E
  v0 = rep(1, nrow(E))
  v0 = v0/sqrt(sum(v0^2))
  vv = power_iteration(X, v0)
  # vv = power_iteration(X, v0, num.trunc=3)
  # vv = power_iteration(X, v0, maxiter=10000) #
  prods = c(prods, abs(v %*% vv))
}
plot(lams, prods, "b")



########## Question 6 ##########
library(phangorn)
phi_generate = function(X, y, e=.1, seed.num=243*6) {
  n = dim(X)[1]
  d = dim(X)[2]
  r = round(d * log(n) / e)
  
  D = sample(c(1,-1), n, replace=T,prob = c(0.5, 0.5))
  DX = apply(X, 2, function(x) D*x)
  Dy = D * y
  
  HDX = apply(DX,2,fhm)
  HDy = fhm(Dy)
  
  index = sample(1:n,r,replace=T)
  X2 = apply(HDX, 2, function(x) sqrt(n/r)*sapply(index, function(p) as.vector(x)[p]))
  y2 = sqrt(n/r)*sapply(index, function(p) HDy[p])

  return(list(X=X2,y=y2))
}
                
set.seed(1)      
x <- matrix(runif(1048576*20), 1048576, 20)
y <- runif(1048576)
                        
time.full = system.time({b.full = solve(crossprod(x,x),crossprod(x,y))})[3]

e.list <- c(0.1,0.05,0.01,0.001)
for (i in length(e.list)) {
  e <- e.list[i]
  phi <- phi_generate(x,y,e)
  time.fast <- system.time({b.fast = solve(crossprod(phi$X,phi$X),crossprod(phi$X,phi$y))})
  time.fast[3]
}


##### Eunseong's new code #####
library(phangorn)
SHD_gen = function(X, y, e=.1, seed.num=243) {
  set.seed(seed.num)
  n = dim(X)[1]
  d = dim(X)[2]
  r = round(d * log(n) / e)
  
  D = sample(c(1,-1), n, replace=T)
  DX = apply(X,2,function(z) D*z)
  Dy = D * y

  HDX = apply(DX,2,fhm)
  HDy = fhm(Dy)

  S = sample(1:n, r, replace=T)
             
  S.fun = function(z) {
    z = as.vector(z)
    result = sqrt(n/r) * sapply(S, function(i) z[i])
    return(result)
  }
                                
  SHDX = apply(HDX,2,S.fun)
  SHDy = S.fun(HDy)

  return(list(X=SHDX, y=SHDy))
}

set.seed(1)
X = matrix(runif(1048576*20,0,1),1048576,20)
y = runif(1048576,0,1)

e_vec = c(.1, .05, .01, .001)
result = matrix(NA, length(e_vec), 3)
rownames(result) = e_vec
colnames(result) = c("Time", "Diff", "Relative Diff")

full.time = system.time({
  b1 = solve(crossprod(X,X),crossprod(X,y))  
})[3]

for (i in 1:length(e_vec)) {
  e = e_vec[i]
  SHD = SHD_gen(X,y,e)
  SHDX = SHD$X
  SHDy = SHD$y
  result[i,1] = system.time({
    b2 = solve(crossprod(SHDX,SHDX),crossprod(SHDX,SHDy))
  })[3]
  result[i,2] = norm(b1-b2,"2")
  result[i,3] = norm(b1-b2,"2")/norm(b1,"2")
}
print(full.time)
print(result)
