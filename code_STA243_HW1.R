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
sketched_OLS = function(X, y, e=.1, seed.num=243*6) {
  n = dim(X)[1]
  d = dim(X)[2]
  r = round(d * log(n) / e)
  S = NULL
  
  set.seed(seed.num)
  index = sample(1:n,r,replace=T)
  for (i in index) {
    S = cbind(S,as.matrix(sqrt(n/r)*replace(rep(0,n),i,1)))
  }
  
  temp = lapply(index, function(i) as.matrix(sqrt(n/r)*replace(rep(0,n),i,1)))
  
  for (i in 1:log(n,base=2)) {
    print(i)
    if (i==1) {
      H = matrix(c(1,1,1,-1),2,2)
    } else {
      H = rbind(cbind(H_old,H_old),cbind(H_old,-H_old))
    }
    H_old = H
  }
  H = 1/sqrt(n) * H
  
  D = diag(sample(c(1,-1),n,replace=T))
  
  X2 = t(S) %*% H %*% D %*% X
  y2 = t(S) %*% H %*% D %*% y
  
  b = solve(crossprod(X2), t(X2) %*%  y2)
  b = as.vector(b)
  
  return(b)
}