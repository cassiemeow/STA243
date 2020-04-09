########## Question 4 ##########
RMM = function(A, B, r, seed.num=243) {
  
  set.seed(seed.num)
  n = dim(A)[2]

  prob = rep(0, n)
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

par(mfrow=c(2,2))
for (i in 1:length(r_vec)) {
  r = r_vec[i]
  title.text = paste0("r=",r)
  image(result_RMM[[i]], main=title.text) 
}
par(mfrow=c(1,1))

image(C)
