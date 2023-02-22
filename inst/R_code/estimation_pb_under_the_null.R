N=1000
P=20

data <- set_data_mococomo(betahat  = rnorm(N),
                          X = matrix(rnorm(P*N), ncol=P),
                          se= rep(1,N))
#fit mococomo model
maxiter   = 100
tol       = 1e-3
max_class = 10
mult      = 2
nullweight     =10

#working init
fit <- init.mococomo(data,
                     max_class = max_class,
                     mult   = mult,
                     upper     = upper

)

#works with two distribution (left and right)
tfit <- fit.mococomo  (data,nullweight = 10)


plot(tfit$post_assignment[,1], tfit$post_assignment[,2])

library(ashr)

tt <- ash( rnorm(N), rep(1,N), mixcompdist = "normal")
tt$fitted_g
tfit$f_list
