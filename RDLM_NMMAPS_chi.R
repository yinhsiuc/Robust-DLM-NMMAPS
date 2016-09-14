

#####################################################################################
setwd("C:/Users/yinhsiuc/Desktop/BIOSTAT-Yin-Hsiu Chen/Bhramar/Dissertation Project 1/Bayesian DLM/NMMAP/")
load("NMMAP-CHI.RData")
require(dlnm)
require(MASS)
require(MCMCpack)
require(splines)
require(rootSolve)
require(CompQuadForm)

## internal function
omega.fn <- function(eta1, eta2) {
	V <- diag(sqrt(exp(eta1*(0:(nlag-1)))))
	M <- diag(exp(eta2*(0:(nlag-1))))
	W <- cov2cor(tcrossprod(M) + tcrossprod((diag(nlag) - M)%*%rep(1, nlag)))
	return(V%*%W%*%V)
}

## data cleaning
nmonth <- 5114
# nlag <- 14
death <- chi.dat[1:nmonth, "death"]+chi.dat[(nmonth+1):(2*nmonth), "death"]+chi.dat[(2*nmonth+1):(3*nmonth), "death"]
time <- seq(-nmonth/2 + 0.5, length = nmonth)
addlag <- data.frame(matrix(NA, nmonth, 7))
colnames(addlag) <- paste("l", 8:14, "pm10tmean", sep = "")
addlag[8:nmonth, 1:7] <- chi.dat[1:(nmonth-7), seq(208, 244, 6)]
dat <- cbind(death, time, chi.dat[1:nmonth, c(3, 13, 17, 198:199, 32, seq(208, 244, 6))], addlag)
dat$dow <- as.factor(dat$dow)
dat1 <- dat[complete.cases(dat), ]

nmethod <- 11
nmonth <- 5114
nlag <- 15
time <- rep(seq(-nmonth/2 + 0.5, length = nmonth), 3)
addlag <- data.frame(matrix(NA, nmonth, 7))
colnames(addlag) <- paste("l", 8:14, "pm10tmean", sep = "")
addlag[8:nmonth, 1:7] <- chi.dat[1:(nmonth-7), seq(208, 244, 6)]
addlags <- rbind(addlag, addlag, addlag)
age2ind <- 1*(chi.dat$agecat == 2)
age3ind <- 1*(chi.dat$agecat == 3)
dat <- cbind(time, chi.dat[, c(3:4, 8, 13, 17, 198:199, 32, seq(208, 244, 6))], addlag, age2ind, age3ind)
dat$dow <- as.factor(dat$dow)
dat$agecat <- as.factor(dat$agecat)
dat1 <- dat[complete.cases(dat), ]
n <- dim(dat1)[1]

coef.est <- coef.cum <- upper <- lower <- upper.cum <- lower.cum <- cover <- cover.cum <- matrix(0, nlag, nmethod)
var.est <- array(0, dim = c(nlag, nlag, nmethod))
var.cum <- matrix(0, nlag, nmethod)

df <- 4
C <- onebasis(0:(nlag-1), fun = "bs", degree = df, intercept = TRUE)
R <- matrix(0, nlag - df - 1, nlag)
for(i in 1:(nlag - df - 1)) {
	R[i, i:(i+df+1)] <- c(1, -5, 10, -10, 5, -1)
}
RTR <- crossprod(R)

tran1 <- matrix(0, 2*nlag, nlag + df + 1)
tran1[1:nlag, 1:nlag] <- diag(nlag)
tran1[(nlag + 1):(2*nlag), (nlag + 1):(nlag + df + 1)] <- C
time_all <- rep(0, 8)

## GLM
# time_tag <- Sys.time()
formu1 <- formula(paste("death ~ dow + agecat + ns(tmpd, 6) + ns(dptp, df = 3) + ns(rmtmpd, df = 6) + ns(rmdptp, df = 3) + ns(time, df = 98) + I(ns(time, df = 14)*age2ind) + I(ns(time, df = 14)*age3ind) + pm10tmean +", 
	paste("l", 1:14, "pm10tmean", sep = "",  collapse = "+")))

mod.GLM <- glm(formu1, data = dat1, family = poisson(link = "log"), control = glm.control(epsilon = 1e-10, maxit = 1000),
               na.action = na.omit)
# Sys.time() - time_tag

coef.est[, 1] <- coef(mod.GLM)[154:168]
coef.cum[, 1] <- cumsum(coef.est[, 1])
var.est[, , 1] <- vcov(mod.GLM)[154:168, 154:168]

## DLM
Z <- as.matrix(dat1[, 9:23])%*%C
dat2 <- cbind(dat1, Z)
formu2 <- formula(paste("death ~ dow + agecat + ns(tmpd, 6) + ns(dptp, df = 3) + ns(rmtmpd, df = 6) + ns(rmdptp, df = 3) + ns(time, df = 98) + I(ns(time, df = 14)*age2ind) + I(ns(time, df = 14)*age3ind) + ", 
	paste("b", 1:(df+1), sep = "",  collapse = "+")))
mod.DLM <- glm(formu2, data = dat2, family = poisson(link = "log"), control = glm.control(epsilon = 1e-10, maxit = 1000),
               na.action = na.omit)
coef.est[, 2] <- C%*%coef(mod.DLM)[154:(154 + df)]
coef.cum[, 2] <- cumsum(coef.est[, 2])
var.est[ , ,2] <- C%*%vcov(mod.DLM)[154:(154 + df), 154:(154 + df)]%*%t(C)

## hypothesis testing on DLM
anova(mod.GLM, mod.DLM, test = "LRT")
U <- svd(crossprod(R))$u
D <- diag(svd(crossprod(R))$d)
Z2 <- as.matrix(dat1[, 9:23])%*%(U[, 1:(nlag - df - 1)])%*%solve(sqrt(D[1:(nlag - df - 1), 1:(nlag - df - 1)]))

Y <- dat1$death
muhat <- exp(predict(mod.DLM))
teststat <- t(Y - muhat)%*%tcrossprod(Z2)%*%(Y - muhat)
M <- t(Z2)%*%diag(muhat)%*%Z2
Eigen <- eigen(M)$values
davies(teststat, lambda = Eigen)$Qq



## EB1
phi <- coef.est[, 1] - coef.est[, 2]
V1 <- var.est[, , 1]
K1.diag <- diag(diag(V1))%*%solve(diag(diag(V1 + tcrossprod(phi))))
coef.est[, 3] <- coef.est[, 1] + K1.diag%*%(coef.est[, 2] - coef.est[, 1])

Gsub3 <- diag(diag(V1)*diag(V1 - tcrossprod(phi))/(diag(V1 + tcrossprod(phi)))^2)
G3 <- cbind(Gsub3, diag(nlag) - Gsub3)

X_GLM <- model.matrix(mod.GLM)
X_DLM <- model.matrix(mod.DLM)

tscore12 <- matrix(0, nlag, df + 1)
for(i in 1:dim(dat1)[1]) {
	tscore12 <- tscore12 + as.numeric((dat1[i, 4] - exp(X_GLM[i, ]%*%coef(mod.GLM)))*(dat1[i, 4] - exp(X_DLM[i, ]%*%coef(mod.DLM))))*
		X_GLM[i, 154:168]%*%t(X_DLM[i, 154:(154 + df)])
}

tscore1 <- matrix(0, nlag, nlag)
for(i in 1:dim(dat1)[1]) {
	tscore1 <- tscore1 + as.numeric((dat1[i, 4] - exp(X_GLM[i, ]%*%coef(mod.GLM)))^2)*tcrossprod(X_GLM[i, 154:168])
}
tscore2 <- matrix(0, df + 1, df + 1)
for(i in 1:dim(dat1)[1]) {
	tscore2 <- tscore2 + as.numeric((dat1[i, 4] - exp(X_DLM[i, ]%*%coef(mod.DLM)))^2)*tcrossprod(X_DLM[i, 154:(154 + df)])
}
i1 <- t(X_GLM[, 154:168])%*%diag(c(exp(X_GLM%*%coef(mod.GLM))))%*%X_GLM[, 154:168]
i2 <- t(X_DLM[, 154:(154 + df)])%*%diag(c(exp(X_DLM%*%coef(mod.DLM))))%*%X_DLM[, 154:(154 + df)]

var.all <- rbind(cbind(solve(i1)%*%tscore1%*%solve(i1), solve(i1)%*%tscore12%*%solve(i2)), 
	cbind(solve(i2)%*%t(tscore12)%*%solve(i1), solve(i2)%*%tscore2%*%solve(i2)))

var.beta <- tran1%*%var.all%*%t(tran1)

V2 <- cbind(diag(nlag), -diag(nlag))%*%var.beta%*%t(cbind(diag(nlag), -diag(nlag)))
K2.diag <- diag(diag(V2))%*%solve(diag(diag(V2 + tcrossprod(phi))))

coef.est[, 4] <- coef.est[, 1] + K2.diag%*%(coef.est[, 2] - coef.est[, 1])

Gsub4 <- diag(diag(V2)*diag(V2 - tcrossprod(phi))/(diag(V2 + tcrossprod(phi)))^2)
G4 <- cbind(Gsub4, diag(nlag) - Gsub4)

var.est[,,3] <- G3%*%var.beta%*%t(G3)
var.est[,,4] <- G4%*%var.beta%*%t(G4)

## GRR
RTR <- crossprod(R)
RTR_expand <- matrix(0, 153+nlag, 153+nlag)
RTR_expand[154:(153+nlag), 154:(153+nlag)] <- RTR
XWX <- t(X_GLM)%*%diag(c(exp(X_GLM%*%coef(mod.GLM))))%*%X_GLM
XW <- t(X_GLM)%*%diag(c(exp(X_GLM%*%coef(mod.GLM))))
XWY <- XW%*%dat1[, 4]

lambda.grid <- seq(1000000, 100000000, 1000000)
aicc.all <- gcv.all <- rep(0, length(lambda.grid))
for(i in 1:length(lambda.grid)) {
	temp1 <- solve(XWX + lambda.grid[i]*RTR_expand)
	temp2 <- temp1%*%XWX
	temp.est <- temp2%*%coef(mod.GLM)	
	temp.trace <- sum(diag(temp2))
	temp.XB <- X_GLM%*%temp.est
	aicc.all[i] <- -2*(sum(-exp(temp.XB)+temp.XB*dat1[, 4])) + 2*temp.trace
	gcv.all[i] <- crossprod(X_GLM%*%temp1%*%XWY - dat1[, 4]) / (1 - temp.trace / n)^2
	# print(c(i, temp.trace))
}
# plot(1:length(lambda.grid), aicc.all, type = "l", pch = 20)
# plot(1:length(lambda.grid), gcv.all, type = "l", pch = 20)

coef.est[, 5] <- (solve(XWX + lambda.grid[which.min(aicc.all)]*RTR_expand)%*%XWX%*%coef(mod.GLM))[154:168]
coef.est[, 6] <- (solve(XWX + lambda.grid[which.min(gcv.all)]*RTR_expand)%*%XWX%*%coef(mod.GLM))[154:168]
var.est[,, 5] <- (solve(XWX + lambda.grid[which.min(aicc.all)]*RTR_expand)%*%XWX%*%solve(XWX + lambda.grid[which.min(aicc.all)]*RTR_expand))[154:168, 154:168]
var.est[,, 6] <- (solve(XWX + lambda.grid[which.min(gcv.all)]*RTR_expand)%*%XWX%*%solve(XWX + lambda.grid[which.min(gcv.all)]*RTR_expand))[154:168, 154:168]


## GRR - choosing parameter to minimize asymptotic MSE
temp_fn <- function(k) {
	sum_mat <- XWX + k*RTR_expand
	sol_sum_mat <- solve(sum_mat)
	gamma1 <- sum(diag(sol_sum_mat%*%XWX%*%sol_sum_mat))
	gamma2 <- k^2*crossprod(sol_sum_mat%*%RTR_expand%*%coef(mod.GLM))
	return(gamma1 + gamma2)
}
# lambda_asymp <- nlm(temp_fn, p = 0)$estimate
lambda.grid2 <- seq(10000, 10000000, 10000)
asym.all <- sapply(lambda.grid2, temp_fn)
coef.est[, 7] <- (solve(XWX + lambda.grid2[which.min(asym.all)]*RTR_expand)%*%XWX%*%coef(mod.GLM))[154:168]
var.est[,, 7] <- (solve(XWX + lambda.grid2[which.min(asym.all)]*RTR_expand)%*%XWX%*%solve(XWX + lambda.grid2[which.min(asym.all)]*RTR_expand))[154:168, 154:168]


## BDLM based on GLM (normal approximation)
grid.len <- 10
omega1.grid <- seq(-0.35, -0.05,, grid.len)
omega2.grid <- seq(-0.37,0,, grid.len)
Omega.all <- solve.Omega.all <- Omega.temp.all <- array(0, dim = c(nlag, nlag, grid.len, grid.len))
det.Omega.all <- matrix(0, grid.len, grid.len)
	for(k in 1:grid.len) {
		for(l in 1:grid.len) {
			Omega.temp.all[,,k, l] <- omega.fn(omega1.grid[k], omega2.grid[l])
			det.Omega.all[k, l] <- det(Omega.temp.all[,,k, l])
			solve.Omega.all[,,k, l] <- solve(Omega.temp.all[,,k, l])
	}
}

theta_hat <- coef(mod.GLM)[154:168]
sigma2 <- 10*vcov(mod.GLM)[154, 154]
# sigma2 <- 0.004^2
Sigma <- vcov(mod.GLM)[154:168, 154:168]

post.mean <- array(0, dim = c(length(omega1.grid), length(omega2.grid), nlag, 1))
post.var <- array(0, dim = c(length(omega1.grid), length(omega2.grid), nlag, nlag))
peta <- matrix(0, length(omega1.grid), length(omega2.grid))
for(i in 1:length(omega1.grid)) { 
	for(j in 1:length(omega2.grid)) {
		omega <- omega.fn(omega1.grid[i], omega2.grid[j])
		post.mean[i, j,,] <- solve(1/sigma2*solve(omega) + solve(Sigma))%*%solve(Sigma)%*%theta_hat
		post.var[i, j,,] <- solve(1/sigma2*solve(omega) + solve(Sigma))
		peta[i, j] <- as.numeric(abs(det(sigma2*omega%*%solve(Sigma)+diag(nlag)))^(-0.5)*
			exp(-0.5*t(theta_hat)%*%(solve(Sigma)-solve(Sigma)%*%solve(solve(Sigma)+1/sigma2*solve(omega))%*%solve(Sigma))%*%theta_hat))
	}
}
weight <- peta/sum(peta)
weighted.post.mean <- array(0, dim = c(length(omega1.grid), length(omega2.grid), nlag, 1))
weighted.post.var <- array(0, dim = c(length(omega1.grid), length(omega2.grid), nlag, nlag))
for(i in 1:length(omega1.grid)) {
	for(j in 1:length(omega2.grid)) {
		weighted.post.mean[i,j,,] <- weight[i, j]*post.mean[i, j,,]
		weighted.post.var[i,j,,] <- weight[i, j]*(tcrossprod(post.mean[i, j,,])+post.var[i, j,,])
	}
}
coef.est[ , 8] <- apply(weighted.post.mean, 3, sum)
var.est[, , 8] <- apply(weighted.post.var, 3:4, sum)-tcrossprod(coef.est[ , 8])


## HB 
acc <- 0
nsim <- 11000
beta_HB_all <- matrix(0, nsim, nlag)
# beta_HB_all[1, ] <- rep(0, nlag)
alpha_HB <- coef(mod.GLM)[1:153]
sigma2_p_all <- rep(0, nsim)
sigma2_p_all[1] <- 0.00001
a0 <- b0 <- 0.001
an <- a0 + (nlag - df)/2
prop_mean <- coef(mod.GLM)[154:168]
prop_var <- vcov(mod.GLM)[154:168, 154:168]
tilda <- 0.5
sol.Sigma <- solve(tilda*prop_var)
beta_HB_all[1, ] <- mvrnorm(n = 1, mu = prop_mean, Sigma = tilda*prop_var)

for(i in 2:nsim) {

	prop_beta <- mvrnorm(n = 1, mu = prop_mean, Sigma = tilda*prop_var)

	log.p1 <- -0.5/sigma2_p_all[i-1]*crossprod(R%*%prop_beta)
	log.p0 <- -0.5/sigma2_p_all[i-1]*crossprod(R%*%beta_HB_all[i - 1, ])
	log.L1 <- -sum(exp(X%*%prop_beta + Z%*%alpha_HB)) + sum((X%*%prop_beta + Z%*%alpha_HB)*Y)
	log.L0 <- -sum(exp(X%*%beta_HB_all[i - 1, ] + Z%*%alpha_HB)) + sum((X%*%beta_HB_all[i - 1, ] + Z%*%alpha_HB)*Y)
	log.J1 <- -0.5*t(prop_beta - prop_mean)%*%sol.Sigma%*%(prop_beta - prop_mean)
	log.J0 <- -0.5*t(beta_HB_all[i-1, ] - prop_mean)%*%sol.Sigma%*%(beta_HB_all[i-1, ] - prop_mean)
	
	r <- exp(log.p1 - log.p0 + log.L1 - log.L0 + log.J0 - log.J1)

	ran.num <- runif(1)	
	if(ran.num < r) {beta_HB_all[i, ] <- prop_beta; acc <- acc + 1}
	else {beta_HB_all[i, ] <- beta_HB_all[i-1, ]}

	# offset.pol <- X%*%beta_HB_all[i, ]
	# temp.mod <- glm(Y ~ -1 + offset(offset.pol) + Z, family = poisson(link = "log"), 
	#  	control = glm.control(epsilon = 1e-10, maxit = 1000), na.action = na.omit)
 	# alpha_HB <- coef(temp.mod)
	
	bn <- b0 + crossprod(R%*%beta_HB_all[i, ])/2
	sigma2_p_all[i] <- rinvgamma(n = 1, shape = an, scale = bn)

	print(c(i, r))
}
coef.est[, 9] <- apply(beta_HB_all[seq(1010, nsim, 10), ], 2, mean)
var.est[, , 9] <- cov(beta_HB_all[seq(1010, nsim, 10), ])


## double shrinkage 1 - BDLM based on GLM (normal approximation)
grid.len <- 10
omega1.grid <- seq(-0.35, -0.05,, grid.len)
omega2.grid <- seq(-0.37,0,, grid.len)
Omega.all <- solve.Omega.all <- Omega.temp.all <- array(0, dim = c(nlag, nlag, grid.len, grid.len))
det.Omega.all <- matrix(0, grid.len, grid.len)
	for(k in 1:grid.len) {
		for(l in 1:grid.len) {
			Omega.temp.all[,,k, l] <- omega.fn(omega1.grid[k], omega2.grid[l])
			det.Omega.all[k, l] <- det(Omega.temp.all[,,k, l])
			solve.Omega.all[,,k, l] <- solve(Omega.temp.all[,,k, l])
	}
}

theta_hat <- coef.est[ , 5]
# sigma2 <- 10*vcov(mod.GLM)[154, 154]
sigma2 <- 10*var.est[1, 1, 5]
Sigma <- var.est[, , 5]

post.mean <- array(0, dim = c(length(omega1.grid), length(omega2.grid), nlag, 1))
post.var <- array(0, dim = c(length(omega1.grid), length(omega2.grid), nlag, nlag))
peta <- matrix(0, length(omega1.grid), length(omega2.grid))
for(i in 1:length(omega1.grid)) { 
	for(j in 1:length(omega2.grid)) {
		omega <- omega.fn(omega1.grid[i], omega2.grid[j])
		post.mean[i, j,,] <- solve(1/sigma2*solve(omega) + solve(Sigma))%*%solve(Sigma)%*%theta_hat
		post.var[i, j,,] <- solve(1/sigma2*solve(omega) + solve(Sigma))
		peta[i, j] <- as.numeric(abs(det(sigma2*omega%*%solve(Sigma)+diag(nlag)))^(-0.5)*
			exp(-0.5*t(theta_hat)%*%(solve(Sigma)-solve(Sigma)%*%solve(solve(Sigma)+1/sigma2*solve(omega))%*%solve(Sigma))%*%theta_hat))
	}
}
weight <- peta/sum(peta)
weighted.post.mean <- array(0, dim = c(length(omega1.grid), length(omega2.grid), nlag, 1))
weighted.post.var <- array(0, dim = c(length(omega1.grid), length(omega2.grid), nlag, nlag))
for(i in 1:length(omega1.grid)) {
	for(j in 1:length(omega2.grid)) {
		weighted.post.mean[i,j,,] <- weight[i, j]*post.mean[i, j,,]
		weighted.post.var[i,j,,] <- weight[i, j]*(tcrossprod(post.mean[i, j,,])+post.var[i, j,,])
	}
}
coef.est[ , 10] <- apply(weighted.post.mean, 3, sum)
var.est[, , 10] <- apply(weighted.post.var, 3:4, sum)-tcrossprod(coef.est[ , 10])

## double shrinkage 2
nfold <- 5
lambda.grid2 <- seq(100, 100000, 100)
beta_hat <- coef.est[ , 5]
Sigma <- var.est[ , , 5]
solSigma <- solve(Sigma)
DD <- omega.fn(-0.25, -0.25)
solveDD <- solve(DD)
ind_fold <-  sample(c(rep(1:3, each = 2134), rep(4:5, each = 2133)))
pred_all <- matrix(0, length(Y), nfold)
mspe_all <- matrix(0, length(lambda.grid2), nfold)

for(i in 1:5) {
	dat_train <- dat2[-which(ind_fold == i), ]
	dat_test <- dat2[which(ind_fold == i), ]
	Z_test <- Z[which(ind_fold == i), ]
	Y_test <- Y[which(ind_fold == i)]
	X_test <- X[which(ind_fold == i), ]

	mod_GLM_train <- glm(formu1, data = dat_train, family = poisson(link = "log"), 
		control = glm.control(epsilon = 1e-10, maxit = 1000), na.action = na.omit)
	coefmod1_train <- coef(mod_GLM_train)[154:168]
	varmod1_train <- vcov(mod_GLM_train)[154:168, 154:168]
	sigma2_train <- (summary(mod_GLM_train)$sigma)^2

	mod_DLM_train <- glm(formu2, data = dat_train, family = poisson(link = "log"), 
		control = glm.control(epsilon = 1e-10, maxit = 1000),na.action = na.omit)

	coefmod2_train <- C%*%(coef(mod_DLM_train)[154:(154+df)])
	varmod2_train <- C%*%vcov(mod_DLM_train)[154:(154+df), 154:(154+df)]%*%t(C)

	phi_train <- coefmod1_train - coefmod2_train
	K1.diag_train <- diag(diag(varmod1_train))%*%solve(diag(diag(varmod1_train + tcrossprod(phi_train))))
	coefmod3_train <- coefmod1_train + K1.diag_train%*%(coefmod2_train - coefmod1_train)

	for(j in 1:length(lambda.grid2)) {
		temp <- solve(solSigma + lambda.grid2[j]*solveDD)
		beta_tilda <- temp%*%solSigma%*%coefmod3_train
		mspe_all[j, i] <- mean((Y_test - exp(X_test%*%beta_tilda + Z_test%*%(coef(mod.GLM)[1:153])))^2)
	}
}
lambda_fin <- lambda.grid2[which.min(apply(mspe_all, 1, mean))]
temp <- solve(solSigma + lambda_fin*solveDD)
beta_tilda <- temp%*%solSigma%*%beta_hat
coef.est[ , 11] <- beta_tilda
var.est[ , , 11] <- temp
# var.est[k,,,11] <- temp%*%solSigma%*%temp


### post-analysis
coef.est2 <- matrix(0, nlag, 12)
var.est2 <- array(0, dim = c(nlag, nlag, 12))
coef.est2[, 1:11] <- coef.est
var.est2[, , 1:11] <- var.est
coef.est2[, 12] <- apply(beta_HB_all[seq(1010, nsim, 10), ], 2, mean)
var.est2[, , 12] <- cov(beta_HB_all[seq(1010, nsim, 10), ])

var.est2_diag <- apply(var.est2, 3, diag)

## cumulative
coef.cum <- apply(coef.est2, 2, cumsum)
var.cum <- matrix(0, nlag, 12)
for(i in 1:nlag) {
	for(j in 1:12) {
		var.cum[i, j] <- rep(1, i)%*%var.est2[(1:i), (1:i), j]%*%rep(1, i)
	}
}

upper <- coef.est2 + qnorm(0.975)*sqrt(apply(var.est2, 3, diag))
lower <- coef.est2 - qnorm(0.975)*sqrt(apply(var.est2, 3, diag))
upper.cum <- coef.cum + qnorm(0.975)*sqrt(var.cum)
lower.cum <- coef.cum - qnorm(0.975)*sqrt(var.cum)

write.csv(coef.est2, file = "NMMAP_coef.est2_7days.csv")
write.csv(var.est2_diag, file = "NMMAP_var.est2_diagv.csv")
write.csv(coef.cum, file = "NMMAP_coef.cum_7days.csv")
write.csv(var.cum, file = "NMMAP_var.cum_7days.csv")
write.csv(upper, file = "NMMAP_upper_7days.csv")
write.csv(lower, file = "NMMAP_lower_7days.csv")
write.csv(upper.cum, file = "NMMAP_upper.cum_7days.csv")
write.csv(lower.cum, file = "NMMAP_lower.cum_7days.csv")

## plotting
coef.est2 <- coef.est
coef.est3 <- coef.est2[, c(1:3, 5, 8:11)]

png('figure1.png')
plot(0:(nlag - 1), coef.est3[, 1], pch = 20, type = "l", lwd = 2, ylim = range(coef.est3), ylab = "Lag Coefficients", xlab = "Lag", 
	main = "Estimated Distributed Lag Function")
for(i in 2:8) {
	lines(0:(nlag - 1), coef.est3[, i], pch = 20, type = "l", lwd = 2, col = i)
}	
abline(h = 0, lty = 2)
legend("topright", legend = c("GLM", "DLM", "EB1", "GRR-AICC", "BDLM", "HB", "HB2-GRR", "HP-GRR"), lwd = 2, col = 1:8)
dev.off()


## cumulative
coef.cum <- apply(coef.est, 2, cumsum)
var.cum <- matrix(0, nlag, 11)
for(i in 1:nlag) {
	for(j in 1:11) {
		var.cum[i, j] <- rep(1, i)%*%var.est[(1:i), (1:i), j]%*%rep(1, i)
	}
}

ind <- c(1:3, 5, 8:11)
coef_est2 <- coef.est[, ind]*1000
coef_std2 <- sqrt(apply(var.est[,, ind], 3, diag))*1000
coef_lower <- coef_est2 - qnorm(0.975)*coef_std2
coef_upper <- coef_est2 + qnorm(0.975)*coef_std2

cum_est <- coef.cum[, ind]*1000
cum_std <- sqrt(var.cum[, ind])*1000
cum_lower <- cum_est - qnorm(0.975)*cum_std
cum_upper <- cum_est + qnorm(0.975)*cum_std

coef_est2[, 5] <- apply(boot_coef[, , 1], 2, median)*1000
coef_lower[, 5] <- apply(boot_coef[, , 1], 2, function(x) quantile(x, 0.025))*1000
coef_upper[, 5] <- apply(boot_coef[, , 1], 2, function(x) quantile(x, 0.975))*1000

cum_est[, 5] <- apply(apply(boot_coef[, , 1], 1, cumsum), 1, median)*1000
cum_lower[, 5] <- apply(apply(boot_coef[, , 1], 1, cumsum), 1, function(x) quantile(x, 0.025))*1000
cum_upper[, 5] <- apply(apply(boot_coef[, , 1], 1, cumsum), 1, function(x) quantile(x, 0.975))*1000

#### cumulative plot

#################################################################################################################
## cumulative plot with CIs
lag_ind <- c(4, 8, 15)
xlim <- length(lag_ind)*(dim(cum_est)[2] + 4)
width <- 0.1
lwd = 2
cex = 1.5
png('figure2.png', width = 600, height = 480)
plot(0, 0, col = "white", xlim = c(0, xlim - 3), ylim = range(cum_lower[lag_ind], cum_upper[lag_ind]) + c(-0.2, 0.2), xlab = "", ylab = "", xaxt='n', main = "Cumulative Effects")
for(i in 1:dim(cum_est)[2]) {
	for(j in seq_along(lag_ind)) {
		pos <- i + j*12 - 12
		lines(pos, cum_est[lag_ind[j], i], col = i, cex = cex, pch = 20, type = "o")
		lines(rep(pos, 2), c(cum_lower[lag_ind[j], i], cum_upper[lag_ind[j], i]), pch = 20, col = i, lwd = lwd)
		lines(pos + width*c(-1, 1), rep(cum_lower[lag_ind[j], i], 2), pch = 20, col = i, lwd = lwd)
		lines(pos + width*c(-1, 1), rep(cum_upper[lag_ind[j], i], 2), pch = 20, col = i, lwd = lwd)
	}
}
abline(h = 0, lty = 2, col = "brown", lwd = lwd)
legend("bottomleft", col = 1:dim(cum_est)[2], lwd = 2, legend = c("GLM", "DLM", "EB1", "GRR-AICC", "BDLM", "HB", "HB2-GRR", "HP-GRR"))
axis(side = 1, at = seq(4.5,, 12, 3), labels = paste("Lag ", 0, "-", lag_ind - 1, sep = ""))
dev.off()


## overall plot
plot(0:14, coef.est[, 1], type = "o", pch = 20, main = "Lag Effects - pm10", ylab = "Lag Effects", 
	xlab = "Lags", ylim = range(coef.est), lwd = 2)
for(i in c(3, 10:11)) {
	lines(0:14, coef.est[, i], type = "o", pch = 20, lwd = 2, col = i)
}
legend("topright", legend = c("GLM", "DLM", "EB1", "EB2", "GRR-AICC", "GRR-GCV", "GRR-Asym", "BDLM", "HB"), col = 1:4, lwd = 2)
abline(h = 0)



## overall plot
setwd('C:/Users/yinhsiuc/Desktop/BIOSTAT-Yin-Hsiu Chen/Bhramar/2015 JSM/')
temp_ind <- c(2:3, 5, 8:9)
png('JSM2015_NMMAPplot.png')
plot(0:14, coef.est[, 1], type = "o", pch = 20, main = "Lag Effects - pm10", ylab = "Lag Effects", 
	xlab = "Lags", ylim = range(coef.est), lwd = 2)
for(i in seq_along(temp_ind)) {
	lines(0:14, coef.est[, temp_ind[i]], type = "o", pch = 20, lwd = 2, col = i + 1)
}
legend("topright", legend = c("GLM", "DLM", "EB", "eGRR", "BDLM", "HB"), col = 1:6, lwd = 2)
abline(h = 0)
dev.off()

## cumulative
coef.cum <- apply(coef.est, 2, cumsum)
var.cum <- matrix(0, nlag, 11)
for(i in 1:nlag) {
	for(j in 1:11) {
		var.cum[i, j] <- rep(1, i)%*%var.est[(1:i), (1:i), j]%*%rep(1, i)
	}
}
cum_est <- coef.cum[c(8, 15), c(1:3, 5, 8:9)]*10000
cum_std <- sqrt(var.cum[c(8, 15), c(1:3, 5, 8:9)])*10000
cum_est - qnorm(0.975)*cum_std
cum_est + qnorm(0.975)*cum_std

################################################################################################################
## HB - Laplace Approximation
nsim <- 11000
beta_HB_all2 <- matrix(0, nsim, nlag)
# beta_HB_all[1, ] <- rep(0, nlag)
alpha_HB <- coef(mod.GLM)[1:153]
sigma2_p_all2 <- rep(0, nsim)
sigma2_p_all2[1] <- 0.000000001
a0 <- b0 <- 0.000000001
an <- a0 + (nlag - df)/2
X <- as.matrix(dat1[, 9:23])
Z <- as.matrix(model.matrix(mod.GLM)[, 1:153])
Y <- dat1[, 4]
tempfn <- function(beta, parms) {
		t(X)%*%(Y - exp(X%*%beta + Z%*%alpha_HB)) - RTR%*%beta / parms
}
for(i in 2:nsim) {

	temp_root <- multiroot(f = tempfn, maxiter = 1000, ctol = 1e-10, atol = 1e-10, start = rnorm(nlag, sd = 0.0001), parms = sigma2_p_all2[i-1])	
	temp_mean <- temp_root$root
	temp_Sigma_inv <- t(X)%*%diag(as.numeric(exp(X%*%temp_mean + Z%*%alpha_HB)))%*%X + RTR / sigma2_p_all2[i-1]

	beta_HB_all2[i, ] <- mvrnorm(n = 1, mu = temp_mean, Sigma = solve(temp_Sigma_inv))
	
	bn <- b0 + crossprod(R%*%beta_HB_all2[i, ])/2
	sigma2_p_all2[i] <- rinvgamma(n = 1, shape = an, scale = bn)

	print(i)
}
coef.est[, 9] <- apply(beta_HB_all2[seq(1010, nsim, 10), ], 2, mean)
var.est[, , 9] <- cov(beta_HB_all2[seq(1010, nsim, 10), ])


###################################################################################################################################



