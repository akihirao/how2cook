#Sim_euler_SPiCT.R

library(spict)#ライブラリーの読み込み
library(frapmfun)
library(tidyverse)

setwd("/Users/akihirao/Desktop/水研/Production.model/SPiCT/sim_dteuler")
data_raw <- read.csv("example1.csv") #データの読み込みの指定
inp  <- get_spict_data(data_raw)

inp$priors$logn <- c(log(1.19), 0.001, 1)
inp$priors$logr <- c(log(0.4), 0.0001, 1)
inp$priors$logK <- c(log(1000), 0.001, 1)
inp$priors$logm <- c(log(134.5543), 0.001, 1)
inp$priors$logbkfrac <- c(log(0.8), 5, 1)
inp$priors$logsdf <- c(log(0.3), 10, 0)
inp$priors$logsdc <- c(log(0.001), 0.001, 1)
inp$priors$logbeta <- c(log(1), 2, 0)
inp$priors$logsdb <- c(log(0.3), 10, 0)
inp$priors$logsdi <- list(c(log(0.3), 10, 0),c(log(0.3), 10, 0))
inp$priors$logalpha <- c(log(1), 2, 0)
set_mapsdi <- 0
mapsdi <- c(1, 2)
inp$priors$logq <- list(c(log(0.8), 2, 0), c(log(0.8), 2, 0))


inp <- check.inp(inp)
inp_sim <- inp

K_true <- 1000　#環境収容力の真の値
r_true <- 0.4　#内的自然増加率の真の値
n_true <- 1.19　#形状パラメータの真の値
m_true <- 134.5543
q_true <- 0.8 #漁具能率
sdb_true <- 0.3　# 資源量のプロセス誤差
sdf_true <- 0.3 # 漁獲率のプロセス誤差

inp$dteuler <- 1/16
inp_sim$ini$logn <- log(n_true)
inp_sim$ini$logK <- log(K_true)
inp_sim$ini$logm <- log(m_true)
inp_sim$ini$logr <- log(r_true)
inp_sim$ini$logsdb <- log(sdb_true)
inp_sim$ini$logsdf <- log(sdf_true)

#sim_r <- data.frame(r=vector(),sim=vector(),dteulter=vector())
#sim_K <- data.frame(K=vector(),sim=vector(),dteulter=vector())
no_sim <- 20
res_all <- as.list(numeric(no_sim))
r_mat <- matrix(ncol = no_sim, nrow = 5)
K_mat <- matrix(ncol = no_sim, nrow = 5)
m_mat <- matrix(ncol = no_sim, nrow = 5)
rownames(r_mat) <- c("e16","e8","e4","e2","e1")
m_mat <- matrix(ncol = no_sim, nrow = 5)


for(i in 1:no_sim){
  sim <- sim.spict(inp_sim)
  sim$priors$logn <- c(log(1.19),0.0001,1)
  sim$priors$logr <- c(log(0.5), 0.5, 1)
  sim$priors$logK <- c(log(1000), 1, 1)
  sim$priors$logm <- c(log(1e-04), 1, 0)
  
  sim_e16 <- sim
  sim_e16$dteuler <- 1/16
  sim_e8 <- sim
  sim_e8$dteuler <- 1/8
  sim_e4 <- sim
  sim_e4$dteuler <- 1/4
  sim_e2  <- sim
  sim_e2$dteuler <- 1/2
  sim_e1  <- sim
  sim_e1$dteuler <- 1
  
  sim_eulter_series <- list(sim_e16,sim_e8,sim_e4,sim_e2,sim_e1)
  
  res_euler_series <- as.list(numeric(5))
  for(j in 1:5){
     res_euler_series[[j]] <- fit.spict(sim_eulter_series[[j]])
     r_pos <- match("r",names(res_euler_series[[j]]$value))
     K_pos <- match("K",names(res_euler_series[[j]]$value))
     m_pos <- match("m",names(res_euler_series[[j]]$value))
     r_mat[j,i] <- c(sim$r,res_euler_series[[j]]$value[r_pos][[1]])
     K_mat[j,i] <- c(sim$K,res_euler_series[[j]]$value[K_pos][[1]])
     m_mat[j,i] <- c(sim$K,res_euler_series[[j]]$value[m_pos][[1]])
  }
  res_all[[i]] <- res_euler_series
}  

sim_series_lab <- paste("sim",seq(1:no_sim),sep="")
rownames(r_mat) <- c("e16","e8","e4","e2","e1")
colnames(r_mat) <- sim_series_lab
rownames(K_mat) <- c("e16","e8","e4","e2","e1")
colnames(K_mat) <- sim_series_lab
rownames(m_mat) <- c("e16","e8","e4","e2","e1")
colnames(m_mat) <- sim_series_lab

par(mfrow=c(3,1))
r_sim_dat <- r_mat %>% as.data.frame() %>% tibble::rownames_to_column(var = "euler") %>% tidyr::pivot_longer(col = -euler, names_to = "sim", values_to = "r")
r_sim_dat$euler <- factor(r_sim_dat$euler,levels=c("e16","e8","e4","e2","e1"))
plot(r_sim_dat$r~r_sim_dat$euler,xlab="Time step",ylab ="r")

K_sim_dat <- K_mat %>% as.data.frame() %>% tibble::rownames_to_column(var = "euler") %>% tidyr::pivot_longer(col = -euler, names_to = "sim", values_to = "K")
K_sim_dat$euler <- factor(K_sim_dat$euler,levels=c("e16","e8","e4","e2","e1"))
plot(K_sim_dat$K~K_sim_dat$euler,xlab="Time step",ylab ="K")

m_sim_dat <- m_mat %>% as.data.frame() %>% tibble::rownames_to_column(var = "euler") %>% tidyr::pivot_longer(col = -euler, names_to = "sim", values_to = "m")
m_sim_dat$euler <- factor(m_sim_dat$euler,levels=c("e16","e8","e4","e2","e1"))
plot(m_sim_dat$m~m_sim_dat$euler,xlab="Time step",ylab ="m")

par(mfrow=c(1,1))
library(tagcloud)
col_sim <- smoothPalette(no_sim, pal="RdBu")

#euler_series <- c(1/16,1/8,1/4,1/2,1)
euler_series <- c(0,1,2,3,4)
r_sim_dat_start <- r_sim_dat %>% filter(sim==sim_series_lab[2])
plot(euler_series, r_sim_dat_start$r, xlim=c(0,4), ylim=c(0,1.2),xlab="Time series", ylab="r", type="l",col=col_sim [1])
par(new=T)
plot(euler_series, r_sim_dat_start$r, xlim=c(0,4), ylim=c(0,1.2),xlab="", ylab="", pch=20, cex=0.3,col=col_sim [1])

#i <- 2 # count for 
for(k in 1:no_sim){
  r_sim_dat_plot <- r_sim_dat %>% filter(sim==sim_series_lab[k])
  par(new=T)
  plot(euler_series, r_sim_dat_plot$r, xlim=c(0,4), ylim=c(0,1.2),xlab="Time series", ylab="r", type="l")
  par(new=T)
  plot(euler_series, r_sim_dat_plot$r, xlim=c(0,4), ylim=c(0,1.2),xlab="", ylab="", pch=20, cex=0.3)
}

