## Simulations illustrating the effect of dteuler value
## Author: Akira Hirao (hirao_akira51@fra.go.jp)
## Date: 2022/06/30

#### Sim_spict_dteuler_contrast_vs_non_contrast.R

#### loading library
library(spict)
library(tidyverse)
library(gridExtra)
library(ggbeeswarm)


#### setting parmeters for simulations
K_true <- 1000  #true K (carring capacity)
r_true <- 0.4 #true r (intrinsic growth rate)
n_true <- 1.19　#true n (shape parameter)
m_true <- r_true*K_true/(n_true^(n_true/(n_true-1)))
q_true <- 0.1 #true catchability
F0 <- 0.3
bkfrac_true <- 0.7 #true initial depletion rate if "norm" type of simulation
sigma_pro <- 0.3　# process error of biomass
sigma_obs <- 0.3　# observed error of cpue
sdf_true <- 0.3 # process error of catch but not used in the simulation
Bmsyd_true <- K_true*n_true^(1/(1 - n_true))
Fmsyd_true <- m_true/Bmsyd_true


range_year <- c(1991,2020)


#### generating simulation data-set based on the discrete production model

sim.dif.pdm <- function(type, range_year, K, r, n, q, sigma_pro, sigma_obs){
  if(type=="non_contrast_type"){ 
    #non_contrast_type: F follows a normal distribution (mean 0.3 +- 0.1 SD)
    years <- seq(range_year[1],range_year[2],by=1)
    nyears <- length(years)
    
    F_vec <- B_vec <- C_vec <- NULL
    B_vec[1] <- K*bkfrac_true # determining initial biomass
    F_vec[1] <- exp(log(F0)+ rnorm(1,0,0.1))
    
    for(t in 2:nyears){
      t_pre <- t - 1
      C_vec[t_pre] <- B_vec[t_pre]*F_vec[t_pre]
      B_vec[t] <- (B_vec[t_pre] + r/(n - 1)*B_vec[t_pre]*(1 - (B_vec[t_pre]/K)^(n - 1)) - F_vec[t_pre] * B_vec[t_pre])*exp(rnorm(1,-sigma_pro^2/2,sigma_pro))
      F_vec[t] <- exp(log(F0)+ rnorm(1,0,0.1)) 
    }
    C_vec[nyears] <- B_vec[nyears]*F_vec[nyears]
    cpue <- exp(rnorm(nyears, log(q*B_vec), sigma_obs))
    
    inp <- list(obsC=C_vec,timeC=years, obsI= cpue,timeI=years)
    BlastBmsyd <- B_vec[nyears]/Bmsyd_true
    FlastFmsyd <- F_vec[nyears]/Fmsyd_true
    
    sim_data <- list(inp, BlastBmsyd, FlastFmsyd)
    return(sim_data)
    
  }else if(type=="contrast_type"){
    #contrast_type: F <- Fmsy*1.2 (during 1st to 20th years); F <- Fmsy*0.5 during 21 th to 30th years
    years <- seq(range_year[1],range_year[2],by=1)
    nyears <- length(years)
    
    F_vec <- B_vec <- C_vec <- NULL
    B_vec[1] <- K*bkfrac_true # determining initial biomass
    F_vec[1] <- exp(log(Fmsyd_true*1.2)+ rnorm(1,0,0.1))
    
    for(t in 2:nyears){
      if(t <= 20){
        t_pre <- t - 1
        C_vec[t_pre] <- B_vec[t_pre]*F_vec[t_pre]
        B_vec[t] <- (B_vec[t_pre] + r/(n - 1)*B_vec[t_pre]*(1 - (B_vec[t_pre]/K)^(n - 1)) - F_vec[t_pre] * B_vec[t_pre])*exp(rnorm(1,-sigma_pro^2/2,sigma_pro))
        F_vec[t] <- exp(log(Fmsyd_true*1.2)+ rnorm(1,0,0.1))
      }else{
        t_pre <- t - 1
        C_vec[t_pre] <- B_vec[t_pre]*F_vec[t_pre]
        B_vec[t] <- (B_vec[t_pre] + r/(n - 1)*B_vec[t_pre]*(1 - (B_vec[t_pre]/K)^(n - 1)) - F_vec[t_pre] * B_vec[t_pre])*exp(rnorm(1,-sigma_pro^2/2,sigma_pro))
        F_vec[t] <- exp(log(Fmsyd_true*0.5)+ rnorm(1,0,0.1))
      }
    }
    C_vec[nyears] <- B_vec[nyears]*F_vec[nyears]
    cpue <- exp(rnorm(nyears, log(q*B_vec), sigma_obs))
    
    inp <- list(obsC=C_vec,timeC=years, obsI= cpue,timeI=years)
    BlastBmsyd <- B_vec[nyears]/Bmsyd_true
    FlastFmsyd <- F_vec[nyears]/Fmsyd_true
    
    sim_data <- list(inp, BlastBmsyd, FlastFmsyd)
    return(sim_data)
    
  }else{
    
  }
}


#========================================
### executing simulations
no_sim <- 100
euler_value_series <- c(1/64,1/16,1/4,1/1)
no_euler_series <- length(euler_value_series)

res_out <- as.list(numeric(no_sim))


sim_convergence_matrix <- matrix(0,ncol=(no_euler_series+1))

colnames(sim_convergence_matrix) <- c("sim","1/64","1/16","1/4","1")
sim_count <- 0
k <- 1
while (k <= no_sim){
  #generating simulation data with "contrast_type"
  sim_data <- sim.dif.pdm("contrast_type",range_year,K_true, r_true, n_true,q_true,sigma_pro, sigma_obs)
  inp <- sim_data[[1]]
  inp$priors$logn <- c(log(1.19),0.0001,1)
  inp$priors$logr <- c(log(0.5),1, 1)
  
  convergence_vec <- numeric(no_euler_series)
  res <- as.list(numeric(no_euler_series))
  for(i in 1:no_euler_series){
    inp$dteuler <- euler_value_series[i]
    inp$stabilise <- set_stabilise <- 1
    inp <- check.inp(inp)
    res_out_temp <- fit.spict(inp)
    res[[i]] <- res_out_temp
    convergence <- res_out_temp$opt$convergence #checking convergence
    finiteness <- all(is.finite(res_out_temp$sd))
    
    if(convergence == 0 && finiteness == TRUE){
      convergence_vec[i] <- 1
    }else{
      convergence_vec[i] <- 0
    } 
  }
  sim_count <- sim_count + 1
  sim_convergence_vec <- c(sim_count, convergence_vec)
  sim_convergence_matrix <- rbind(sim_convergence_matrix,sim_convergence_vec)
  
  if(sum(convergence_vec) >= no_euler_series){
    BlBmsyd_true <- sim_data[[2]]
    FlFmsyd_true <- sim_data[[3]]  
    res_out[[k]] <- list(res,BlBmsyd_true,FlFmsyd_true)
    k <- k + 1
  }else{
    k <- k
  }
}

#========================================
# Plotting convergence rate
sim_convergence_matrix <- sim_convergence_matrix[-1,]
df_sim_convergence <- data.frame(sim=sim_convergence_matrix[,1],dteuler64=sim_convergence_matrix[,2],dteuler16=sim_convergence_matrix[,3],dteuler4=sim_convergence_matrix[,4],dteuler1=sim_convergence_matrix[,5])
sim_convergence <- df_sim_convergence %>% pivot_longer(cols = c(dteuler64,dteuler16,dteuler4,dteuler1),names_to ="dteuler", values_to="convergence")
sim_convergence$dteuler <- sub("dteuler64", "1/64", sim_convergence$dteuler)
sim_convergence$dteuler <- sub("dteuler16", "1/16", sim_convergence$dteuler)
sim_convergence$dteuler <- sub("dteuler4", "1/4", sim_convergence$dteuler)
sim_convergence$dteuler <- sub("dteuler1", "1", sim_convergence$dteuler)

sim_convergence$dteuler <- factor(sim_convergence$dteuler,levels=c("1/64","1/16","1/4","1"))

convergence_summary <- sim_convergence %>% group_by(dteuler) %>% dplyr::summarise(convergence_mean=mean(convergence))

convergence_barplot <- ggplot(convergence_summary, aes(x=dteuler, y=convergence_mean, fill=dteuler)) + labs(x = "Time step", y = "Convergence rate")
convergence_barplot <- convergence_barplot + geom_bar(stat = "identity") + coord_cartesian(ylim = c(0,1))

plot(convergence_barplot)
#========================================



###========================================
####summarizing parameters

#### output of results
r_estimate <- vector()
K_estimate <- vector()
m_estimate <- vector()
Fmsys_estimate <- vector()
Bmsys_estimate <- vector()
MSYs_estimate <- vector()
sdb_estimate <- vector()
bkfrac_estimate <- vector()
sim_lab <- vector()
time_step <- vector()
BlBmsy_estimate <- vector()
FlFmsy_estimate <- vector()
BlBmsyd_true <- vector()
FlFmsyd_true <- vector()
biased_BlBmsyd <- vector()
biased_FlFmsyd <- vector()

count_vec <- 1
for(m in 1:no_sim){
  
  res_output <- res_out[[m]][[1]]
  
  for(n in 1:no_euler_series){
    
    pos_r <- match("r",names(res_output[[n]]$value))
    pos_K <- match("K",names(res_output[[n]]$value))
    pos_m <- match("m",names(res_output[[n]]$value))
    pos_sdb <- match("sdb",names(res_output[[n]]$value))
    pos_Bmsys <- match("Bmsys",names(res_output[[n]]$value))
    pos_Fmsys <- match("Fmsys",names(res_output[[n]]$value))
    pos_MSYs <- match("MSYs",names(res_output[[n]]$value))
    pos_logbkfrac <- match("logbkfrac", names(res_output[[n]]$value))
    
    pos_logBlBmsy <- match("logBlBmsy", names(res_output[[n]]$value))
    pos_logFlFmsy <- match("logFlFmsy", names(res_output[[n]]$value))
    
    r_estimate[count_vec] <- res_output[[n]]$value[pos_r][[1]]
    K_estimate[count_vec] <- res_output[[n]]$value[pos_K][[1]]
    m_estimate[count_vec] <- res_output[[n]]$value[pos_m][[1]]
    sdb_estimate[count_vec] <- res_output[[n]]$value[pos_sdb][[1]]
    Bmsys_estimate[count_vec] <- res_output[[n]]$value[pos_Bmsys][[1]]
    Fmsys_estimate[count_vec] <- res_output[[n]]$value[pos_Fmsys][[1]]
    MSYs_estimate[count_vec] <- res_output[[n]]$value[pos_MSYs][[1]]
    bkfrac_estimate[count_vec] <- exp(res_output[[n]]$value[pos_logbkfrac][[1]])
    
    BlBmsy_estimate[count_vec] <- exp(res_output[[n]]$value[pos_logBlBmsy][[1]])
    FlFmsy_estimate[count_vec] <- exp(res_output[[n]]$value[pos_logFlFmsy][[1]])
    BlBmsyd_true[count_vec] <- res_out[[m]][[2]]
    FlFmsyd_true[count_vec] <- res_out[[m]][[3]]
    
    sim_lab[count_vec] <- m
    time_step[count_vec] <- euler_value_series[n]
    
    count_vec <- count_vec + 1
  }
}



param <- as.tibble(data.frame(r=r_estimate, K=K_estimate, m=m_estimate, sdb=sdb_estimate, Bmsys=Bmsys_estimate, Fmsys=Fmsys_estimate, MSYs=MSYs_estimate, BlBmsy=BlBmsy_estimate, BlBmsy_true = BlBmsyd_true,FlFmsy=FlFmsy_estimate,FlFmsy_true=FlFmsyd_true, bkfrac= bkfrac_estimate, sim=sim_lab, time_step_value=time_step))
param <- param %>% mutate(deviate_BlBmsy=((BlBmsy-BlBmsy_true)/BlBmsy_true),deviate_FlFmsy=((FlFmsy-FlFmsy_true)/FlFmsy_true),diff_BlBmsy=(BlBmsy-BlBmsy_true),diff_FlFmsy=(FlFmsy-FlFmsy_true)) %>% mutate(time_step=time_step_value)

param$time_step <- sub(0.015625, "1/64",param$time_step)
param$time_step <- sub("0.0625", "1/16",param$time_step)
param$time_step <- sub("0.25", "1/4",param$time_step)
param$time_step <- sub("1", "1",param$time_step)
param$time_step <- factor(param$time_step, levels=c("1/64","1/16","1/4","1"))

param_summary <- param %>% group_by(time_step) %>% dplyr::summarise(r_median=format(round(median(r),2),nsmall=2), 
                                                                    K_median=format(round(median(K),1),nsmall=1),
                                                                    m_median=format(round(median(m),1),nsmall=1),
                                                                    sdb_median=format(round(median(sdb),3),nsmall=3),
                                                                    Bmsys_median=format(round(median(Bmsys),1),nsmall=1),
                                                                    Fmsys_median=format(round(median(Fmsys),3),nsmall=3),
                                                                    MSYs_median=format(round(median(MSYs),1),nsmall=1),
                                                                    bkfrac_median=format(round(median(bkfrac),2),nsmall=2),
                                                                    deviate_BlBmsy_median=format(round(median(deviate_BlBmsy, na.rm=T),1),nsmall=1),
                                                                    deviate_FlFmsy_median=format(round(median(deviate_FlFmsy),1),nsmall=1), 
                                                                    r_CV=format(round(sd(r)/mean(r),3),nsmall=3),
                                                                    K_CV=format(round(sd(K)/mean(K),3),nsmall=3),
                                                                    m_CV=format(round(sd(m)/mean(m),3),nsmall=3),
                                                                    sdb_CV=format(round(sd(sdb)/mean(sdb),3),nsmall=3),
                                                                    Bmsys_CV=format(round(sd(Bmsys)/mean(Bmsys),3),nsmall=3),
                                                                    Fmsys_CV=format(round(sd(Fmsys)/mean(Fmsys),3),nsmall=3),
                                                                    MSYs_CV=format(round(sd(MSYs)/mean(MSYs),3),nsmall=3),
                                                                    bkfrac_CV=format(round(sd(bkfrac)/mean(bkfrac),3),nsmall=3),
                                                                    deviate_BlBmsy_CV=format(round(sd(deviate_BlBmsy, na.rm=T)/mean(deviate_BlBmsy, na.rm=T),3),nsmall=3),
                                                                    deviate_FlFmsy_CV=format(round(sd(deviate_FlFmsy, na.rm=T)/mean(deviate_FlFmsy, na.rm=T),3),nsmall=3)
                                                                    )
print(param_summary)


#==============================================================================
####Violinplot for each of parameters
r_violinplot <- ggplot(param, aes(x=time_step, y=r, color=time_step)) + labs(x = "Time step")
r_violinplot <- r_violinplot + geom_violin() 
r_violinplot <- r_violinplot + geom_quasirandom(size=0.01, alpha=0.5) + geom_hline(yintercept = r_true, linetype="dashed", col="red")

K_violinplot <- ggplot(param, aes(x=time_step, y=K, color=time_step)) + labs(x = "Time step")
K_violinplot <- K_violinplot + geom_violin() 
K_violinplot <- K_violinplot + geom_quasirandom(size=0.01, alpha=0.5) + geom_hline(yintercept = K_true, linetype="dashed", col="red")

m_violinplot <- ggplot(param, aes(x=time_step, y=m, color=time_step)) + labs(x = "Time step")
m_violinplot <- m_violinplot + geom_violin() 
m_violinplot <- m_violinplot + geom_quasirandom(size=0.01, alpha=0.5) + geom_hline(yintercept = m_true, linetype="dashed", col="red")

bkfrac_violinplot <- ggplot(param, aes(x=time_step, y=bkfrac, color=time_step)) + labs(x = "Time step")
bkfrac_violinplot <- bkfrac_violinplot + geom_violin() 
bkfrac_violinplot <- bkfrac_violinplot + geom_quasirandom(size=0.01, alpha=0.5) + geom_hline(yintercept = bkfrac_true, linetype="dashed", col="red")

Bmsys_violinplot <- ggplot(param, aes(x=time_step, y=Bmsys, color=time_step)) + labs(x = "Time step")
Bmsys_violinplot <- Bmsys_violinplot + geom_violin() 
Bmsys_violinplot <- Bmsys_violinplot + geom_quasirandom(size=0.01, alpha=0.5)

Fmsys_violinplot <- ggplot(param, aes(x=time_step, y=Fmsys, color=time_step)) + labs(x = "Time step")
Fmsys_violinplot <- Fmsys_violinplot + geom_violin() 
Fmsys_violinplot <- Fmsys_violinplot + geom_quasirandom(size=0.01, alpha=0.5)

diff_BlBmsys_violinplot <- ggplot(param, aes(x=time_step, y=diff_BlBmsy, color=time_step)) + labs(x = "Time step", y="B/Bmsy(Estimate-true)")
diff_BlBmsys_violinplot <- diff_BlBmsys_violinplot + geom_violin() 
diff_BlBmsys_violinplot <- diff_BlBmsys_violinplot + geom_quasirandom(size=0.01, alpha=0.5)

diff_FlFmsys_violinplot <- ggplot(param, aes(x=time_step, y=diff_FlFmsy, color=time_step)) + labs(x = "Time step", y="F/Fmsy(Estimate-true)")
diff_FlFmsys_violinplot <- diff_FlFmsys_violinplot + geom_violin() 
diff_FlFmsys_violinplot <- diff_FlFmsys_violinplot + geom_quasirandom(size=0.01, alpha=0.5)


deviate_BlBmsys_violinplot <- ggplot(param, aes(x=time_step, y=deviate_BlBmsy, color=time_step)) + labs(x = "Time step", y="B/Bmsy [(Estimate-true)/true]")
deviate_BlBmsys_violinplot <- deviate_BlBmsys_violinplot + geom_violin() 
deviate_BlBmsys_violinplot <- deviate_BlBmsys_violinplot + geom_quasirandom(size=0.01, alpha=0.5)

deviate_FlFmsys_violinplot <- ggplot(param, aes(x=time_step, y=deviate_FlFmsy, color=time_step)) + labs(x = "Time step", y="F/Fmsy [(Estimate-true)/true]")
deviate_FlFmsys_violinplot <- deviate_FlFmsys_violinplot + geom_violin() 
deviate_FlFmsys_violinplot <- deviate_FlFmsys_violinplot + geom_quasirandom(size=0.01, alpha=0.5)

params_violinplot <- grid.arrange(m_violinplot, r_violinplot, K_violinplot, bkfrac_violinplot,ncol=2)
plot(params_violinplot)

#deviate_BlBmsy_FlFmsy_violinplot <- grid.arrange(diff_BlBmsys_violinplot, deviate_BlBmsys_violinplot, diff_FlFmsys_violinplot, deviate_FlFmsys_violinplot,ncol=2)
deviate_BlBmsy_FlFmsy_violinplot <- grid.arrange(deviate_BlBmsys_violinplot, deviate_FlFmsys_violinplot,ncol=1)
plot(deviate_BlBmsy_FlFmsy_violinplot)

#BlBmsy_FlFmsy_summary <- param %>% group_by(time_step) %>% dplyr::summarise(deviate_BlBmsy_median=format(round(median(deviate_BlBmsy,na.rm=TRUE),2),nsmall=2), deviate_FlFmsy_median=format(round(median(deviate_FlFmsy,na.rm=TRUE),2),nsmall=2),deviate_BlBmsy_CV=format(round(sd(deviate_BlBmsy,na.rm=TRUE)/mean(deviate_BlBmsy,na.rm=TRUE),3),nsmall=3),deviate_FlFmsy_CV=format(round(sd(deviate_FlFmsy,na.rm=TRUE)/mean(deviate_FlFmsy,na.rm=TRUE),3),nsmall=3))
BlBmsy_FlFmsy_summary <- param %>% group_by(time_step) %>% dplyr::summarise(deviate_BlBmsy_median=format(round(median(deviate_BlBmsy,na.rm=TRUE),2),nsmall=2), deviate_FlFmsy_median=format(round(median(deviate_FlFmsy,na.rm=TRUE),2),nsmall=2),deviate_BlBmsy_CV=format(round(sd(deviate_BlBmsy,na.rm=TRUE)/mean(deviate_BlBmsy,na.rm=TRUE),3),nsmall=3))
print(BlBmsy_FlFmsy_summary)

BlBmsy_FlFmsy_potinplot <- ggplot(param, aes(x=deviate_BlBmsy,y=deviate_FlFmsy, color = time_step)) + geom_point() + labs(x = "B/Bmsy [(Estimate-true)/true]", y="F/Fmsy [(Estimate-true)/true]")
plot(BlBmsy_FlFmsy_potinplot)

#==============================================================================
####correlation between estimated parameters (m, r, K, bkfrac, and sigma_pro)
corr_eqn <- function(x,y, digits = 2) {
  corr_coef <- round(cor(x, y), digits = digits)
  corr_p_value <- cor.test(x, y)$p.value
  paste("italic(r) == ", corr_coef,"~';'~italic(P)==", corr_p_value)
}

r_K_labels = data.frame(x =1.4, y = 14000, label = corr_eqn(param$r, param$K))
r_K_cor_out <- cor.test(param$r,param$K)
r_K_plot <- ggplot(param) +geom_point(aes(x=r,y=K,color=time_step)) + labs(x="r") + labs(y="K") + geom_text(data = r_K_labels, aes(x = x, y = y,label = label), parse = TRUE) 

r_sdb_labels = data.frame(x =0.4, y = 3.5, label = corr_eqn(param$r, param$sdb))
r_sdb_cor_out <- cor.test(param$r,param$sdb)
r_sdb_plot <- ggplot(param) +geom_point(aes(x=sdb,y=r,color=time_step)) + labs(x="sigma_B", y="r") + geom_text(data = r_sdb_labels, aes(x = x, y = y,label = label), parse = TRUE) 

K_sdb_labels = data.frame(x =0.4, y = 14000, label = corr_eqn(param$sdb, param$K))
K_sdb_cor_out <- cor.test(param$sdb,param$K)
K_sdb_plot <- ggplot(param) +geom_point(aes(x=sdb,y=K,color=time_step)) + labs(x="sigma_B") + labs(y="K") + geom_text(data = K_sdb_labels, aes(x = x, y = y,label = label), parse = TRUE) 


bk_sdb_labels = data.frame(x =0.4, y = 2.5, label = corr_eqn(param$sdb, param$bkfrac))
bk_sdb_cor_out <- cor.test(param$sdb,param$bkfrac)
bk_sdb_plot <- ggplot(param) +geom_point(aes(x=sdb,y=bkfrac,color=time_step)) + labs(x="sigma_B", y="bkfrac") + geom_text(data = bk_sdb_labels, aes(x = x, y = y,label = label), parse = TRUE) 

m_sdb_labels = data.frame(x =0.4, y = 4000, label = corr_eqn(param$sdb, param$m))
m_sdb_cor_out <- cor.test(param$sdb,param$m)
m_sdb_plot <- ggplot(param) +geom_point(aes(x=sdb,y=m,color=time_step)) + labs(x="sigma_B") + labs(y="m") + geom_text(data = m_sdb_labels, aes(x = x, y = y,label = label), parse = TRUE) 


K_bk_labels = data.frame(x =6000, y = 2.5, label = corr_eqn(param$K, param$bkfrac))
K_bk_cor_out <- cor.test(param$K,param$bkfrac)
K_bk_plot <- ggplot(param) +geom_point(aes(x=K,y=bkfrac,color=time_step)) + labs(x="K", y="bkfrac") + geom_text(data = K_bk_labels, aes(x = x, y = y,label = label), parse = TRUE) 

cor_combine_plot <- grid.arrange(r_K_plot, r_sdb_plot, K_sdb_plot, m_sdb_plot)
plot(cor_combine_plot)

cor_bk_K_combine_plot <- grid.arrange(bk_sdb_plot, K_bk_plot,ncol=2)
plot(cor_bk_K_combine_plot)
