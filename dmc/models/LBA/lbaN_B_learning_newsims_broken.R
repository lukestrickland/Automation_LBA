# Specifying a PM 'trigger failure' model where the PM accumulator
#is not present on a percentage of trials. Reactive inhibition of
#ongoing task accumulation must also fail on those same trials.
#This requires substituting the nonPM accumulation rate into
#pm trials.
#Edited simulate.dmc to make a second p.df with nonPM accumulation rates
#which random accepts
#added a function two_p_lists which gets a second p.list with nonPM
#accumulation rates overwriting PM trial accumulation rates
#Currently this requires hand-coding each design

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications



#reps out matrix
rep_matrix <- function(theta_col, nrow) {
  
  matrix(rep(theta_col, times=nrow), nrow=length(theta_col))
  
}


transform.dmc <- function(par.df)
{
  #Regression using number of stimulus repetitions
  # attr(par.df,"cvs") should have the number of times 
  # that a PM stimulus has previously been repeated.
  #For single target PM this number is divided by 8 (to match the multi target condition)
  
  stimreps <- attr(par.df,"cvs")
  par.df$b <- par.df$B + par.df$A
  #!is.data.frame for likelihood
  #is.data.frame for random
  
  #LIKELIHOOD
  if (!is.data.frame(par.df)) {
    #convert the PM failure rate from the probit scale
    par.df$pmf <-
      apply(par.df$pmf, MARGIN = c(1, 2) , FUN = pnorm)
    
    #regression method for covariates
    #MVs is mean_v + the number of times PM has been repeated times by the slope.
    #Note the slope value will be specified as 0 for all non-PM trials
    MVs <- par.df$mean_v + stimreps * par.df$slope
    
    #Get the slope for response inhibition based on slope for mvs
    #Basically take the slope for the PM and times by par.df
    inh_slope = par.df$slope * par.df$inh_ex_rel
    inh_slope$W = inh_slope$P
    inh_slope$N = inh_slope$P
    inh_slope$P = 0
    INH <- par.df$inh + (stimreps * inh_slope)
    
    A <- par.df$A
    B <- par.df$B
    
    #Learning on the Bs is parameterized differently because we need B min
    #to be 0 (i.e., b can't be <A)
    #
    #  slope = (total diff in Bs / diff in stimreps)
    
    #get the difference in number of stimreps between first and last trial
    #for learning:
    # total_x_diff = max(stimreps) - min(stimreps)
    
    #Here I just hard coded as 7 to avoid calculating
    total_x_diff=7
    
    B$P <- ((par.df$finalB$P - par.df$B$P)/total_x_diff) * stimreps + par.df$B$P
    
  } else 
    
    #RANDOM
    
    {

    par.df$pmf <- pnorm(par.df$pmf)
    
    MVs <- rep_matrix(par.df$mean_v, length(stimreps)) + stimreps*rep_matrix(
      par.df$slope, length(stimreps))
    
    # browser()
    # 
    inh_slope = par.df$slope * par.df$inh_ex_rel
    inh_slope[1] = inh_slope[3]
    inh_slope[2] = inh_slope[3]
    inh_slope[3] = 0

    INH <- rep_matrix(par.df$inh, length(stimreps)) + par.df$inh_ex_rel[1]*stimreps*rep_matrix(
      inh_slope, length(stimreps))
    
    A <-  rep_matrix(par.df$A, length(stimreps))
    
    B_int <- rep_matrix(par.df$B, length(stimreps))
    B_fin <- rep_matrix(par.df$finalB, length(stimreps))
    
    B <- B_int
    B[3,] <- ((B_fin[3,] - B_int[3,])/7) * stimreps + B_int[3,]    
    
    
    
    
  }
  
# if (any(par.df$inh!=0)) browser()
  list(A = A,b = B + A, t0 = par.df$t0,
       mean_v = MVs,sd_v = par.df$sd_v,
       st0 = par.df$st0, N = par.df$N,pmf = par.df$pmf,
       inh=INH
  )
}


random.dmc <- function(n, p.df, model)
{
  #p.df2 is taking from simulate. On non-PM trials,
  #p.df2 should be the same as p.df
  #on PM trials, p.df2 will include the parameters for non-PM trials
  
  ##if else statement so that still works if there are 100% pm trigger failures
  if (p.df$pmf[1] < 1) {
    #determine number of PM trigger failures
    pmfs <- sum(rbinom(n, 1, p.df$pmf[1]))
    # pmfails <-  suppressWarnings(
    #   rlba.norm(
    #     pmfs,
    #     A = t(p.df$A[,1:2]),
    #     b = t(p.df$b[,1:2]),
    #     mean_v = t(p.df$mean_v[,1:2]),
    #     sd_v = p.df$sd_v[1:2],
    #     t0 = p.df$t0[1],
    #     st0 = p.df$st0[1],
    #     posdrift = attr(model, "posdrift")
    #   )
    # )
    browser()
    #trials where PM didn't fail -
    #simulate 1:3 accumulators from p.df
    nonfails <- rlba.norm(
      n - pmfs,
      A = p.df$A[1:p.df$N[1],],
      b = p.df$b[1:p.df$N[1],],
      mean_v = p.df$mean_v[1:p.df$N[1],] - p.df$inh[1:p.df$N[1],],
      sd_v = p.df$sd_v[1:p.df$N[1]],
      t0 = p.df$t0[1],
      st0 = p.df$st0[1],
      posdrift = attr(model, "posdrift")
    )
    
    # out <- rbind(pmfails, nonfails)
    out <- rbind(nonfails)
  } else {
    #if all PM failures, always use p.df2 and 
    #2 accumulators
    out <- rlba.norm(
      n,
      A = p.df$A[1:2,],
      b = p.df$b[1:2,],
      mean_v = p.df$mean_v[1:2,],
      sd_v = p.df$sd_v[1:2],
      t0 = p.df$t0[1],
      st0 = p.df$st0[1],
      posdrift = attr(model, "posdrift")
    )
  }
  out
}



# data=dm; min.like = 1e-10
likelihood.dmc <- function(p.vector, data, min.like = 1e-10)
  # Returns vector of likelihoods for each RT in data (in same order)
{
  
  p.list <- p.list.dmc(p.vector,model=attributes(data)$model,n1order=TRUE,
      cells=attributes(data)$cells,
      n1.index=attr(data,"n1.index"),
      cvs=data[,attr(attributes(data)$model,"cvs")]
  )
    
  likelihood <- numeric(dim(data)[1])
  
  #control block
  is2 <- p.list$N[, 1] == 2
  #not control block, not a PM trial, not a PM response
  is3OT <- !is2 & !data$R == "P" 
  #PM response
  is3PM <- !is2 & data$R == "P"

  
  #control block likelihood - L(two accumulators)
  likelihood[is2] <- n1PDFfixedt0.norm(
    dt = data$RT[is2] - p.list$t0[is2, 1],
    A = p.list$A[is2, 1:2],
    b = p.list$b[is2, 1:2],
    mean_v = p.list$mean_v[is2, 1:2],
    sd_v = p.list$sd_v[is2, 1:2],
    posdrift = attr(attr(data, "model"), "posdrift")
  )
  
  #not control block, not a PM trial, not a PM response likelihood
  #pmf * L(two accumulators) + (1-pmf) * L(three accumulators)
  likelihood[is3OT] <- (p.list$pmf[is3OT, 1]) * n1PDFfixedt0.norm(
    dt = data$RT[is3OT] - p.list$t0[is3OT, 1],
    A = p.list$A[is3OT, 1:2],
    b = p.list$b[is3OT, 1:2],
    mean_v = p.list$mean_v[is3OT, 1:2],
    sd_v = p.list$sd_v[is3OT, 1:2],
    posdrift = attr(attr(data, "model"), "posdrift")
  ) +
    (1 - p.list$pmf[is3OT, 1]) * n1PDFfixedt0.norm(
      dt = data$RT[is3OT] - p.list$t0[is3OT, 1],
      A = p.list$A[is3OT, 1:3],
      b = p.list$b[is3OT, 1:3],
      mean_v = p.list$mean_v[is3OT, 1:3] - p.list$inh[is3OT, 1:3],
      sd_v = p.list$sd_v[is3OT, 1:3],
      posdrift = attr(attr(data, "model"), "posdrift")
    )
  
  
  #PM response likelihood
  #(1-pmf) * L(3 accumulators)
  likelihood[is3PM] <- (1 - p.list$pmf[is3PM, 1]) * n1PDFfixedt0.norm(
    dt = data$RT[is3PM] - p.list$t0[is3PM, 1],
    A = p.list$A[is3PM, 1:3],
    b = p.list$b[is3PM, 1:3],
    mean_v = p.list$mean_v[is3PM, 1:3] - p.list$inh[is3PM, 1:3],
    sd_v = p.list$sd_v[is3PM, 1:3],
    posdrift = attr(attr(data, "model"), "posdrift")
  )
  
  pmax(likelihood, min.like)
}
# 
# p.df <- c()
# p.df$A <- c(0.5, 0.5, 0.5)
# p.df$b <- c(1, 1, 1)
# p.df$mean_v <- c(2,1, 2)
# p.df$sd_v <- c(1,1, 1)
# p.df$t0 <- c(0.1,0.1,0.1)
# p.df$st0 <- c(0,0,0)
# p.df$pmf <- c(0.1,0.1,0.1)
# p.df$N <- c(3,3,3)
# p.df$inh <- c(1,1,1)
# p.df <- as.data.frame(p.df)
# n=1e5
# model <- "model"
# attr(model, 'posdrift') = TRUE
# 
# sim = random.dmc(n=n, p.df, model)
# #test likelihoods
# names(sim) <- c("RT", "R")
# par(mfrow=c(1,3))
# dns <- plot.cell.density(sim,xlim=c(0,4),save.density=TRUE)
# print(pc)
# dt=dns$correct$x
# d <-  1 * n1PDFfixedt0.norm(
#     dt = dt - p.list$t0[is3PM, 1],
#     A = p.list$A[is3PM, 1:3],
#     b = p.list$b[is3PM, 1:3],
#     mean_v = p.list$mean_v[is3PM, 1:3] - p.list$inh[is3PM, 1:3],
#     sd_v = p.list$sd_v[is3PM, 1:3],
#     posdrift = attr(attr(data, "model"), "posdrift")
#   )
#   
#   
#   
#   n1Wald(dt,A=A,v=v,B=B,t0=t0)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")
# dt=dns$error$x
# d <- n1Wald(dt,A=A[2:1],v=v[2:1],B=B[2:1],t0=t0)
# plot(dns$error$x,dns$error$y,lty=2,type="l")
# lines(dns$error$x,d,col="red",lty=2)
# 
# # 
# # # Check go failure 
# # n=1e5
# # v=c(2,1); B=c(1,1); A=c(2,2); t0=1; gf=.2
# # sim <- rWaldRace(n=n,A=A,v=v,B=B,t0=t0,gf=gf)
# # par(mfrow=c(1,3))
# # dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# # pc <- integrate(n1Wald,lower=0,upper=Inf,A=A,v=v,B=B,t0=t0)$value
# # print(pc)
# # dt=dns$correct$x
# # d <- n1Wald(dt,A=A,v=v,B=B,gf=gf,t0=t0)
# # # red=black?
# # plot(dns$correct$x,dns$correct$y,type="l")
# # lines(dns$correct$x,d,col="red")
# # dt=dns$error$x
# # d <- n1Wald(dt,A=A[2:1],v=v[2:1],B=B[2:1],gf=gf,t0=t0)
# # plot(dns$error$x,dns$error$y,lty=2,type="l")
# # lines(dns$error$x,d,col="red",lty=2)
# 
# # # Check (effectively) single accumulator case 
# # n=1e5
# # v=c(2,-1); B=c(1,1); A=c(2,2); t0=1
# # sim <- rWaldRace(n=n,A=A,v=v,B=B,t0=t0)
# # par(mfrow=c(1,4))
# # dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# # pc <- integrate(n1Wald,lower=0,upper=Inf,A=A,v=v,B=B,t0=t0)$value
# # print(pc)
# # dt=dns$correct$x
# # d <- n1Wald(dt,A=A,v=v,B=B,t0=t0)
# # # red=black?
# # plot(dns$correct$x,dns$correct$y,type="l")
# # lines(dns$correct$x,d,col="red")
# # 
# # v=c(-1,2); B=c(1,1); A=c(2,2); t0=1
# # sim <- rWaldRace(n=n,A=A,v=v,B=B,t0=t0)
# # dns <- plot.cell.density(sim,C=1,xlim=c(0,4),save.density=TRUE)
# # pe <- integrate(n1Wald,lower=0,upper=Inf,A=A,v=v,B=B,t0=t0)$value
# # print(pe)
# # dt=dns$error$x
# # d <- n1Wald(dt,A=A[2:1],v=v[2:1],B=B[2:1],gf=gf,t0=t0)
# # # red=black?
# # plot(dns$error$x,dns$error$y,type="l")
# # lines(dns$error$x,d,col="red")
