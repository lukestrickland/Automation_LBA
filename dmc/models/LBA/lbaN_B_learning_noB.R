# Template setup for n-choice LBA, B=b-A parameterization
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# Allows for a variable number of accumulators by passing parameter for number
# of accumulators in p.df$N[1]
#
# This ncl (no cell likelihood) version uses p.list.dmc in the likelihood to 
# avoid looping over design cells.


# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

t_matrix <- function(theta_col, nrow) {
  
  t(matrix(rep(theta_col, times=nrow), nrow=length(theta_col)))
  
}

transform.dmc <- function(par.df)
{
  #Regression using number of stimulus repetitions
  # attr(par.df,"cvs") should have the number of times 
  # that a PM stimulus has previously been repeated.
  #For single target PM this number is divided by 8 (to match the multi target condition
  stimreps <- attr(par.df,"cvs")
  
  #Need to change the covariate format if calling from model.dmc to
  #dodge weird error
  if(is.data.frame(stimreps)) stimreps <- stimreps[,1]
  
  #Make finalB parameters meaningless if B is fixed
  
  

  #!is.data.frame for likelihood
  #is.data.frame for random
  #LIKELIHOOD
  if (!is.data.frame(par.df)) {
    #regression method for covariates
    #MVs is mean_v + the number of times PM has been repeated times by the slope.
    #Note the slope value will be specified as 0 for all non-PM trials
    MVs <- par.df$mean_v + stimreps * par.df$slope

    #Get the slope for response inhibition based on slope for mvs
    #Basically take the slope for the PM and times by inh_ex_rel
    inh_slope = par.df$slope * par.df$inh_ex_rel
    #Structure the data frame appropriately:
    #extra inhibition for w/nonword accumulators and no inhibition of PM
    inh_slope["W"] = inh_slope["P"]
    inh_slope["N"] = inh_slope["P"]
    inh_slope["P"] = 0
    
    #calculate level of inhibition given current stimreps
    INH <- par.df$inh + (stimreps * inh_slope)

    A <- par.df$A
    B <- par.df$B
    
  } else {
    
    #RANDOM

    #
    MVs <- t_matrix(par.df$mean_v, length(stimreps)) + stimreps*t_matrix(
      par.df$slope, length(stimreps))

    inh_slope = par.df$slope * par.df$inh_ex_rel
    inh_slope[1] = inh_slope[3]
    inh_slope[2] = inh_slope[3]
    inh_slope[3] = 0
    
    INH <- t_matrix(par.df$inh, length(stimreps)) + stimreps*t_matrix(
      inh_slope, length(stimreps))

    A <-  t_matrix(par.df$A, length(stimreps))
    B <- t_matrix(par.df$B, length(stimreps))
    
  }

  list(A = A,b = B + A, t0 = par.df$t0,
       mean_v = MVs,sd_v = par.df$sd_v,
       st0 = par.df$st0, N = par.df$N,
       inh=INH
  )
}

random.dmc<- function(n,p.df,model)
{
  rlba.norm(n,
            
    A = t(p.df$A[,1:p.df$N[1]]),
    b = t(p.df$b[,1:p.df$N[1]]),
    mean_v = t(p.df$mean_v[,1:p.df$N[1]] - p.df$inh[,1:p.df$N[1]]),
    sd_v=p.df$sd_v[1:p.df$N[1]],
    t0=p.df$t0[1], 
    st0=p.df$st0[1],
    posdrift = attr(model,"posdrift"))
}


likelihood.dmc <- function(p.vector,data,min.like=1e-10) 
# Returns vector of likelihoods for each RT in data (in same order)
{
   
  p.list <- p.list.dmc(p.vector,model=attributes(data)$model,n1order=TRUE,
      cells=attributes(data)$cells,
      n1.index=attr(data,"n1.index"),
      cvs=data[,attr(attributes(data)$model,"cvs")]
  )
  
  
  
  likelihood <- numeric(dim(data)[1])
  is2 <- p.list$N[,1]==2
  
  
  likelihood[is2] <- n1PDFfixedt0.norm(
    dt=data$RT[is2]-p.list$t0[is2,1],
    A=p.list$A[is2,1:2],
    b=p.list$b[is2,1:2],
    mean_v=p.list$mean_v[is2,1:2],
    sd_v=p.list$sd_v[is2,1:2],
    posdrift=attr(attr(data,"model"),"posdrift"))
  
  likelihood[!is2] <- n1PDFfixedt0.norm(
    dt=data$RT[!is2]-p.list$t0[!is2,1],
    A=p.list$A[!is2,1:3],
    b=p.list$b[!is2,1:3],
    mean_v= p.list$mean_v[!is2, 1:3] - p.list$inh[!is2, 1:3],
    sd_v=p.list$sd_v[!is2,1:3],
    posdrift=attr(attr(data,"model"),"posdrift"))

  pmax(likelihood,min.like)
}

#Test that transform gives expected inputs to random
#and likelihood

#test that rlba can handle the inputs the way they are provided

#test that likelihood can handle the inputs the way they are provided

