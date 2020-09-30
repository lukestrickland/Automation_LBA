#LBA model with potential learning 
#effects on PM rates, inhibition of OT and thresholds


#Function for use in simulation
#that converts all parameter values
#to long matrices. Then I use
#covariates to modify each matrix
#appropriately

t_matrix <- function(theta_col, nrow) {
  
  t(matrix(rep(theta_col, times=nrow), nrow=length(theta_col)))
  
}

transform.dmc <- function(par.df)
{
  #Regression using number of stimulus repetitions
  # attr(par.df,"cvs") should have the number of times 
  # that a PM stimulus has previously been repeated.
  #For single target PM this number is divided by 8 (to match the multi target condition)
  stimreps <- attr(par.df,"cvs")
  
  #Need to change the covariate format if calling from model.dmc to
  #dodge weird error
  #Also need to change alpha if calling from model.dmc
  ALPH <- par.df$ALPH[1]
  if(is.data.frame(stimreps)) {
    stimreps <- stimreps[,1]
    ALPH <- par.df$ALPH
  }
  
  #!is.data.frame for likelihood
  #is.data.frame for random
  
  if (!is.data.frame(par.df)) {
  #LIKELIHOOD  
    #regression method for covariates
    #Note the slope value will be specified as 0 for all non-PM trials
    MVs <- par.df$mean_v - par.df$LRN*exp(-ALPH*stimreps)[,1]

    #Get the slope for response inhibition based on slope for mvs
    #linear function with inh_ex_rel determining how strong the relationship is 
    #between excitation learning and inhibition learning
    inh_LRN = par.df$LRN * par.df$inh_ex_rel
    
    #Structure the inhibition appropriately:
    #inhibition for w/nonword accumulators and no inhibition of PM
    inh_LRN["W"] = inh_LRN["P"]
    inh_LRN["N"] = inh_LRN["P"]
    inh_LRN["P"] = 0
    
    #calculate level of inhibition given current stimreps (i.e. number
    #of times stimulus has been presented)
    INH <- par.df$inh - inh_LRN * exp(-ALPH*stimreps)[,1]
      
  } else {
    
  #RANDOM
    MVs <- t_matrix(par.df$mean_v, length(stimreps)) - t_matrix(
      par.df$LRN, length(stimreps)) * exp(-ALPH*stimreps)

    inh_LRN = par.df$LRN * par.df$inh_ex_rel
    inh_LRN[1] = inh_LRN[3]
    inh_LRN[2] = inh_LRN[3]
    inh_LRN[3] = 0
    
    INH <- t_matrix(par.df$inh, length(stimreps)) - t_matrix(
      inh_LRN, length(stimreps)) * exp(-ALPH*stimreps)


  }
  # if (any (par.df$inh!=0)) browser()
  list(A = par.df$A,b = par.df$B + par.df$A, t0 = par.df$t0,
       mean_v = MVs,sd_v = par.df$sd_v,
       st0 = par.df$st0, N = par.df$N,
       inh=INH
  )
}

random.dmc<- function(n,p.df,model)
{
  rlba.norm(n,
    A = p.df$A[1:p.df$N[1]],
    b = p.df$b[1:p.df$N[1]],
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
