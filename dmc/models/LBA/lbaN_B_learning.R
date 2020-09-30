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
  if(is.data.frame(stimreps)) stimreps <- stimreps[,1]
  
  #!is.data.frame for likelihood
  #is.data.frame for random
  
  if (!is.data.frame(par.df)) {
  #LIKELIHOOD  
    #regression method for covariates
    #Note the slope value will be specified as 0 for all non-PM trials
    MVs <- par.df$mean_v + stimreps * par.df$slope

    #Get the slope for response inhibition based on slope for mvs
    #linear function with inh_ex_rel determining how strong the relationship is 
    #between excitation learning and inhibition learning
    inh_slope = par.df$slope * par.df$inh_ex_rel
    
    #Structure the inhibition appropriately:
    #inhibition for w/nonword accumulators and no inhibition of PM
    inh_slope["W"] = inh_slope["P"]
    inh_slope["N"] = inh_slope["P"]
    inh_slope["P"] = 0
    
    #calculate level of inhibition given current stimreps (i.e. number
    #of times stimulus has been presented)
    INH <- par.df$inh + (stimreps * inh_slope)

    A <- par.df$A
    
    #Learning on the PM thresholds:
    #parameterized differently because we need B min
    #to be 0 (i.e., b can't be <A)
    #
    #  m  = (total diff in Bs / diff in stimreps)
    # minB = PM threshold on first presentation of items (stimreps=0)
    # PM threshold = m * stimreps + minB

    #Here I just hard coded total diff in stimreps as 7 to avoid calculating
    total_x_diff=7
   
    B <- par.df$B
    B[,3] <- ((par.df$finalB[,3] - par.df$B[,3])/total_x_diff) * stimreps + par.df$B[,3]
    
      # like.ps <<- list(A = A,b = B + A, t0 = par.df$t0,
      #  mean_v = MVs,sd_v = par.df$sd_v,
      #  st0 = par.df$st0, N = par.df$N,
      #  inh=INH  )
    
  } else {
    
  #RANDOM
    MVs <- t_matrix(par.df$mean_v, length(stimreps)) + stimreps*t_matrix(
      par.df$slope, length(stimreps))

    inh_slope = par.df$slope * par.df$inh_ex_rel
    inh_slope[1] = inh_slope[3]
    inh_slope[2] = inh_slope[3]
    inh_slope[3] = 0
    
    INH <- t_matrix(par.df$inh, length(stimreps)) + stimreps*t_matrix(
      inh_slope, length(stimreps))

    A <-  t_matrix(par.df$A, length(stimreps))
    
    B_int <- t_matrix(par.df$B, length(stimreps))
    B_fin <- t_matrix(par.df$finalB, length(stimreps))
    B <- B_int
    B[,3] <- ((B_fin[,3] - B_int[,3])/7) * stimreps + B_int[,3]

    # random.ps <<- list(A = A,b = B + A, t0 = par.df$t0,
    #    mean_v = MVs,sd_v = par.df$sd_v,
    #    st0 = par.df$st0, N = par.df$N,
    #    inh=INH)
  }
  # if (any (par.df$inh!=0)) browser()
  list(A = A,b = B + A, t0 = par.df$t0,
       mean_v = MVs,sd_v = par.df$sd_v,
       st0 = par.df$st0, N = par.df$N,
       inh=INH
  )
}

random.dmc<- function(n,p.df,model)
{
  # browser()
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
