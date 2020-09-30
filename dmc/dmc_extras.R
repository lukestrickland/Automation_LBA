#Miscellaneous functions that have accumulated over time
#Many from early days of working on this stuff, need refactoring
#test submodule updating
#test 2
#check
require("gridExtra")
require("lme4")
# require("plyr")
require("dplyr")
require("data.table")

singlerep.mean.sd <- function(rsamples, fun) {
  inference <- list()
  for (i in 1:length(rsamples)) {
    thetas <- rsamples[[i]]$theta
    inference [[i]] <- fun (thetas)
  }
  inf2 <- unlist(inference)
  dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
  dim(inf2) <- dim3
  test<- apply(inf2, 4, function(x) c(mean(x),sd(x)))
  test<- apply(test,1,mean); names(test) <- c("M", "SD")
  paste(round(test[1],3), " (", round(test[2],3), ")", sep="")
}



get.95.50 <- function (rsamples, p.vector, fun){ 
  test <- c(h.check.function.recovery.dmc(
  rsamples, p.vector, fun), h.check.function.recovery.dmc(
  rsamples, p.vector, fun, CI=c(0.25,0.75)))
  paste(round(test[1],2), "% / ", round(test[2],2), "%", sep="")
  }

h.check.function.recovery.dmc <- function(rsamples, p.vector, fun, CI= c(0.025, 0.975)) {
#Get the posterior mean and credible intervals for each replicate
  inference <- list()
  for (i in 1:length(rsamples)) {
    thetas <- rsamples[[i]]$theta
    inference [[i]] <- fun (thetas)
  }
  inf2 <- unlist(inference)
  dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
  dim(inf2) <- dim3
  test<- apply(inf2, 4, function(x) c(mean(x), quantile(x, probs=CI)))
#Hacky way to reshape p.vector same way to get functions in same way
#I should probably make fun work with sim.p.vector rather than thetas for future work,
#simpler, but already too far down the rabbit hole with this way for this paper
  pnames <- names(p.vector)
  dim(p.vector) <- c(1,length(p.vector),1); colnames(p.vector) <- pnames
  actual <- fun(p.vector)
#percentage of time CI contains actual  
  sum(apply(test,2, function(x) {x[2]<actual & x[3] >actual})) / length(test[1,]) *100
}

get.cors <- function(thetas, data) {
  
  CORS <- apply(thetas, c(1,2,3), function(x) cor(x, data))

  RAV <- apply(CORS, c(2), postRav, n=length(data), kappa=1, spacing=.01)
  post_medians <- apply(RAV, c(2), postRav.mean)
  post_LCI <- apply(RAV, c(2), function(x) postRav.ci(x)[1])
  post_HCI <- apply(RAV, c(2), function(x) postRav.ci(x)[2])
  out <- list(post_medians, post_LCI, post_HCI)
  names(out) <- c("medians", "LCI", "HCI")
  out
}

get.subj.effects.m <- function (PPS, fun, names) {
  EFFECTS <- lapply(PPS, get.subj.effects, fun=fun)
  for (i in 1:length(names)) EFFECTS[[i]]$model <- names[i]
  do.call(rbind, EFFECTS)
}

get.subj.effects <- function(PP, fun) {
  for (i in 1:length(PP)){
    effects<-get.subj.pp.MCI (PP[[i]], fun)
    cat(i)
    if (i ==1) out <- effects else out <- rbind(out,effects)
  }
  
  ENS <- rownames(out)
  OUT <- data.frame(out)
  OUT$effect <- ENS
  colnames(OUT) <- c("mean", "lower", "upper", "data", "effect")
  OUT
}


subj.meanthetas <- function (samples){
    samps <- lapply(samples, function(x) x["theta"])
    ## thetas into big array for apply
    samps2<- unlist(samps)
    dim3 <- c(dim(samps[[1]]$theta), length(samps2)/prod(dim(samps[[1]]$theta)))
    dim(samps2) <- dim3
    samps3<- apply(samps2, c(4,2), mean)
    ## back to a theta list after applied
    colnames(samps3)<- colnames(samps[[1]]$theta)
    df <- cbind(names(samples), data.frame(samps3))
    names(df)[1] <- "s"
    df
}

get.subj.pp.MCI <- function(sim, fun) {
 
  data <- attr(sim, "data")
  nreps=max(sim$reps)
  data.effects<- fun(data)
  noutput <- length(data.effects)

  sim.effects <- matrix(NA, nrow= nreps, ncol=noutput+1)
  sim.effects[,noutput+1] <- 1:nreps

  colnames(sim.effects) <- c(names(data.effects), "n.rep")

  for (j in 1:nreps) {

    currentsim.effects <- sim[sim$reps==j,]
    sim.effects[j,1:noutput] <- fun(currentsim.effects)

  }
  out <- apply(sim.effects, 2,function(x) c(mean(x, na.rm=T), quantile(x, probs=c(0.025,0.975), na.rm=T))) 
  cbind(t(out[,(!colnames(out) %in% "n.rep")]), data.effects)
  }


convert.magic <- function(obj,types){
    for (i in 1:length(obj)){
        FUN <- switch(types[i],character = as.character, 
                                   numeric = as.numeric, 
                                   factor = as.factor)
        obj[,i] <- FUN(obj[,i])
    }
    obj
}


tabtoAPA <- function (tab){
  for (i in seq (1, length(tab), 2))  {
   newtab<- paste(tab[,i], " (", tab[,i+1], ")", sep="")
  if (i == 1) { new.tab <- newtab } else {new.tab <- cbind (new.tab, newtab)}
   i = i+1
  }
   rownames(new.tab)<- rownames (tab)
   new.tab
  
}

get.theta.array <- function(hsamples){
  samps <- lapply(hsamples, function(x) x["theta"])
  samps2<- unlist(samps)
  dim3 <- c(dim(samps[[1]]$theta), length(samps2)/prod(dim(samps[[1]]$theta)))
  dim(samps2) <- dim3
  colnames(samps2) <- colnames(hsamples[[1]]$theta)
  samps2
}


GET.fitgglist.dmc <- function (
  PP, factors=NA, noR = FALSE,
  quantiles.to.get = c(0.1, 0.5, 0.9), CI= c(0.025, 0.975),
  acc.fun=function(x){as.numeric(x$S)==as.numeric(x$R)},
  correct.only=FALSE,error.only=FALSE
  
) {
 
  sim <- do.call(rbind, PP)
  # Do the same for the data
  data <- lapply(PP, function(x) attr(x, "data"))
  data <- do.call(rbind, data)
  get.fitgglist.dmc (sim,data, factors=factors, noR=noR, quantiles.to.get=quantiles.to.get,
                     CI=CI, acc.fun=acc.fun, correct.only=correct.only, error.only=
                       error.only)
  
}


get.ns.dmc<- function(samples) {
model <- attributes(samples$data)$model
facs <- names(attr(model,"factors"))
table(samples$data[,facs],dnn=facs)}

ggplot.recov <- function(ggdf, ncol=10) {
  ggplot(ggdf, aes(reps, M))+
    geom_point(size=0.5)+
    geom_hline(aes(yintercept=data), linetype=1, size=1.5)+ylab("")+ 
    geom_ribbon(data=ggdf,aes(ymin=LCI,ymax=HCI), alpha=0.3)+ facet_wrap(~param, ncol=ncol) 
}


get.ggdf.recov <- function(post_summaries, msds, grepchr="B") {
  tmp <- list()
  j=0
  for (i in grep(grepchr, colnames(post_summaries[[1]]))){
    j = j + 1
    tmp [[j]] <- data.frame(cbind(post_summaries[[1]][,i], 
                                   post_summaries[[2]][,i], post_summaries[[3]][,i]))
    name <- rownames(msds)[i]
    colnames(tmp[[j]]) <- c("M", "LCI", "HCI")
    tmp[[j]]$reps <- 1:100
    tmp[[j]]$param <- name
    tmp[[j]]$data <- msds$M[i]
  }
  do.call(rbind, tmp)
}


get.participant.median.CIs <- function(hsamples) {
  samps <- lapply(hsamples, function(x) x["theta"])
  samps2<- unlist(samps)
  dim3 <- c(dim(samps[[1]]$theta), length(samps2)/prod(dim(samps[[1]]$theta)))
  dim(samps2) <- dim3
  post_medians <- apply(samps2, c(4,2), mean)
  post_LCI <- apply(samps2, c(4,2), quantile, prob=0.025)
  post_HCI <- apply(samps2, c(4,2), quantile, prob=0.975)
  colnames(post_medians) <- colnames(hsamples[[1]]$theta)
  colnames(post_LCI) <- colnames(hsamples[[1]]$theta)
  colnames(post_HCI) <- colnames(hsamples[[1]]$theta)
  out <- list(post_medians, post_LCI, post_HCI)
  names(out) <- c("medians", "LCI", "HCI")
  out
}

# A few functions for posterior predicctive p values and z scores.
minp <- function (effect) min(ecdf(effect)(0), 1-ecdf(effect)(0))

zandp <- function(samples, fun){
    effect<- group.inference.dist(samples, fun)
    Z <- mean(effect)/sd(effect)
    p <- minp(effect)
    paste(round(Z,2), " (", round(p,3), ")", sep="")
}

mean.sd <- function(samples, fun){
    effect<- group.inference.dist(samples, fun)
    M <- mean(effect)
    SD <- sd(effect)
    paste(round(M,3), " (", round(SD,3), ")", sep="")
}

  Z.p.acrossexp <- function(samples1,samples2, fun){
    effect1<- group.inference.dist(samples1, fun)
    effect2 <- group.inference.dist(samples2, fun)
    effect<- effect1 - effect2
    Z <- mean(effect)/sd(effect)
    p <- minp(effect)
    paste(round(Z,2), " (", round(p,3), ")", sep="")
  }

##accepts a function and does it to the thetas for each subject then averages after
group.inference.dist <- function (hsamples, fun) {
  inference <- list()
  for (i in 1:length(hsamples)) {
    thetas <- hsamples[[i]]$theta
    inference [[i]] <- fun (thetas)
  }
  inf2 <- unlist(inference)
  dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
  dim(inf2) <- dim3
  apply(inf2, c(1,2,3), mean)
}

## Below funciton averages parameters across conditions by grepping out
# from av.posts. So if you set av.posts to match anything containing mean_v,
# it would average all rates together and replace all values with the avg before
# simming. We use it to parse the effects of rates/thresholds on the manifests
# in terms of block and cond.

# samples=samples[[3]]
# probs=c(1:99)/100;random=TRUE
# bw="nrd0";report=10;save.simulation=TRUE;factors=NA
# av.posts<-av.posts.threscond
avps.post.predict.dmc = function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
                                    bw="nrd0",report=10,save.simulation=TRUE,factors=NA, av.posts=c())
  # make list of posterior preditive density, quantiles and response p(robability)
{


  get.dqp <- function(sim,facs,probs,n.post=NA) {

    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...) {
      out <- quantile(x,probs=probs,na.rm=na.rm,type=type,names=FALSE,...)
      names(out) <- probs*100
      if (all(is.na(out))) NULL else out
    }

    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)
    #     qs <- apply(qs,1:length(dim(qs)),function(x){
    #       if ( is.null(x[[1]]) || all(is.na(x[[1]]))) NULL else x[[1]]})

    # get probabilities given not na
    simOK <- sim[!is.na(sim$RT),]
    p <- tapply(simOK$RT,simOK[,c(facs,"R")],length)
    p[is.na(p)] <- 0 # In case some cells are empty
    np <- rep(apply(p,1:length(facs),sum),times=length(levels(simOK$R)))
    p <- p/np
    # get p NA
    pNA <- tapply(sim$RT,sim[,c(facs,"R")],length)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    npNA <- rep(apply(pNA,1:length(facs),sum),times=length(levels(sim$R)))
    pNA <- tapply(is.na(sim$RT),sim[,c(facs,"R")],sum)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    pNA <- pNA/npNA

    # For a simulation get proability replicates
    if (!is.na(n.post)) {
      repfac <- rep(1:n.post,each=sum(ns))
      repfac <- repfac[!is.na(sim$RT)]
      ps <- tapply(simOK$RT,cbind(simOK[,c(facs,"R")],rep=repfac),length)
      ps[is.na(ps)] <- 0 # In case some cells are empty
      ps <- n.post*ps/np
      # and NA replicates
      pNAs <- tapply(is.na(sim$RT),cbind(sim[,c(facs,"R")],rep=rep(1:n.post,each=sum(ns))),sum)
      pNAs[is.na(pNAs)] <- 0 # In case some cells are empty
      pNAs <- n.post*pNAs/npNA
    } else {
      ps=NULL
      pNAs=NULL
    }

    # cell names
    cell.names <- dimnames(qs)[[1]]
    n.cell <- length(facs)+1
    for ( i in 2:n.cell )
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    cell.names <- as.vector(cell.names)
    # Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R")],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01,na.rm=TRUE) & x<=quantile(x,.99,na.rm=TRUE)]
        if (length(x[!is.na(x)])<2) NULL else density(x[!is.na(x)],bw=bw)
      }
    })
    for (i in 1:length(p)) if ( is.finite(p[i]) && !(p[i]==0) )
    {
      if (!is.null(qs[i][[1]])) {
        names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
        attr(qs[i][[1]],"cell.name") <- cell.names[i]
      }
      if (!is.null(dens[i][[1]]) ) {
        dens[i][[1]]$y <- dens[i][[1]]$y*p[i]
        attr(dens[i][[1]],"cell.name") <- cell.names[i]
      }
    }
    dnd <- dimnames(dens)
    dnq <- dimnames(qs)
    dens <- apply(dens,1:length(facs),function(x){x})
    qs <- apply(qs,1:length(facs),function(x){x})
    if ( is.null(dim(dens)) ) {
      dens <- array(dens,dim=c(length(dens)))
      dimnames(dens) <- dnd[-length(dnd)]
      qs <- array(qs,dim=c(length(qs)))
      dimnames(qs) <- dnq[-length(dnq)]
    }
    list(pdf=dens,cdf=qs,p=p,ps=ps,pNA=pNA,pNAs=pNAs)
  }

  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  
  cvs <- samples$data[,attr(model,"cvs"), drop=FALSE]
  attr(cvs,"row.facs") <- apply(apply(
  samples$data[,facs,drop=FALSE],2,as.character),1,paste,collapse=".")
  
  if (any(is.na(factors))) factors <- facs
  if (!all(factors %in% facs))
    stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]


  cat("Below is how I'm averaging (each row is averaged). If this is wrong, adjust your
      av.posts to grep correctly.")
  for (q in 1:length(av.posts)) print(colnames(posts[, grep(av.posts[q], colnames(posts))]))
  ###tweak to average av.posts
  q=1

  if(length(av.posts)!= 0) {
    ### loop through all the averaged posts
    for (q in 1:length(av.posts)) {

      num.params <- dim(posts[, grep(av.posts[q], colnames(posts))])[2]
      average.params <- rowMeans(posts[, grep(av.posts[q], colnames(posts))])
      posts[, grep(av.posts[q], colnames(posts))] <- matrix(average.params,nrow=length(average.params),ncol=num.params,byrow=F)

    }
  }

  ########
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns)) 
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data) %in% c("RT",names(cvs),"R2")]
  } else {
    # Assumes last two are SSD and RT! FIX ME. EG WONT WORK IF THERE ARE CVS
    if ( is.null(facs) ) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity)) 
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data) %in% c("RT","SSD")]
    # leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(posts[i,],model,ns,SSD=SSD,cvs=cvs)
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="")
    sim$R[sim$R2=="DK"] <- "DK"
    sim$R <- factor(sim$R)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="")
    samples$data$R[samples$data$R2=="DK"] <- "DK"
    samples$data$R <- factor(samples$data$R)
  }
  if (save.simulation) {
    sim <- cbind(reps=rep(1:n.post,each=dim(samples$data)[1]),sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,factors,probs,n.post)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    c(sim.dqp,dat.dqp)
  }
}


#lapplys the above function on everybody
avps.h.post.predict.dmc<- function(samples,n.post=100,probs=c(1:99)/100,
                                          bw="nrd0",
                                     save.simulation=FALSE, av.posts=c())
  # apply lost.predict to each subject
{
  lapply(samples,avps.post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
         save.simulation=save.simulation, av.posts=av.posts)
}

# PPs = E1PP
# fun = block.effects.E1A4
# lower=.025
# upper=.975
get.effects.dmc <- function (PPs, fun = function (x) {mean (x)}, lower=.025, upper=.975) {

  simdata<- do.call(rbind, PPs)
  data <- lapply(PPs, function(x) attr(x, "data"))
  data <- do.call(rbind, data)
  nreps=max(PPs[[1]]$reps)


  data.effects<- fun(data)
  noutput <- length(data.effects)

  sim.effects <- matrix(NA, nrow= nreps, ncol=noutput+1)
  sim.effects[,noutput+1] <- 1:nreps

  colnames(sim.effects) <- c(names(data.effects), "n.rep")
  ######

  ##calculate effects separately for each rep
  for (j in 1:nreps) {

    currentsim.effects <- simdata[simdata$reps==j,]
    sim.effects[j,1:noutput] <- fun(currentsim.effects)

  }
  ##Get a ggplot df with posterior mean, lower, and upper.
  effects.ggdf <-  t(apply(sim.effects, c(2), function(x) c(mean(x),
                                                quantile(x, probs=c(lower,upper)))))
  effects.ggdf <- data.frame(effects.ggdf)
  effects.ggdf <- effects.ggdf[(!rownames(effects.ggdf) %in% "n.rep"),]
  colnames(effects.ggdf) <- c("mean", "lower", "upper")
  contrast <- rownames(effects.ggdf)
  effects.ggdf$data<-as.vector(data.effects)
  attr(effects.ggdf, "post.effects.samples") <- sim.effects
  effects.ggdf
}

ss.get.effects.dmc <- function (PPs, fun = function (x) {mean (x)}, lower=.025, upper=.975) {

  simdata<- PPs
  data <- attr(PPs, "data")
  nreps=max(PPs$reps)


  data.effects<- fun(data)
  noutput <- length(data.effects)

  sim.effects <- matrix(NA, nrow= nreps, ncol=noutput+1)
  sim.effects[,noutput+1] <- 1:nreps

  colnames(sim.effects) <- c(names(data.effects), "n.rep")
  ######

  ##calculate effects separately for each rep
  for (j in 1:nreps) {

    currentsim.effects <- simdata[simdata$reps==j,]
    sim.effects[j,1:noutput] <- fun(currentsim.effects)

  }

  ##Get a ggplot df with posterior mean, lower, and upper.
  effects.ggdf <-  t(apply(sim.effects, c(2), function(x) c(mean(x, na.rm=T),
                                                quantile(x, probs=c(lower,upper), na.rm=T))))
  effects.ggdf <- data.frame(effects.ggdf)
  effects.ggdf <- effects.ggdf[(!rownames(effects.ggdf) %in% "n.rep"),]
  colnames(effects.ggdf) <- c("mean", "lower", "upper")
  contrast <- rownames(effects.ggdf)
  effects.ggdf$data<-as.vector(data.effects)
  attr(effects.ggdf, "post.effects.samples") <- sim.effects
  effects.ggdf
}


#The below function picks certain parameters from a list called pickps_set,
# and replaces them with pickps_other, before performing posterior prediciton.
# We use it to turn control mechanisms off in the model. To turn off proactive
# control, we set the ongoing task thresholds equal in the PM block to the control
#threshols. To turn off reactive, we set the ongiong rates on PM trials (PM block)
# to the ongoing rates on non-PM trials (PM block)

pickps.post.predict.dmc = function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
                                   bw="nrd0",report=10,save.simulation=TRUE,factors=NA, pickps_others, pickps_set,
                                   special_model=NA)
  # make list of posterior preditive density, quantiles and response p(robability)
{
  
  
  get.dqp <- function(sim,facs,probs,n.post=NA) {
    
    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...) {
      out <- quantile(x,probs=probs,na.rm=na.rm,type=type,names=FALSE,...)
      names(out) <- probs*100
      if (all(is.na(out))) NULL else out
    }
    
    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)
    #     qs <- apply(qs,1:length(dim(qs)),function(x){
    #       if ( is.null(x[[1]]) || all(is.na(x[[1]]))) NULL else x[[1]]})
    
    # get probabilities given not na
    simOK <- sim[!is.na(sim$RT),]
    p <- tapply(simOK$RT,simOK[,c(facs,"R")],length)
    p[is.na(p)] <- 0 # In case some cells are empty
    np <- rep(apply(p,1:length(facs),sum),times=length(levels(simOK$R)))
    p <- p/np
    # get p NA
    pNA <- tapply(sim$RT,sim[,c(facs,"R")],length)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    npNA <- rep(apply(pNA,1:length(facs),sum),times=length(levels(sim$R)))
    pNA <- tapply(is.na(sim$RT),sim[,c(facs,"R")],sum)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    pNA <- pNA/npNA
    
    # For a simulation get proability replicates
    if (!is.na(n.post)) {
      repfac <- rep(1:n.post,each=sum(ns))
      repfac <- repfac[!is.na(sim$RT)]
      ps <- tapply(simOK$RT,cbind(simOK[,c(facs,"R")],rep=repfac),length)
      ps[is.na(ps)] <- 0 # In case some cells are empty
      ps <- n.post*ps/np
      # and NA replicates
      pNAs <- tapply(is.na(sim$RT),cbind(sim[,c(facs,"R")],rep=rep(1:n.post,each=sum(ns))),sum)
      pNAs[is.na(pNAs)] <- 0 # In case some cells are empty
      pNAs <- n.post*pNAs/npNA
    } else {
      ps=NULL
      pNAs=NULL
    }
    
    # cell names
    cell.names <- dimnames(qs)[[1]]
    n.cell <- length(facs)+1
    for ( i in 2:n.cell )
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    cell.names <- as.vector(cell.names)
    # Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R")],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01,na.rm=TRUE) & x<=quantile(x,.99,na.rm=TRUE)]
        if (length(x[!is.na(x)])<2) NULL else density(x[!is.na(x)],bw=bw)
      }
    })
    for (i in 1:length(p)) if ( is.finite(p[i]) && !(p[i]==0) )
    {
      if (!is.null(qs[i][[1]])) {
        names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
        attr(qs[i][[1]],"cell.name") <- cell.names[i]
      }
      if (!is.null(dens[i][[1]]) ) {
        dens[i][[1]]$y <- dens[i][[1]]$y*p[i]
        attr(dens[i][[1]],"cell.name") <- cell.names[i]
      }
    }
    dnd <- dimnames(dens)
    dnq <- dimnames(qs)
    dens <- apply(dens,1:length(facs),function(x){x})
    qs <- apply(qs,1:length(facs),function(x){x})
    if ( is.null(dim(dens)) ) {
      dens <- array(dens,dim=c(length(dens)))
      dimnames(dens) <- dnd[-length(dnd)]
      qs <- array(qs,dim=c(length(qs)))
      dimnames(qs) <- dnq[-length(dnq)]
    }
    list(pdf=dens,cdf=qs,p=p,ps=ps,pNA=pNA,pNAs=pNAs)
  }
  
  if (any(is.na(special_model))) model <- attributes(samples$data)$model else model <-
      special_model
  facs <- names(attr(model,"factors"))
  if (any(is.na(factors))) factors <- facs
  if (!all(factors %in% facs))
    stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]
  #more robust
  ###Replace some parameter vlaues with others.
  
  #if to handle the case where only 1 post value
  if(length(pickps_others==1)) {
    posts[,pickps_others] <-  posts[,pickps_set] 
  } else {
    posts[,colnames(posts) %in% pickps_others][,pickps_others] <- 
    posts[,colnames(posts) %in% pickps_set][,pickps_set] 
    
  }

  
  ########
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns))
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data)=="RT"]
  } else {
    # Assumes last two are SSD and RT! FIX ME.
    if (is.null(facs)) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity))
    leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(posts[i,],model,ns,SSD=SSD)
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="")
    sim$R[sim$R2=="DK"] <- "DK"
    sim$R <- factor(sim$R)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="")
    samples$data$R[samples$data$R2=="DK"] <- "DK"
    samples$data$R <- factor(samples$data$R)
  }
  if (save.simulation) {
    sim <- cbind(reps=rep(1:n.post,each=dim(samples$data)[1]),sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,factors,probs,n.post)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    c(sim.dqp,dat.dqp)
  }
}

#lapply the above to the whole samples object.
pickps.h.post.predict.dmc<- function(samples,n.post=100,probs=c(1:99)/100,
                                   bw="nrd0",
                                   save.simulation=FALSE, pickps_set, pickps_others,
                                   special_model=NA)
  # apply lost.predict to each subject
{
  lapply(samples,pickps.post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
         save.simulation=save.simulation, pickps_set=pickps_set, pickps_others=pickps_others,
         special_model=special_model)
}



get.hdata.dmc <- function(hsamples) {
  list.wind <-
    lapply(seq_along(hsamples), function(samples, n, i)
      cbind(n[[i]], samples[[i]]$data),
      samples = hsamples, n = names(hsamples))
  out <- do.call(rbind, list.wind)
  names(out)[1] <- "s"
  out
}

fixedeffects.meanthetas <- function (samples) {
  
#handle different nmcs  
  nmcs<- sapply(samples, function(x) x$nmc)
  nmc <- min(nmcs)
#Different numbers of nmc for each participant... use the min number and then
  #for participants wtih more randomly sample out that many
  for (i in 1:length(samples)) if (nmcs[i] > nmc) samples[[i]]$theta <- 
    samples[[i]]$theta[,,sample(1:dim(samples[[i]]$theta)[3], nmc)]
####
  
    
  samps <- lapply(samples, function(x)
    x["theta"])
  ##
  ## thetas into big array for apply
  samps2 <- unlist(samps)
  dim3 <-
    c(dim(samps[[1]]$theta), length(samps2) / prod(dim(samps[[1]]$theta)))
  dim(samps2) <- dim3
  samps3 <- apply(samps2, c(1, 2, 3), mean)
  ## back to a theta list after applied
  colnames(samps3) <- colnames(samps[[1]]$theta)
  samps5 <- list(samps3)
  attributes(samps5) <- attributes(samps[[1]])
  samps5
}

get.pinf.subjects <- function(funs=list(mean), hsamples, eff.names= c()) {
  inference<- list()
   for (i in 1:length(hsamples)) {
     thetas <- hsamples[[i]]$theta
     effects<-list()
     for(j in 1:length(funs)){
       effects [[j]] <- funs[[j]](thetas)
       # names(effects) <- eff.names
     }
     inference[[i]] <- effects
     names(inference[[i]]) <- eff.names
   } 
   final.effects <- list()
   for (k in 1:length(funs)) {
     this.inf <- lapply(inference, function(x) x[[k]])
     inf2 <- unlist(this.inf)
     dim3 <- c(dim(this.inf[[1]]), length(inf2)/prod(dim(this.inf[[1]])))
     dim(inf2) <- dim3
     final.effects[[k]] <- inf2
   }
   final.effects
}


get.msds <- function(samples) {
  av.thetas <-fixedeffects.meanthetas(samples)[[1]]
  msds <- cbind(apply(av.thetas, 2, mean), apply(av.thetas, 2, sd))
  colnames(msds) <- c("M", "SD")
  msds <- data.frame(msds)
  msds
}


cors.plot <- function(out) {
  tmp <- data.frame(cbind(out[[1]],
                                   out[[2]], out[[3]]))
  colnames(tmp)<- c("M", "LCI", "HCI")
  tmp$param <- names(out[[1]])
  ggdf<-tmp
  ggplot(ggdf, aes(param, M))+
  geom_point(size=3)+
  ylab("")+ 
  geom_hline(aes(yintercept=0), linetype=2) +
  geom_errorbar(data=ggdf,aes(ymin=LCI,ymax=HCI)) 
  
}


#averages thetas together that match a certain regex expression
#for a single participant's thetas
get_average_thetas <- function (regex, thetas) {
  apply(thetas[, grep(regex, colnames(thetas)),], c(1,3),
      mean)
}

#averages thetas together that match a certain regex expression
#for .list of participants
h.get_average_thetas <- function (regex, samples) {
  lapply(samples, function(x) get_average_thetas(regex, x$theta))
}




get.subj.effects.m <- function (PPS, fun, names) {
  EFFECTS <- lapply(PPS, get.subj.effects, fun=fun)
  for (i in 1:length(names)) EFFECTS[[i]]$model <- names[i]
  do.call(rbind, EFFECTS)
}

get.subj.effects <- function(PP, fun) {
  for (i in 1:length(PP)){
    effects<-get.subj.pp.MCI (PP[[i]], fun)
    cat(i)
    if (i ==1) out <- effects else out <- rbind(out,effects)
  }
  
  ENS <- rownames(out)
  OUT <- data.frame(out)
  OUT$effect <- ENS
  colnames(OUT) <- c("mean", "lower", "upper", "data", "effect")
  OUT
}


