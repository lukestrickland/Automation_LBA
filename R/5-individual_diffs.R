source("dmc/dmc.R")
source("dmc/dmc_extras.R")
load("samples_data/CA_top_samples.RData")
load_model("LBA", "lba_B.R")

theme_set(theme_simple())

#First get plot of individual participant inhibition/excitation levels

#apply a function to an individual participant's thetas
# and get mean + 95% CIs

individual_summary <- function(theta_FUN, samples){
  individual_thetas <- samples$theta
  effect_summary <- theta_FUN(individual_thetas)
  c(mean= mean(effect_summary), LCI= quantile(effect_summary, 0.025),
    HCI= quantile(effect_summary, 0.975))
}

# Deploy individual_summary to get mean + 95% CIs of a function of thetas
# for participant list

get_individual_summaries <- function(theta_FUN, hsamples, effect_name){
  
  effect_df<- as.data.frame(
    t(
      sapply(FUN=function(x)  individual_summary(x, theta_FUN=theta_FUN), 
        hsamples)
      )
  )
  effect_df$effect <- effect_name
  effect_df$participant <- rownames(effect_df)
  colnames(effect_df) <- c("M", "LCI", "HCI", "effect", "participant")
  effect_df
}


#Inhibition and excitation functions- averaged over the effects

inhibition <- function(x) (
                     (
                        #non-conflict, non-fail, false reduction   
                        (x[,"mean_v.nn.M.nonf.false",,drop=FALSE] - x[,"mean_v.nn.A.nonf.false",,drop=FALSE]) +
                        #non-conflict, fail, true reduction
                        (x[,"mean_v.nn.M.fail.true",,drop=FALSE] - x[,"mean_v.nn.A.fail.true",,drop=FALSE]) +
                        #conflict, non-fail, false reduction  
                        (x[,"mean_v.cc.M.nonf.false",,drop=FALSE] - x[,"mean_v.cc.A.nonf.false",,drop=FALSE]) +
                        #conflict, fail, true reduction
                        (x[,"mean_v.cc.M.fail.true",,drop=FALSE] - x[,"mean_v.cc.A.fail.true",,drop=FALSE]) 
                     )/4
)


excitation <- function(x) (
                     (
                        #non-conflict, non-fail, true increase   
                        (x[,"mean_v.nn.A.nonf.true",,drop=FALSE] - x[,"mean_v.nn.M.nonf.true",,drop=FALSE]) +
                        #non-conflict, fail, false increase
                        (x[,"mean_v.nn.A.fail.false",,drop=FALSE] - x[,"mean_v.nn.M.fail.false",,drop=FALSE]) +
                        #conflict, non-fail, true increase
                        (x[,"mean_v.cc.A.nonf.true",,drop=FALSE] - x[,"mean_v.cc.M.nonf.true",,drop=FALSE]) +
                        #conflict, fail, false increase
                        (x[,"mean_v.cc.A.fail.false",,drop=FALSE] - x[,"mean_v.cc.M.fail.false",,drop=FALSE]) 
                     )/4
)
               




effects_inhibition <- get_individual_summaries(
  theta_FUN = inhibition, hsamples=CA_top_samples,
                         effect_name= "Inhibition")

effects_excitation <- get_individual_summaries(
  theta_FUN = excitation, hsamples=CA_top_samples,
                         effect_name= "Excitation")


all_effects <- rbind(effects_excitation,
                            effects_inhibition)

all_effects$participant <- factor(as.numeric(all_effects$participant))

ggplot(all_effects, aes(participant, M)) +
    geom_point(stat = "identity", size=2.5) +
    geom_errorbar(aes(ymax = HCI , ymin = LCI, width = 0.3)) +
    geom_hline(aes(yintercept=0), linetype=2)+
    facet_grid(.~effect) +xlab("Participant") +ylab ("Effect")


# Next look at correlations across participants between inhibition/excitation
# and automation accuracy effects (benefits on automation correct, 
# costs on automation incorrect)

#Define a function that obtains this measure

automation_effects <- function (currentsim) {

  benefit=NA;cost=NA
  
  nonfail_accuracy_manual <- mean(substr(currentsim$S[currentsim$cond=="M" & currentsim$failtrial=="nonf"],2,2)==
                                    tolower(currentsim$R[currentsim$cond=="M" & currentsim$failtrial=="nonf"]))  
  
  nonfail_accuracy_auto <- mean(substr(currentsim$S[currentsim$cond=="A" & currentsim$failtrial=="nonf"],2,2)==
                                  tolower(currentsim$R[currentsim$cond=="A" & currentsim$failtrial=="nonf"]))  
  
  fail_accuracy_manual <- mean(substr(currentsim$S[currentsim$cond=="M" & currentsim$failtrial=="fail"],2,2)==
                                  tolower(currentsim$R[currentsim$cond=="M" & currentsim$failtrial=="fail"]))  

  fail_accuracy_auto <- mean(substr(currentsim$S[currentsim$cond=="A" & currentsim$failtrial=="fail"],2,2)==
                                  tolower(currentsim$R[currentsim$cond=="A" & currentsim$failtrial=="fail"]))  
  
  out <- c(nonfail_accuracy_auto-nonfail_accuracy_manual, fail_accuracy_manual- fail_accuracy_auto)
  names(out) <- c("benefit", "cost")
  out
  
}

#get list of participant data

data<- get.hdata.dmc(CA_top_samples)

#make data frame containing relevant effects

for (i in unique(data$s)) {
  effects<- automation_effects(data[data$s==i,])
  if (i ==1) out <- effects else out <- rbind(out,effects)
}

#same functions as up above but cor.plausible function collapses all dimensions
# before calculation

inhibition_corplausible <- function(x) (
                     (
                        #non-conflict, non-fail, false reduction   
                        (x["mean_v.nn.M.nonf.false"] - x["mean_v.nn.A.nonf.false"]) +
                        #non-conflict, fail, true reduction
                        (x["mean_v.nn.M.fail.true"] - x["mean_v.nn.A.fail.true"]) +
                        #conflict, non-fail, false reduction  
                        (x["mean_v.cc.M.nonf.false"] - x["mean_v.cc.A.nonf.false"]) +
                        #conflict, fail, true reduction
                        (x["mean_v.cc.M.fail.true"] - x["mean_v.cc.A.fail.true"]) 
                     )/4
)


excitation_corplausible <- function(x) (
                     (
                        #non-conflict non-fail true increase   
                        (x["mean_v.nn.A.nonf.true"] - x["mean_v.nn.M.nonf.true"]) +
                        #non-conflict fail false increase
                        (x["mean_v.nn.A.fail.false"] - x["mean_v.nn.M.fail.false"]) +
                        #conflict non-fail true increase
                        (x["mean_v.cc.A.nonf.true"] - x["mean_v.cc.M.nonf.true"]) +
                        #conflict fail false increase
                        (x["mean_v.cc.A.fail.false"] - x["mean_v.cc.M.fail.false"]) 
                     )/4
)

#See dmc tutorial "plausible" https://osf.io/pbwx8/


#inhibition, costs and benefits

get_corplausible_MCI <- function(samples, cv, p.name, n, fun, kappa=1){
  
  cor.r <- cor.plausible(samples,fun=fun,
                         cv=cv, p.name=p.name)
  
  dens.r <- postRav(r=cor.r,n=n,spacing=.01,kappa=kappa) 
  
  MCIs <- c(postRav.mean(dens.r), postRav.ci(dens.r,interval=c(.025,.975)))
  names(MCIs) <- c("M", "LCI", "HCI")
  MCIs
}

get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible,
                     cv=as.data.frame(out), 
                     p.name="benefit", n=24)

get_corplausible_MCI(CA_top_samples, 
                     fun=inhibition_corplausible,
                     cv=as.data.frame(out), 
                     p.name="cost", n=24)

get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible,
                     cv=as.data.frame(out), 
                     p.name="benefit", n=24)

get_corplausible_MCI(CA_top_samples, 
                     fun=excitation_corplausible,
                     cv=as.data.frame(out), 
                     p.name="cost", n=24)


load("samples_data/CA_top_samples_pp.RData")
load("samples_data/postexp.RData")

subject_postexp<- get.subj.effects.m(list(pp, noinh, noex), automation_effects , 
                          c("Full Model", "Inhibition Removed", "Excitation Removed"))

subject_postexp$effect <- factor(subject_postexp$effect, levels=c("benefit", "cost"),
                      labels = c("Automation Benefit", "Automation Cost"))

subject_postexp$model <- factor(subject_postexp$model, levels=c("Full Model", "Inhibition Removed",
                                          "Excitation Removed"))

ggplot(subject_postexp, aes(data, mean)) + geom_point(size=1) + geom_abline(slope=1, intercept=0) +
facet_grid(effect~model, scales = "free") + geom_errorbar(aes(ymax = upper, ymin = lower), alpha=0.3) +
  ylab("Model") + xlab("Observed")+ theme(text = element_text(size = 20))


