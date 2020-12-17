#This Script sets up the final reported model set for the paper which are then
# dispatched on a grid system using pbs pro (see grid_dispatch.R)

##Create model-ready data frame and load dmc functions

#load dmc
source("dmc/dmc.R")

create_model_data <- function(file) {
  load(file)
  cleandats <- cleandats[!colnames(cleandats) %in% "C"]
  cleandats <- as.data.frame(cleandats)
  cleandats$cond <- factor(cleandats$cond, levels=c("AUTO", "MANUAL"),
                           labels=c("A", "M"))
  
  cleandats$S <- factor(cleandats$S, levels=c("n", "c"),
                           labels=c("nn", "cc"))
  
  cleandats$R <- factor(cleandats$R, levels=c("N", "C"))
  cleandats$s<- factor(cleandats$s)
  cleandats
}

cleandats <- create_model_data("img/cleandats.RData")


#Full classic LBA model of the experiment
load_model("LBA", "lba_B.R")

CA_top_model <- model.dmc(
  p.map = list(
    A = "1",B = c("cond", "sess", "R"), t0 = "1", mean_v = c("S", "cond", "failtrial", "M"),
    sd_v = c("M"), st0 = "1"),
  match.map = list(
    M = list(nn = "N", cc="C")
  ),
  factors = list(
    S = c("nn", "cc"), cond = c("A", "M"), sess = c("1", "2"),
    failtrial=c("nonf", "fail")
  ),
  constants = c(st0 = 0, sd_v.false = 1
  ),
  responses = c("N", "C"),type = "norm"
)


CA_top_p.vector  <- c(t0=0.3,A=3,
                                sd_v.true = 1,
               
  B.A.1.N=2, B.M.1.N=2,             
   B.A.2.N=2, B.M.2.N=2, B.A.1.C=2,             
  B.M.1.C=2, B.A.2.C=2,   B.M.2.C=2, 
  
  mean_v.nn.A.nonf.true=1,  mean_v.cc.A.nonf.true=1, 
 mean_v.nn.M.nonf.true=1,  mean_v.cc.M.nonf.true=1,
 mean_v.nn.A.fail.true=1, mean_v.cc.A.fail.true=1, 
 mean_v.nn.M.fail.true=1,  mean_v.cc.M.fail.true=1, 
 mean_v.nn.A.nonf.false=0, mean_v.cc.A.nonf.false=0,
 mean_v.nn.M.nonf.false=0,mean_v.cc.M.nonf.false=0, 
 mean_v.nn.A.fail.false=0, mean_v.cc.A.fail.false=0,
 mean_v.nn.M.fail.false=0, mean_v.cc.M.fail.false=0
 )

check.p.vector(CA_top_p.vector, CA_top_model)

CA_top_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(CA_top_p.vector)),
  p1=CA_top_p.vector,                           
  p2=c(1,1,1,rep(1, 8), rep(2, 16)),
  lower=c(0.1, 0,0, rep(0, 8), rep(NA, 16)),
  upper=c(5,10, rep(Inf, length(CA_top_p.vector)-2))
)

CA_top_dm <- data.model.dmc(cleandats,
                                   CA_top_model)

CA_top_samples <- h.samples.dmc(nmc = 180,
                                          CA_top_p.prior,
                                          CA_top_dm, thin=20)

save(CA_top_samples, file="CA_top_samples.RData")

# Reviewer comment: Inhibition/excitation mechanisms by day (practice effects?)

CA_sess_model <- model.dmc(
  p.map = list(
    A = "1",B = c("cond", "sess", "R"), t0 = "1", mean_v = c("S", "cond", "sess", "failtrial", "M"),
    sd_v = c("M"), st0 = "1"),
  match.map = list(
    M = list(nn = "N", cc="C")
  ),
  factors = list(
    S = c("nn", "cc"), cond = c("A", "M"), sess = c("1", "2"),
    failtrial=c("nonf", "fail")
  ),
  constants = c(st0 = 0, sd_v.false = 1
  ),
  responses = c("N", "C"),type = "norm"
)


CA_sess_p.vector  <- c(t0=0.3,A=3,
                                sd_v.true = 1,
               
  B.A.1.N=2, B.M.1.N=2,             
   B.A.2.N=2, B.M.2.N=2, B.A.1.C=2,             
  B.M.1.C=2, B.A.2.C=2,   B.M.2.C=2, 
  
  mean_v.nn.A.1.nonf.true=1,  mean_v.cc.A.1.nonf.true=1, 
 mean_v.nn.M.1.nonf.true=1,  mean_v.cc.M.1.nonf.true=1,
 mean_v.nn.A.1.fail.true=1, mean_v.cc.A.1.fail.true=1, 
 mean_v.nn.M.1.fail.true=1,  mean_v.cc.M.1.fail.true=1, 
 mean_v.nn.A.1.nonf.false=0, mean_v.cc.A.1.nonf.false=0,
 mean_v.nn.M.1.nonf.false=0,mean_v.cc.M.1.nonf.false=0, 
 mean_v.nn.A.1.fail.false=0, mean_v.cc.A.1.fail.false=0,
 mean_v.nn.M.1.fail.false=0, mean_v.cc.M.1.fail.false=0,
 
 
 mean_v.nn.A.2.nonf.true=1,  mean_v.cc.A.2.nonf.true=1, 
 mean_v.nn.M.2.nonf.true=1,  mean_v.cc.M.2.nonf.true=1,
 mean_v.nn.A.2.fail.true=1, mean_v.cc.A.2.fail.true=1, 
 mean_v.nn.M.2.fail.true=1,  mean_v.cc.M.2.fail.true=1, 
 mean_v.nn.A.2.nonf.false=0, mean_v.cc.A.2.nonf.false=0,
 mean_v.nn.M.2.nonf.false=0,mean_v.cc.M.2.nonf.false=0, 
 mean_v.nn.A.2.fail.false=0, mean_v.cc.A.2.fail.false=0,
 mean_v.nn.M.2.fail.false=0, mean_v.cc.M.2.fail.false=0
 
 
 )

check.p.vector(CA_sess_p.vector, CA_sess_model)

CA_sess_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(CA_sess_p.vector)),
  p1=CA_sess_p.vector,                           
  p2=c(1,1,1,rep(1, 8), rep(2, 32)),
  lower=c(0.1, 0,0, rep(0, 8), rep(NA, 32)),
  upper=c(5,10, rep(Inf, length(CA_sess_p.vector)-2))
)

CA_sess_dm <- data.model.dmc(cleandats,
                                   CA_sess_model)

CA_sess_samples <- h.samples.dmc(nmc = 180,
                                          CA_sess_p.prior,
                                          CA_sess_dm, thin=20)

save(CA_sess_samples, file="CA_sess_samples.RData")

# Reviewer comment: Threshold shifts to approximate a quick information boost

# Thresholds can be affected by automation via a threshold multiplier
# use of multiplier assures we can keep the constraint that b>A

load_model("LBA", "lba_B_autothres.R")

tmap <-
  empty.map(list(
         S = c("nn", "cc"), cond = c("A", "M"), sess = c("1", "2"),
            failtrial=c("nonf", "fail"),
             R = c("N", "C")
  ), 
    levels=c(
             "MAN",
             "AnNS","AnNF",
             "AnCS", "AnCF",
              "AcNS","AcNF",
             "AcCS", "AcCF"
             ))

tmap[1:32] <- c(
  "AnNS","AcNS","MAN", "MAN",
  "AnNS","AcNS","MAN", "MAN",
  
  "AnNF","AcNF","MAN", "MAN",
  "AnNF","AcNF","MAN", "MAN",
  
  "AnCS","AcCS","MAN", "MAN",
  "AnCS","AcCS","MAN", "MAN",
  
  "AnCF","AcCF","MAN", "MAN",
  "AnCF","AcCF","MAN", "MAN"
  
)


CA_top_thresholdsmult_model <- model.dmc(
  p.map = list(
    A = "1",B = c("sess", "R"), t0 = "1", mean_v = c("S", "cond", "failtrial", "M"),
    sd_v = c("M"), st0 = "1", tb="TMAP"),
  match.map = list(
    M = list(nn = "N", cc="C"),
    TMAP=tmap
  ),
  factors = list(
    S = c("nn", "cc"), cond = c("A", "M"), sess = c("1", "2"),
    failtrial=c("nonf", "fail")
  ),
  constants = c(st0 = 0, sd_v.false = 1, tb.MAN=1
  ),
  responses = c("N", "C"),type = "norm"
)


CA_top_thresholdsmult_p.vector  <- c(t0=0.3,A=3,
                                sd_v.true = 1,
               
  B.1.N=2,B.2.N=2, B.1.C=2, B.2.C=2, 
          
 tb.AnNS=1,                tb.AnNF=1,                tb.AnCS=1,               
tb.AnCF=1,                tb.AcNS=1,                tb.AcNF=1,               
 tb.AcCS=1,                tb.AcCF=1,   
  
  mean_v.nn.A.nonf.true=1,  mean_v.cc.A.nonf.true=1, 
 mean_v.nn.M.nonf.true=1,  mean_v.cc.M.nonf.true=1,
 mean_v.nn.A.fail.true=1, mean_v.cc.A.fail.true=1, 
 mean_v.nn.M.fail.true=1,  mean_v.cc.M.fail.true=1, 
 mean_v.nn.A.nonf.false=0, mean_v.cc.A.nonf.false=0,
 mean_v.nn.M.nonf.false=0,mean_v.cc.M.nonf.false=0, 
 mean_v.nn.A.fail.false=0, mean_v.cc.A.fail.false=0,
 mean_v.nn.M.fail.false=0, mean_v.cc.M.fail.false=0
 )

check.p.vector(CA_top_thresholdsmult_p.vector, CA_top_thresholdsmult_model)

CA_top_thresholdsmult_p.prior <- prior.p.dmc(
  dists = rep("tnorm", length(CA_top_thresholdsmult_p.vector)),
  p1=CA_top_thresholdsmult_p.vector,                           
  p2=c(1,1,1,rep(1, 12), rep(2, 16)),
  lower=c(0.1, 0,0, rep(0, 12), rep(NA, 16)),
  upper=c(5,10, rep(Inf, length(CA_top_thresholdsmult_p.vector)-2))
)

CA_top_thresholdsmult_dm <- data.model.dmc(cleandats,
                                   CA_top_thresholdsmult_model)

CA_top_thresholdsmult_samples <- h.samples.dmc(nmc = 180,
                                          CA_top_thresholdsmult_p.prior,
                                          CA_top_thresholdsmult_dm, thin=20)

save(CA_top_thresholdsmult_samples, file="CA_top_thresholdsmult_samples.RData")








