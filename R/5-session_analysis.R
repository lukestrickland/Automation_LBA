source("dmc/dmc.R")
source("dmc/dmc_extras.R")
load("samples_data/CA_sess_samples_finished.RData")

msds <- get.msds(CA_sess_samples)

Vs <- msds[grep("mean_v", rownames(msds)),]

Vs$Cond <- "Manual" 
Vs$Cond[grep("A", rownames(Vs))] <- "Automation"
Vs$Auto <- "Automation Correct"
Vs$Auto[grep("fail", rownames(Vs))] <- "Automation Incorrect"
Vs$S <- "Conflict"
Vs$S[grep("nn", rownames(Vs))] <- "Non-conflict"
Vs$match <- "Match"
Vs$match[grep("false", rownames(Vs))] <- "Mismatch"

names(Vs)[names(Vs)=="Cond"] <- "Condition"

Vs$exinh <- NA
Vs$exinh[Vs$Auto=="Automation Correct" & Vs$match=="Match"] <- "Excitation"
Vs$exinh[Vs$Auto=="Automation Incorrect" & Vs$match=="Match"] <- "Inhibition"
Vs$exinh[Vs$Auto=="Automation Correct" & Vs$match=="Mismatch"] <- "Inhibition"
Vs$exinh[Vs$Auto=="Automation Incorrect" & Vs$match=="Mismatch"] <- "Excitation"

Vs_1 <- Vs[grep("1", rownames(Vs)),]


ggplot(Vs_1, aes(factor(Auto),M)) + 
  geom_point(stat = "identity",aes(shape=Condition, col=Condition), size=2.5) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.3, col=Condition))+ 
  ylab("Accumulation Rate") + xlab("")+
  facet_grid(S ~ match,scales = "free", space = "free") + ggtitle("Session One") +
    theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5))
            

Vs_2 <- Vs[grep("2", rownames(Vs)),]

ggplot(Vs_2, aes(factor(Auto),M)) + 
  geom_point(stat = "identity",aes(shape=Condition, col=Condition), size=2.5) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.3, col=Condition))+ 
  ylab("Accumulation Rate") + xlab("")+
  facet_grid(S ~ match,scales = "free", space = "free") + ggtitle("Session Two") +
    theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5))


Bs <- msds[grep("B", rownames(msds)),]

Bs$R <- "Non-conflict"
Bs$R[grep("C", rownames(Bs))] <- "Conflict"
Bs$Cond <- "Automation"
Bs$Cond[grep("M", rownames(Bs))] <- "Manual"
Bs$Session <- "Session One"
Bs$Session[grep("2", rownames(Bs))] <- "Session Two"

names(Bs)[names(Bs)=="Cond"] <- "Condition"

ggplot(Bs, aes(factor(R),M)) + 
  geom_point(stat = "identity",aes(col=Condition, shape=Condition), size=2.5) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.3, col=Condition))+ 
  ylab("Threshold") + xlab("Accumulator")+
  facet_grid(.~Session)
