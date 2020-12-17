
rm(list=ls())
load("samples_data/CA_top_thresholdsmult_samples_finished.RData")

source("dmc/dmc.R")
source("dmc/dmc_extras.R")

load_model("LBA", "lba_B_autothres.R")

sapply(CA_top_thresholdsmult_samples, function(x) attr(x, "auto"))

#Drop participant who doesn't converge

CA_top_thresholdsmult_samples <- CA_top_thresholdsmult_samples[!(
  names(CA_top_thresholdsmult_samples)=="16")]

sapply(CA_top_thresholdsmult_samples, function(x) attr(x, "auto"))

msds <- get.msds(CA_top_thresholdsmult_samples)

#Threshold multiplier parameters:
#Coded as follows:
# A = automation condition (they're all for automation condition)
# first (lower) letter = stimulus type
# second (upper) letter = response
# last letter (S/F) = Automation success or failure
# Thus, the first one AnNS = the multiplier (applied to manual thresholds)
# that determines the non-conflict accumulator threshold on non-conflict trials
# on instances that automation succeeds (recommends non-conflcit).
# If threshold updates were adaptive, or participants gained an information
# boost initially, you would expect thresholds to be lower in this condition than
# in manual (hence multiplier <1)


msds[grep("tb", rownames(msds)),]

#Examine inhibition/excitation to see if inferences affected by
# this extra mechanism - they're not

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


ggplot(Vs, aes(factor(Auto),M)) + 
  geom_point(stat = "identity",aes(shape=Condition, col=Condition), size=2.5) +
  geom_errorbar(aes(ymax = M + SD, ymin = M - SD, width = 0.3, col=Condition))+ 
  ylab("Accumulation Rate") + xlab("")+
  facet_grid(S ~ match,scales = "free", space = "free") +
    theme(text = element_text(size = 12))



#Check - does fit look ok? Yes looks about the same as model reported in text
pp <- h.post.predict.dmc(CA_top_thresholdsmult_samples, save.simulation = TRUE, cores=8)

theme_set(theme_simple())

rescore_column <- function(df) {
  df$R <- factor(as.character(toupper(substr(df$S,1,1))==df$R))
  new_data <- attr(df, "data")
  new_data$R <- factor(as.character(toupper(substr(new_data$S,1,1))==new_data$R))
#  new_data <- new_data %>% select(C, everything())
  #%>% select(-R)
  attr(df, "data") <- new_data
#  df %>% select(reps, C, everything())
  #%>% select(-R)
  df
}

pp1 <- lapply(pp, rescore_column)


fitlist <- GET.fitgglist.dmc(pp1, factors=c("cond", "failtrial"))
save(pp, fitlist,
     file="samples_data/CA_top_thresholdsmult_samples_pp.RData")

load("samples_data/CA_top_thresholdsmult_samples_pp.RData")



accs <- fitlist$pps %>% filter(R=="TRUE") %>% select(-R)

accs$cond <- factor(accs$cond, levels=c("A", "M"), labels =
                      c("Automation Condition", "Manual Condition"))

accs$failtrial <- factor(accs$failtrial, levels=c("nonf", "fail"), labels =
                      c("Automation Success", "Automation Failure"))

plot1 <- ggplot.RP.dmc(accs, xaxis="cond") +xlab("") +ylab("Accuracy")

corRTs <- fitlist$RTs %>% filter(R=="TRUE") %>% select(-R)

corRTs$cond <- factor(corRTs$cond, levels=c("A", "M"), labels =
                      c("Automation Condition", "Manual Condition"))

corRTs$failtrial <- factor(corRTs$failtrial, levels=c("nonf", "fail"), labels =
                      c("Automation Success", "Automation Failure"))

plot2 <- ggplot.RT.dmc(corRTs, xaxis="cond") +xlab("") +ylab("Correct RT")



errRTs <- fitlist$RTs %>% filter(R=="FALSE") %>% select(-R)

errRTs$cond <- factor(errRTs$cond, levels=c("A", "M"), labels =
                      c("Automation Condition", "Manual Condition"))

errRTs$failtrial <- factor(errRTs$failtrial, levels=c("nonf", "fail"), labels =
                      c("Automation Success", "Automation Failure"))

plot3 <- ggplot.RT.dmc(errRTs, xaxis="cond") +xlab("") +ylab("Error RT")


grid.arrange(plot1,plot2,plot3)