library(tidyverse)

extract_one <- function(ppt, sess, block, cond, logdirstr='XML_logs/') {
  
  logFileName <-
    paste(logdirstr,"ppt", ppt, "_sess", sess,
          '_block', block, '_cond_', cond,
          '.xml.log', sep = '')
  
  rawData <- readLines(logFileName)

  ### extract lines containing experimental condition ###
  condition_raw <- rawData[grep("<task_id>cond", rawData)]
  length(condition_raw)
  
  ### extract condition ###
  condition_crop <-
    sapply(strsplit(condition_raw, "<task_id>cond_"), "[", 2)
  condition <- sapply(strsplit(condition_crop, "_"), "[", 1)
  condition <- toupper(condition)
  
  ###################################################################################################
  ################################################################## EXTRACT RESPONSE & RT DATA #####
  ###################################################################################################
  response_lines <- grep(".*</event><elapsed>", rawData)

  ### create object with lines containing raw responses ###
  responses_raw <- rawData[response_lines]
  
  rt_raw <- rawData[response_lines + 1]
  rkey_raw <- rawData[response_lines+2]
  
  rkey <- gsub(".*<key>", "",gsub("</key>.*", "",  rkey_raw))

  ### create object with lines containing raw RTs ###
  RT <- gsub(".*<elapsed>", "",gsub("</elapsed>.*", "",  rt_raw))
  
  ### create object with lines containing missed responses ###
  misses <- rt_raw[grep("decisions", rt_raw)]

  miss_line <- grep("decisions", rt_raw)
  RT[miss_line] <- NA
  RT <- as.numeric(RT)
  ###################################################################################################
  ########################################################### EXCTRACT TRIAL-BY-TRIAL RESPONSES #####
  ###################################################################################################
  
  ### extract responses ###
  responses_crop <- sapply(strsplit(responses_raw, "<event>"), "[", 2)
  response <- sapply(strsplit(responses_crop, "</event>"), "[", 1)
  response <- factor(response)

  ### create trial number vector ###
  trial <- 1:length(response)


  ###################################################################################################
  ############################################ EXTRACT TRIAL-BY-TRIAL AND CUMULATIVE SCORE DATA #####
  ###################################################################################################
  
  ### extract trial-by-trial scores ###
  score_parts <- function(responses_raw) {
    m <- regexec("(<score>([[:print:]]+)</score>)", responses_raw)
    parts <- do.call(rbind,
                     lapply(regmatches(responses_raw, m), `[`, c(3L)))
    colnames(parts) <- c("score")
    parts
  }
  score <- as.integer(score_parts(responses_raw))

  cumulative_score_parts <- function(responses_raw) {
    m <-
      regexec("(<grand_total_score>([[:print:]]+)</grand_total_score>)",
              responses_raw)
    parts <- do.call(rbind,
                     lapply(regmatches(responses_raw, m), `[`, c(3L)))
    colnames(parts) <- c("cumulative_score")
    parts
  }
  
  cumulative_score <-
    as.integer(cumulative_score_parts(responses_raw))
  
  maps_and_aircraft <-read.csv(
    paste('mapac_data/exp_vars_p', ppt, '_s', sess, "_",
          cond, '.csv', sep = ''), header=TRUE
  )[,-1]
  
  maps_and_aircraft[order(maps_and_aircraft$presOrder),]

  cbind(tibble(ppt, sess, block, cond, response, rkey, RT, score, cumulative_score),
        maps_and_aircraft[order(maps_and_aircraft$presOrder),])
}