
source("dmc/dmc.R")

load_model("LBA", "lba_B_autothres.R")
run.grid.dmc("CA_top_thresholdsmult_samples",model.dir ="LBA",
             model.file="lba_B_autothres.R",user="ljs392",
             n.add=60, wall.hours = 300,
             GB = 3, max.try=30)