t1 = Sys.time()
library(rmarkdown)
# render("Norm_only.Rmd")
# gc()
# render("Norm_Gpd.Rmd")
# gc()
# render("mix_only.Rmd")
# gc()
render("mix_gpd.Rmd")
t2 = Sys.time()
t2-t1

