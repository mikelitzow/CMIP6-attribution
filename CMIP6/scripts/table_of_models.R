# create table of included CMIP6 models for appendix
library(tidyverse)


models.245 <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp245") %>%
  str_remove_all("_245.nc")
  
models.585 <- list.files("./CMIP6/CMIP6_outputs/1850-2099_runs/ssp585") %>%
  str_remove_all(".nc")

models1 <- data.frame(Model = models.585,
                      "SSP5-8.5" = "X")

models2 <- data.frame(Model = models.245,
                      "SSP2-4.5" = "X")

model.table <- left_join(models1, models2)

write.csv(model.table, "appendix_model_table.csv", row.names = F)
