suppressPackageStartupMessages(library(tidyverse))

read_zip <- function(tar) {
  files <- unzip(tar, list = TRUE)$Name
  p <- progressr::progressor(along = 1:length(files))
  purrr::map(files, function(file) {
    p()
    con <- unz(tar, file)
    read.csv(con)
  })
}

basepath <- "_research/data/sim_binary_onestep_"

res <- map_dfr(1:4, \(x) bind_rows(read_zip(paste0(basepath, x, ".zip"))))

mutate(res, 
       abs_bias = case_when(estimand == "nie" ~ abs(psi - 0.034), 
                            TRUE ~ abs(psi - 0.099)), 
       coverage = case_when(estimand == "nie" ~ 
                              map2_lgl(conf_low, conf_high, \(l, h) between(0.034, l, h)), 
                            TRUE ~ map2_lgl(conf_low, conf_high, \(l, h) between(0.099, l, h)))) |> 
  group_by(n, spec, estimand) |> 
  summarise(abs_bias = mean(abs_bias), 
            coverage = mean(coverage)) |> 
  mutate(rootn_bias = abs_bias * sqrt(n)) |> 
  arrange(spec, estimand, n)
