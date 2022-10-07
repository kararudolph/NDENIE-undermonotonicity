suppressPackageStartupMessages(library(tidyverse))
library(kableExtra)

read_zip <- function(tar) {
  files <- unzip(tar, list = TRUE)$Name
  p <- progressr::progressor(along = 1:length(files))
  purrr::map(files, function(file) {
    p()
    con <- unz(tar, file)
    read.csv(con)
  })
}

nie <- 0.034
nde <- 0.099

local({
  source("_research/dgp_cont.R")
  nde <<- truth_dgp_cont()["nde"]
  nie <<- truth_dgp_cont()["nie"]
})

basepath <- "_research/data/sim_binary_gcomp_"

res <- map_dfr(1:4, \(x) bind_rows(read_zip(paste0(basepath, x, ".zip"))))

res <- filter(res, estimand == "nie") |> 
  group_by(n, spec) |> 
  summarise(abs_bias = abs(mean(psi) - nie), 
            coverage = mean(map2_lgl(conf_low, conf_high, \(l, h) between(nie, l, h)))) |> 
  mutate(rootn_bias = sqrt(n) * abs_bias, 
         estimand = "NIE") |> 
  select(estimand, spec, n, everything()) |> 
  arrange(spec) |> 
  bind_rows({
    filter(res, estimand == "nde") |> 
      group_by(n, spec) |> 
      summarise(abs_bias = abs(mean(psi) - nde), 
                coverage = mean(map2_lgl(conf_low, conf_high, \(l, h) between(nde, l, h)))) |> 
      mutate(rootn_bias = sqrt(n) * abs_bias, 
             estimand = "NDE") |> 
      select(estimand, n, spec, everything()) |> 
      arrange(spec)
  })

kbl(res, 
    format = "latex", 
    align = "lllccc",
    col.names = c("Estimand", "Model misspec.", "N", "$|Bias|$", 
                  "95\\% cov.", "$\\sqrt{N} \\times |Bias|$"),
    booktabs = TRUE, 
    linesep = c(rep(c("", "\\addlinespace"), 3), "", "\\midrule"), 
    escape = FALSE, 
    digits = 2)
