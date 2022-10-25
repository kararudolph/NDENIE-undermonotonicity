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

dgp <- "mult"
estimator <- "gcomp"

local({
  source(glue::glue("_research/dgp_{dgp}.R"))
  nde <<- truth()["nde"]
  nie <<- truth()["nie"]
  bounds <<- efficiency_bound(1e6)
})

basepath <- glue::glue("_research/data/sim_{dgp}_{estimator}_")

if (estimator == "onestep") {
  res <- map_dfr(c(1:2, 4:5), \(x) bind_rows(read_zip(paste0(basepath, x, ".zip"))))
  
  res <- mutate(res, 
                std_error = case_when(
                  spec == 5 ~ sqrt(var),
                  TRUE ~ (psi - conf_low) / qnorm(0.975),
                  estimator == "gcomp" ~ ((conf_high - conf_low) / 2) / qnorm(0.975) 
                ))
} else {
  res <- map_dfr(1:4, \(x) bind_rows(read_zip(paste0(basepath, x, ".zip"))))
  res <- mutate(res, std_error = ((conf_high - conf_low) / 2) / qnorm(0.975))
}

res <- filter(res, estimand == "nie") |> 
  group_by(n, spec) |> 
  summarise(abs_bias = abs(mean(psi) - nie), 
            relse = mean(std_error / sd(psi)),
            relsd = mean(std_error / sqrt(bounds$nie / n)),
            relmse = mean((psi - nie)^2 / bounds$nie),
            coverage = mean(map2_lgl(conf_low, conf_high, \(l, h) between(nie, l, h)))) |> 
  mutate(rootn_bias = sqrt(n) * abs_bias, 
         relmse = relmse * n,
         estimand = "NIE") |> 
  select(estimand, spec, n, everything()) |> 
  arrange(spec) |> 
  bind_rows({
    filter(res, estimand == "nde") |> 
      group_by(n, spec) |> 
      summarise(abs_bias = abs(mean(psi) - nde), 
                relse = mean(std_error) / sd(psi),
                relsd = mean(std_error / sqrt(bounds$nde / n)),
                relmse = mean((psi - nde)^2 / bounds$nde),
                coverage = mean(map2_lgl(conf_low, conf_high, \(l, h) between(nde, l, h)))) |> 
      mutate(rootn_bias = sqrt(n) * abs_bias, 
             relmse = relmse * n,
             estimand = "NDE") |> 
      select(estimand, n, spec, everything()) |> 
      arrange(spec)
  })

res <- mutate(res, 
       spec = case_when(
         spec == 1 ~ "-", 
         spec == 2 ~ "$\\g, \\q, \\e, \\rr$", 
         spec == 4 ~ "$\\mu, \\rho, \\q$", 
         spec == 5 ~ "$\\g, \\q, \\mu$", 
         spec == "correct" ~ "-", 
         spec == "m incorrect" ~ "$\\mu, \\q$", 
         spec == "y incorrect" ~ "$\\rho, \\q$", 
         spec == "z incorrect" ~ "$\\mu, \\rho$"
       ))

kbl(select(res, -estimand), 
    format = "latex", 
    align = "cccccccc",
    col.names = c("Correctly specified", "N", "$|Bias|$", 
                  "relse", "relsd", "relmse", 
                  "95\\% cov.", "$\\sqrt{N} \\times |Bias|$"),
    booktabs = TRUE, 
    linesep = c(rep(c("", "\\addlinespace"), 3), "", "\\midrule"), 
    escape = FALSE, 
    digits = 2) |> 
  pack_rows("NIE", 1, 8) %>%
  pack_rows("NDE", 9, 16)
