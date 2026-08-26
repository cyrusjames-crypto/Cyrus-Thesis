# ============================================================================
# EpiSewer starter - RSV effective reproduction number (Rt) for ONE NWSS site
# Week-1 goal: prove the pipeline end-to-end on a single site, THEN scale to N.
# Built on the EpiSewer getting-started README (package v0.0.4, Oct 2025).
# Numbers marked <CONFIRM> are placeholders - set them from the literature with
# your advisor / Bayesian expert before trusting any output.
# ============================================================================

# ---- 0. One-time setup (run once per machine) ------------------------------
#remotes::install_github("adrian-lison/EpiSewer", dependencies = TRUE)
#cmdstanr::check_cmdstan_toolchain()
# cmdstanr::install_cmdstan(cores = 2)
#EpiSewer::sewer_compile()      # compiles the Stan models (once per install/update)
# 1. the missing piece - install remotes from CRAN
#install.packages("remotes")

# 2. cmdstanr is not on CRAN; get it from the Stan R-universe
#install.packages("cmdstanr",
#                 repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

# 3. now this will work
#remotes::install_github("adrian-lison/EpiSewer", dependencies = TRUE)
#cmdstanr::check_cmdstan_toolchain()   # checks your C++ toolchain first
#cmdstanr::install_cmdstan(cores = 2)  # downloads + compiles CmdStan — takes a few minutes
#EpiSewer::sewer_compile()             # compiles EpiSewer's models — the real green light
# 4. and this
library(EpiSewer)
library(EpiSewer)
library(data.table)
library(ggplot2)

# ---- 1. Load ONE site's wrangled NWSS data ---------------------------------
# EpiSewer expects two tables:
#   measurements: date (Date), concentration (numeric, gc/mL)
#   flows:        date (Date), flow          (numeric, mL/day)
#
# UNIT / FORMAT CHECKS before trusting anything:
#  - concentration and flow MUST share a volume unit (gc/mL pairs with mL/day).
#    NWSS often reports gc/L (or gc/g dry weight for solids) - convert to gc/mL.
#  - every measured date needs a flow value on that date (impute flow or drop).
#  - DO NOT KNN-impute non-detects. EpiSewer handles missing days natively, and
#    has a dedicated LOD + dPCR noise model for non-detects (see section 7).

library(EpiSewer); library(data.table); library(ggplot2)

# filter to RSV plant-influent composites
rsv <- dt[pcr_target == "rsv" & sample_location == "wwtp" & sample_type %like% "composite"]

# pick the site
site <- "fl_1702_wwtp"
one  <- rsv[key_plot == site &
              sample_collect_date >= as.Date("2024-08-01") &
              sample_collect_date <= as.Date("2025-06-30")]

# measurements — now using the REAL concentration (0 = non-detect), not the LOD/2 column
measurements <- one[, .(date = as.Date(sample_collect_date),
                        concentration = pcr_target_avg_conc)][order(date)]

# constant flow (this site's flow is mostly NA)
flow_const <- mean(one$capacity_mgd, na.rm = TRUE) * 3785411.784
flows <- measurements[, .(date, flow = flow_const)]

ww_data <- sewer_data(measurements = measurements, flows = flows)
plot(measurements$date, measurements$concentration, type = "b",
     xlab = "Date", ylab = "RSV copies/L", main = "fl_1702_wwtp 2024-25")

gen_rsv  <- get_discrete_gamma_shifted(gamma_mean = 7.5, gamma_sd = 2.5)  # RSV serial interval ~7-9 d (CONFIRM)
shed_rsv <- get_discrete_gamma(gamma_shape = 5.29, gamma_scale = 0.8695652)      # RSV shedding mean ~4.5 d (Kenyan study)

one[, .(n = .N, first = min(sample_collect_date), last = max(sample_collect_date))]
one[, .(n = .N, first = min(sample_collect_date), last = max(sample_collect_date)),
    by = concentration_method]
one <- rsv[key_plot == site &
             sample_collect_date >= as.Date("2024-08-01") &
             sample_collect_date <= as.Date("2025-06-30") &
             concentration_method == "none"]        # 218 samples, full season
library(EpiSewer); library(data.table); library(ggplot2)

# filter to RSV plant-influent composites
rsv <- dt[pcr_target == "rsv" & sample_location == "wwtp" & sample_type %like% "composite"]

# pick the site + ONE concentration method (none = 218 samples, full season)
site <- "fl_1702_wwtp"
one  <- rsv[key_plot == site &
              sample_collect_date >= as.Date("2024-08-01") &
              sample_collect_date <= as.Date("2025-06-30") &
              concentration_method == "none"]

# measurements — real concentration (0 = non-detect), not the LOD/2 column
measurements <- one[, .(date = as.Date(sample_collect_date),
                        concentration = pcr_target_avg_conc)][order(date)]

# guard: stop if any date still has >1 row
stopifnot(nrow(measurements) == uniqueN(measurements$date))

# constant flow (this site's flow_rate is mostly NA)
flow_const <- mean(one$capacity_mgd, na.rm = TRUE) * 3785411.784
flows <- measurements[, .(date, flow = flow_const)]

# build the EpiSewer data object
ww_data <- sewer_data(measurements = measurements, flows = flows)

# eyeball the signal
plot(measurements$date, measurements$concentration, type = "b",
     xlab = "Date", ylab = "RSV copies/L", main = "fl_1702_wwtp 2024-25 (none)")
options(mc.cores = 4)
fit_rsv <- EpiSewer(data = ww_data, assumptions = assume_rsv,
                    fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                      iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42)))

fit_rsv$fitted$diagnostic_summary()   # want num_divergent = 0
plot_R(fit_rsv)
# assumptions
gen_rsv  <- get_discrete_gamma_shifted(gamma_mean = 7.5, gamma_sd = 2.5)   # RSV serial interval ~7-9 d (CONFIRM)
shed_rsv <- get_discrete_gamma(gamma_shape = 5.29, gamma_scale = 0.8695652) # RSV shedding mean ~4.6 d, CV ~0.43
assume_rsv <- sewer_assumptions(generation_dist = gen_rsv,
                                shedding_dist = shed_rsv,
                                shedding_reference = "infection")
# assumptions

#Generation time distribution: Distribution of the time between a primary 
# infection and its resulting secondary infections
gen_rsv  <- get_discrete_gamma_shifted(gamma_mean = 7.5, gamma_sd = 2.1)
# Shedding load distribution: Distribution of the load shed by an average 
# individual over time. We also specify that our distribution is relative 
# to the day of symptom onset.
# COVID: gamma_shape = 0.929639, gamma_scale = 7.241397
# mean = 6.7, SD = 6.98
# RSV: days 4--12, mean=8\
# from Houston paper, mean =4.6, sd=2
mean.shed <- 4.6
sd.shed <- 2
ga.shape<-(mean.shed/sd.shed)^2
ga.scale<-sd.shed^2/mean.shed
c(ga.shape,ga.scale)
shed_rsv <- get_discrete_gamma(gamma_shape = ga.shape, gamma_scale = ga.scale)
assume_rsv <- sewer_assumptions(generation_dist = gen_rsv, 
                                shedding_dist = shed_rsv, 
                                shedding_reference = "infection")

#try new location in mass to potentially compare to state data
rsv[wwtp_jurisdiction == "ma", .N, by = key_plot][order(-N)]
ma_site <- "..."   # a well-covered key_plot from the list above

ma_site <- "ma_2032_wwtp"

rsv[key_plot == ma_site,
    .(n = .N, first = min(sample_collect_date), last = max(sample_collect_date)),
    by = concentration_method][order(-n)]

rsv[key_plot == ma_site, .(flow_present = sum(!is.na(flow_rate)), total = .N)]

rsv[key_plot == ma_site, .N, by = .(pcr_target_units, sample_matrix, pcr_type)]
#check for any liquid influent]
rsv[wwtp_jurisdiction == "ma", .N, by = .(key_plot, pcr_target_units, sample_matrix)][order(-N)]
#check to see whether there was flow data for this site
ma_site <- "ma_2539_wwtp"
#double check what type of flow data compared to florida site
rsv[key_plot == ma_site,
    .(n = .N, first = min(sample_collect_date), last = max(sample_collect_date)),
    by = concentration_method][order(-n)]

rsv[key_plot == ma_site,
    .(n = .N, first = min(sample_collect_date), last = max(sample_collect_date)),
    by = concentration_method][order(-n)]
ma_site <- "ma_2539_wwtp"
#compare sites florida to mass, using same 
one_ma <- rsv[key_plot == ma_site &
                sample_collect_date >= as.Date("2024-08-01") &
                sample_collect_date <= as.Date("2025-06-30") &
                concentration_method == "ceres nanotrap"]

measurements <- one_ma[, .(date = as.Date(sample_collect_date),
                           concentration = pcr_target_avg_conc)][order(date)]
stopifnot(nrow(measurements) == uniqueN(measurements$date))

flows <- one_ma[, .(flow = mean(flow_rate, na.rm = TRUE) * 3785411.784),
                by = .(date = as.Date(sample_collect_date))][order(date)]

ww_data_ma <- sewer_data(measurements = measurements, flows = flows)

plot(measurements$date, measurements$concentration, type = "b",
     xlab = "Date", ylab = "RSV copies/L", main = "ma_2539_wwtp 2024-25")

options(mc.cores = 4)
fit_ma <- EpiSewer(data = ww_data_ma, assumptions = assume_rsv,
                   fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                     iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42)))
fit_ma$fitted$diagnostic_summary()
plot_R(fit_ma)
#double check to see whether assumptions are the same shedding distributions
assume_rsv$shedding$shedding_dist   # the shedding distribution passed to both fits

str(fit_rsv$job$assumptions$shedding)   # Florida
str(fit_ma$job$assumptions$shedding)    # Massachusetts
str(assume_rsv, max.level = 2)
#run overlap to see how they fit
fl <- fit_rsv$summary$R[, .(date, median, site = "FL (fl_1702)")]
ma <- fit_ma$summary$R[,  .(date, median, site = "MA (ma_2539)")]
ggplot(rbind(fl, ma), aes(date, median, colour = site)) +
  geom_line(linewidth = 1) + geom_hline(yintercept = 1, linetype = 2) +
  labs(y = "Rt", colour = NULL) + theme_minimal()

#check to see when MA site starts recording data
range(fit_ma$summary$R$date)          # where MA's estimates actually begin
range(measurements$date)              # MA's actual sample dates (if measurements is still MA's)
trim <- function(x, days = 10) x[date >= min(date) + days]
fl <- trim(fit_rsv$summary$R)[, .(date, median, site = "FL (fl_1702)")]
ma <- trim(fit_ma$summary$R)[,  .(date, median, site = "MA (ma_2539)")]
ggplot(rbind(fl, ma), aes(date, median, colour = site)) +
  geom_line(linewidth = 1) + geom_hline(yintercept = 1, linetype = 2) +
  labs(y = "Rt", colour = NULL) + theme_minimal()
plot(measurements$date, measurements$concentration, type = "b",
     xlab = "Date", ylab = "RSV copies/L", main = "ma_2539 2024-25")
#check what interval columns exist
names(fit_ma$summary$R)
#check uncertainty 
prep <- function(fit, label) {
  fit$summary$R[seeding == FALSE,
                .(date, median, lower_0.95, upper_0.95, site = label)]
}
fl <- prep(fit_rsv, "FL (fl_1702)")
ma <- prep(fit_ma,  "MA (ma_2539)")
both <- rbind(fl, ma)

ggplot(both, aes(date, median, colour = site, fill = site)) +
  geom_ribbon(aes(ymin = lower_0.95, ymax = upper_0.95), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(y = expression(R[t]), x = NULL, colour = NULL, fill = NULL) +
  theme_minimal()
#load Zurich RSV data as sanity check to compare to MA data
library(data.table)
zrsv <- fread(file.choose())   # pick your Zurich RSV csv
names(zrsv)
head(zrsv)
str(zrsv)
zrsv[, .N, by = data_role]
zmeas[, .N, by = date][N > 1]     # empty = one per date; rows = replicates present
zmeas <- zrsv[data_role == "measurement" & !is.na(gc_per_mlww),
              .(date = as.Date(date), concentration = gc_per_mlww)][order(date)]
nrow(zmeas)
zmeas[, .N, by = date][N > 1]
zmeas[, replicate := rowid(date)]
ww_zurich <- sewer_data(measurements = zmeas,
                        flows = zmeas[, .(date, flow = 1e9)],
                        replicate_col = "replicate")

options(mc.cores = 4)
fit_zurich <- EpiSewer(data = ww_zurich, assumptions = assume_rsv,
                       fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                         iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42)))
fit_zurich$fitted$diagnostic_summary()
plot_R(fit_zurich)

# collapse replicates -> one mean concentration per date
zmeas <- zrsv[data_role == "measurement" & !is.na(gc_per_mlww),
              .(concentration = mean(gc_per_mlww)), by = .(date = as.Date(date))][order(date)]

# guard: confirm one row per date now
stopifnot(nrow(zmeas) == uniqueN(zmeas$date))

# flows: one row per date
zflows <- zmeas[, .(date, flow = 1e9)]

ww_zurich <- sewer_data(measurements = zmeas, flows = zflows)


options(mc.cores = 4)
fit_zurich <- EpiSewer(data = ww_zurich, assumptions = assume_rsv,
                       fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                         iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42)))
fit_zurich$fitted$diagnostic_summary()
plot_R(fit_zurich)
#plot florida, zurich, mass 
prep <- function(fit, lab) fit$summary$R[seeding == FALSE][, .(t = as.integer(date - min(date)), median, lower_0.95, upper_0.95, series = lab)]
all3 <- rbind(
  prep(fit_rsv,    "Florida 2024/25"),
  prep(fit_ma,     "Massachusetts 2024/25"),
  prep(fit_zurich, "Zurich 2023/24")
)
ggplot(all3, aes(t, median, colour = series, fill = series)) +
  geom_ribbon(aes(ymin = lower_0.95, ymax = upper_0.95), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) + geom_hline(yintercept = 1, linetype = 2) +
  labs(x = "Days since season start (~Aug 1)", y = expression(R[t]), colour = NULL, fill = NULL) +
  theme_minimal()
#checking to see if there are multiple ma sites
rsv[, season := ifelse(month(sample_collect_date) >= 8,
                       paste0(year(sample_collect_date), "/", year(sample_collect_date)+1),
                       paste0(year(sample_collect_date)-1, "/", year(sample_collect_date)))]
rsv[key_plot %in% c("fl_1702_wwtp","ma_2539_wwtp","ca_352_wwtp","fl_1631_wwtp"),
    .N, by = .(key_plot, season)][order(key_plot, season)]
#check to see which concentration methods are used over the season's 
rsv[key_plot == "fl_1702_wwtp",
    .(n = .N, methods = paste(unique(concentration_method), collapse = "; "),
      units = paste(unique(pcr_target_units), collapse = "; ")),
    by = season][order(season)]
#check flow for all sites 
rsv[key_plot %in% c("fl_1702_wwtp","fl_1631_wwtp","ca_352_wwtp","ma_2539_wwtp") &
      pcr_target_units == "copies/l wastewater",
    .(n = .N,
      flow_present = sum(!is.na(flow_rate)),
      flow_frac = round(mean(!is.na(flow_rate)), 2)),
    by = .(key_plot, season)][order(key_plot, season)]

#check to see what kind of pcr is used qpcr vs dpcr and then what model to apply
rsv[key_plot %in% c("fl_1702_wwtp", "fl_1631_wwtp", "ma_2539_wwtp", "ca_352_wwtp"),
    .N, by = .(key_plot, pcr_type)][order(key_plot)]
rsv[key_plot == "fl_1631_wwtp", .N, by = .(pcr_type, lab_id, concentration_method)]
#checking the full range of dates
# actual sample span + where signal actually exists, per site
measurements[, .(first = min(date), last = max(date), n = .N)]   # run for each site's measurements table
#sensitivity analysis for three sites 
cv <- 0.8
shed_means <- c(4.5, 6, 8, 10, 12, 13.5)

# the three fitted data objects you already built
site_data <- list("Florida" = ww_data, "Massachusetts" = ww_data_ma, "Zurich" = ww_zurich)

options(mc.cores = 4)
sweep3 <- rbindlist(lapply(names(site_data), function(site) {
  rbindlist(lapply(shed_means, function(m) {
    a <- sewer_assumptions(
      generation_dist    = gen_rsv,
      shedding_dist      = get_discrete_gamma(gamma_shape = 1/cv^2, gamma_scale = m * cv^2),
      shedding_reference = "infection")
    f <- EpiSewer(data = site_data[[site]], assumptions = a,
                  fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                    iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42)))
    r <- f$summary$R[seeding == FALSE]
    data.table(site = site, shedding_mean = m,
               peak_t   = r[which.max(median), as.integer(date - min(r$date))],
               peak_Rt  = r[which.max(median), median])
  }))
}))
sweep3
#Sensitivity analsis
cv <- 0.8
shed_means <- c(4.5, 6, 8, 10, 12, 13.5)
site_data <- list("Florida" = ww_data, "Massachusetts" = ww_data_ma, "Zurich" = ww_zurich)

options(mc.cores = 4)
sweep3 <- rbindlist(lapply(names(site_data), function(site) {
  rbindlist(lapply(shed_means, function(m) {
    tryCatch({
      a <- sewer_assumptions(generation_dist = gen_rsv,
                             shedding_dist = get_discrete_gamma(gamma_shape = 1/cv^2, gamma_scale = m * cv^2),
                             shedding_reference = "infection")
      f <- EpiSewer(data = site_data[[site]], assumptions = a,
                    fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                      iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42)))
      r <- f$summary$R[seeding == FALSE]
      data.table(site = site, shedding_mean = m,
                 peak_t  = r[which.max(median), as.integer(date - min(r$date))],
                 peak_Rt = r[which.max(median), median], status = "ok")
    }, error = function(e) data.table(site = site, shedding_mean = m,
                                      peak_t = NA_integer_, peak_Rt = NA_real_, status = paste("ERR:", conditionMessage(e))))
  }))
}))
sweep3
#Graph showing peaks
ggplot(sweep3, aes(shedding_mean, peak_t, colour = site)) +
  geom_line() + geom_point() +
  labs(x = "Assumed shedding mean (days)", y = "Day of peak Rt", colour = NULL) + theme_minimal()

ggplot(sweep3, aes(shedding_mean, peak_Rt, colour = site)) +
  geom_line() + geom_point() +
  labs(x = "Assumed shedding mean (days)", y = "Peak Rt", colour = NULL) + theme_minimal()

#graph showing peak relative to baseline
sweep3[, peak_t_shift := peak_t - peak_t[shedding_mean == 4.5], by = site]
ggplot(sweep3, aes(shedding_mean, peak_t_shift, colour = site)) +
  geom_line() + geom_point() +
  labs(x = "Assumed shedding mean (days)",
       y = "Shift in peak-Rt day (relative to 4.5 d)", colour = NULL) + theme_minimal()

sweep3[, peak_Rt_shift := peak_Rt - peak_Rt[shedding_mean == 4.5], by = site]
ggplot(sweep3, aes(shedding_mean, peak_Rt_shift, colour = site)) +
  geom_line() + geom_point() + geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Assumed shedding mean (days)",
       y = "Shift in peak Rt (relative to 4.5 d)", colour = NULL) + theme_minimal()
#check 4th site 
ca <- "ca_352_wwtp"

# units + matrix (must be liquid copies/L, not solids)
rsv[key_plot == ca, .N, by = .(pcr_target_units, sample_matrix, pcr_type)]

# method split + coverage per season
rsv[key_plot == ca,
    .(n = .N, first = min(sample_collect_date), last = max(sample_collect_date)),
    by = .(season, concentration_method)][order(season, -n)]

# flow: real or NA?
rsv[key_plot == ca, .(flow_present = sum(!is.na(flow_rate)), total = .N)]
#check california site 
one_ca <- rsv[key_plot == "ca_352_wwtp" &
                sample_collect_date >= as.Date("2024-08-01") &
                sample_collect_date <= as.Date("2025-06-30") &
                concentration_method == "none"]

meas_ca <- one_ca[, .(concentration = mean(pcr_target_avg_conc, na.rm = TRUE)),
                  by = .(date = as.Date(sample_collect_date))][order(date)]
stopifnot(nrow(meas_ca) == uniqueN(meas_ca$date))

flow_ca <- mean(one_ca$capacity_mgd, na.rm = TRUE) * 3785411.784
if (!is.finite(flow_ca)) flow_ca <- 1e9   # fallback if capacity is NA too
flows_ca <- meas_ca[, .(date, flow = flow_ca)]

ww_ca <- sewer_data(measurements = meas_ca, flows = flows_ca)

plot(meas_ca$date, meas_ca$concentration, type = "b",
     xlab = "Date", ylab = "RSV copies/L", main = "ca_352 2024-25")

options(mc.cores = 4)
fit_ca <- EpiSewer(data = ww_ca, assumptions = assume_rsv,
                   fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                     iter_warmup = 1000, iter_sampling = 1000, chains = 4, seed = 42)))
fit_ca$fitted$diagnostic_summary()
plot_R(fit_ca)
#check diagnostics
fit_ca$fitted$diagnostic_summary()   # num_divergent = 0? ebfmi > 0.3 on all chains?
#check to see if WAIC is in Episewer
# does EpiSewer expose log-likelihood or a built-in information criterion?
grep("waic|loo|log_lik|elpd", ls("package:EpiSewer"), value = TRUE, ignore.case = TRUE)
# and look at what the fit object carries
str(fit_rsv, max.level = 2)
#check loo 
loo_rsv <- fit_rsv$fitted$loo()
print(loo_rsv)
grep("loo|waic|elpd|compare|lik", ls("package:EpiSewer"), value = TRUE, ignore.case = TRUE)
#check to see how many people are served in my chosen catchments
rsv[key_plot %in% c("fl_1702_wwtp", "ma_2539_wwtp", "ca_352_wwtp"),
    .(population = unique(population_served),
      capacity_mgd = unique(capacity_mgd)),
    by = key_plot]
#check to see unique values
rsv[key_plot %in% c("fl_1702_wwtp", "ma_2539_wwtp", "ca_352_wwtp"),
    .N, by = .(key_plot, population_served)][order(key_plot, -N)]
#check 
rsv[key_plot %in% c("fl_1702_wwtp", "ma_2539_wwtp", "ca_352_wwtp"),
    .(population = as.numeric(names(sort(table(population_served), decreasing = TRUE))[1]),
      capacity_mgd = as.numeric(names(sort(table(capacity_mgd), decreasing = TRUE))[1])),
    by = key_plot]
#run MA site with more iterations due to low BFMI
options(mc.cores = 4)
fit_ma <- EpiSewer(data = ww_data_ma, assumptions = assume_rsv,
                   fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                     iter_warmup = 1000, iter_sampling = 1000, chains = 4, seed = 42)))
fit_ma$fitted$diagnostic_summary()
exists("fit_rsv"); exists("fit_ma"); exists("fit_ca"); exists("fit_zurich")
#four site overlay 
prep <- function(fit, lab) fit$summary$R[seeding == FALSE][, .(t = as.integer(date - min(date)), median, lower_0.95, upper_0.95, series = lab)]
all4 <- rbind(
  prep(fit_rsv,    "Florida 2024/25"),
  prep(fit_ma,     "Massachusetts 2024/25"),
  prep(fit_ca,     "California 2024/25"),
  prep(fit_zurich, "Zurich 2023/24")
)
ggplot(all4, aes(t, median, colour = series, fill = series)) +
  geom_ribbon(aes(ymin = lower_0.95, ymax = upper_0.95), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 1) + geom_hline(yintercept = 1, linetype = 2) +
  coord_cartesian(xlim = c(0, 220)) +
  labs(x = "Days since season start (~Aug 1)", y = expression(R[t]), colour = NULL, fill = NULL) +
  theme_minimal()
#table for all 4 sites
all4[, .(peak_t = t[which.max(median)], peak_Rt = round(max(median), 2)), by = series]]
#Double check florida season start
all4[t >= 30, .(peak_t = t[which.max(median)], peak_Rt = round(max(median), 2)), by = series]
#check crossing of 1 as optimal metric
crossing <- all4[order(series, t)][, {
  below <- which(median < 1)
  .(cross_t = if (length(below)) t[below[1]] else NA_integer_)
}, by = series]
crossing
grep("^plot", ls("package:EpiSewer"), value = TRUE)
#posterior check
library(ggplot2)
ppc_plot <- function(fit, observed, label) {
  pred <- fit$summary$concentration
  ggplot() +
    geom_ribbon(data = pred, aes(date, ymin = lower_0.95, ymax = upper_0.95), fill = "#4C72B0", alpha = 0.18) +
    geom_ribbon(data = pred, aes(date, ymin = lower_0.5,  ymax = upper_0.5),  fill = "#4C72B0", alpha = 0.30) +
    geom_line(data = pred, aes(date, median), colour = "#4C72B0", linewidth = 0.8) +
    geom_point(data = observed, aes(as.Date(date), concentration), size = 1, alpha = 0.55) +
    labs(x = NULL, y = "RSV concentration", title = paste0("Posterior predictive check — ", label)) +
    theme_minimal()
}
ppc_plot(fit_ca, meas_ca, "California")
#rescale graph
plot_concentration(fit_ca) + coord_cartesian(ylim = c(0, 1e5))
#check each site
plot_concentration(fit_rsv)     # Florida
plot_concentration(fit_ma)      # Massachusetts
plot_concentration(fit_ca)      # California
plot_concentration(fit_zurich)  # Zurich
#overlayed graph 
library(ggplot2); library(data.table)

ppc_all <- rbindlist(list(
  fit_rsv$summary$concentration[,    .(date, median, lower_0.95, upper_0.95, site = "Florida")],
  fit_ma$summary$concentration[,     .(date, median, lower_0.95, upper_0.95, site = "Massachusetts")],
  fit_ca$summary$concentration[,     .(date, median, lower_0.95, upper_0.95, site = "California")],
  fit_zurich$summary$concentration[, .(date, median, lower_0.95, upper_0.95, site = "Zurich")]
))

obs_all <- rbindlist(list(
  meas_fl[, .(date = as.Date(date), concentration, site = "Florida")],
  meas_ma[, .(date = as.Date(date), concentration, site = "Massachusetts")],
  meas_ca[, .(date = as.Date(date), concentration, site = "California")],
  zmeas[,   .(date = as.Date(date), concentration, site = "Zurich")]
))

ggplot() +
  geom_ribbon(data = ppc_all, aes(date, ymin = lower_0.95, ymax = upper_0.95), fill = "#4C72B0", alpha = 0.2) +
  geom_line(data = ppc_all, aes(date, median), colour = "#4C72B0", linewidth = 0.7) +
  geom_point(data = obs_all, aes(date, concentration), size = 0.8, alpha = 0.5) +
  facet_wrap(~ site, scales = "free", ncol = 2) +
  labs(x = NULL, y = "RSV concentration") + theme_minimal()
meas_fl <- one[, .(concentration = mean(pcr_target_avg_conc, na.rm = TRUE)),
               by = .(date = as.Date(sample_collect_date))][order(date)]

meas_ma <- one_ma[, .(concentration = mean(pcr_target_avg_conc, na.rm = TRUE)),
                  by = .(date = as.Date(sample_collect_date))][order(date)]
nrow(meas_fl); nrow(meas_ma)   # expect ~218 and ~84
#all 4 panels
ppc_all <- rbindlist(list(
  fit_rsv$summary$concentration[,    .(date, median, lower_0.95, upper_0.95, site = "Florida")],
  fit_ma$summary$concentration[,     .(date, median, lower_0.95, upper_0.95, site = "Massachusetts")],
  fit_ca$summary$concentration[,     .(date, median, lower_0.95, upper_0.95, site = "California")],
  fit_zurich$summary$concentration[, .(date, median, lower_0.95, upper_0.95, site = "Zurich")]
))

obs_all <- rbindlist(list(
  meas_fl[, .(date = as.Date(date), concentration, site = "Florida")],
  meas_ma[, .(date = as.Date(date), concentration, site = "Massachusetts")],
  meas_ca[, .(date = as.Date(date), concentration, site = "California")],
  zmeas[,   .(date = as.Date(date), concentration, site = "Zurich")]
))

ggplot() +
  geom_ribbon(data = ppc_all, aes(date, ymin = lower_0.95, ymax = upper_0.95), fill = "#4C72B0", alpha = 0.2) +
  geom_line(data = ppc_all, aes(date, median), colour = "#4C72B0", linewidth = 0.7) +
  geom_point(data = obs_all, aes(date, concentration), size = 0.8, alpha = 0.5) +
  facet_wrap(~ site, scales = "free", ncol = 2) +
  labs(x = NULL, y = "RSV concentration") + theme_minimal()
#check variability 
cover_check <- function(fit, obs, label) {
  pred <- fit$summary$concentration[, .(date = as.Date(date), lo = lower_0.95, hi = upper_0.95)]
  m <- merge(obs[, .(date = as.Date(date), concentration)], pred, by = "date")
  data.table(site = label, n = nrow(m),
             coverage = round(mean(m$concentration >= m$lo & m$concentration <= m$hi), 3))
}
rbind(
  cover_check(fit_rsv,    meas_fl, "Florida"),
  cover_check(fit_ma,     meas_ma, "Massachusetts"),
  cover_check(fit_ca,     meas_ca, "California"),
  cover_check(fit_zurich, zmeas,   "Zurich")
)
#check coverage splits
cover_split <- function(fit, obs, label) {
  pred <- fit$summary$concentration[, .(date = as.Date(date), lo = lower_0.95, hi = upper_0.95)]
  m <- merge(obs[, .(date = as.Date(date), concentration)], pred, by = "date")
  m[, detect := concentration > 0]
  m[, .(site = label, n = .N,
        coverage = round(mean(concentration >= lo & concentration <= hi), 3)), by = detect]
}
rbind(
  cover_split(fit_rsv,    meas_fl, "Florida"),
  cover_split(fit_ma,     meas_ma, "Massachusetts"),
  cover_split(fit_ca,     meas_ca, "California"),
  cover_split(fit_zurich, zmeas,   "Zurich")
)
#rt curves for all sites
library(ggplot2); library(data.table)

rt_all <- rbindlist(list(
  fit_rsv$summary$R[seeding == FALSE,    .(date, median, lower_0.95, upper_0.95, site = "Florida")],
  fit_ma$summary$R[seeding == FALSE,     .(date, median, lower_0.95, upper_0.95, site = "Massachusetts")],
  fit_ca$summary$R[seeding == FALSE,     .(date, median, lower_0.95, upper_0.95, site = "California")],
  fit_zurich$summary$R[seeding == FALSE, .(date, median, lower_0.95, upper_0.95, site = "Zurich")]
))

p_rt <- ggplot(rt_all, aes(date, median)) +
  geom_ribbon(aes(ymin = lower_0.95, ymax = upper_0.95), fill = "#4C72B0", alpha = 0.2) +
  geom_line(colour = "#4C72B0", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ site, scales = "free_x", ncol = 2) +
  labs(x = NULL, y = expression(R[t])) +
  theme_minimal()

p_rt

#run RT sensitivity check
# timing sensitivity (the relative-shift version)
sweep3[, peak_t_shift := peak_t - peak_t[shedding_mean == 4.5], by = site]
p_timing <- ggplot(sweep3, aes(shedding_mean, peak_t_shift, colour = site)) +
  geom_line() + geom_point() +
  labs(x = "Assumed shedding mean (days)",
       y = "Shift in peak-Rt day (relative to 4.5 d)", colour = NULL) + theme_minimal()
ggsave("~/Desktop/thesis_figs/fig3_sensitivity_timing.png", plot = p_timing, width = 8, height = 5, dpi = 300)

# magnitude sensitivity
sweep3[, peak_Rt_shift := peak_Rt - peak_Rt[shedding_mean == 4.5], by = site]
p_mag <- ggplot(sweep3, aes(shedding_mean, peak_Rt_shift, colour = site)) +
  geom_line() + geom_point() + geom_hline(yintercept = 0, linetype = 2) +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  labs(x = "Assumed shedding mean (days)",
       y = "Shift in peak Rt (relative to 4.5 d)", colour = NULL) + theme_minimal()
ggsave("~/Desktop/thesis_figs/fig4_sensitivity_magnitude.png", plot = p_mag, width = 8, height = 5, dpi = 300)
#print plots
p_timing
p_mag
p_rt
unique(sweep3$site)
site_data <- list("Florida" = ww_data, "Massachusetts" = ww_data_ma,
                  "California" = ww_ca, "Zurich" = ww_zurich)
cv <- 0.8
shed_means <- c(4.5, 6, 8, 10, 12, 13.5)

options(mc.cores = 4)
sweep3 <- rbindlist(lapply(names(site_data), function(site) {
  rbindlist(lapply(shed_means, function(m) {
    tryCatch({
      a <- sewer_assumptions(generation_dist = gen_rsv,
                             shedding_dist = get_discrete_gamma(gamma_shape = 1/cv^2, gamma_scale = m * cv^2),
                             shedding_reference = "infection")
      f <- EpiSewer(data = site_data[[site]], assumptions = a,
                    fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                      iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42)))
      r <- f$summary$R[seeding == FALSE]
      data.table(site = site, shedding_mean = m,
                 peak_t  = r[which.max(median), as.integer(date - min(r$date))],
                 peak_Rt = r[which.max(median), median], status = "ok")
    }, error = function(e) data.table(site = site, shedding_mean = m,
                                      peak_t = NA_integer_, peak_Rt = NA_real_, status = paste("ERR:", conditionMessage(e))))
  }))
}))
sweep3
#check estimable window
# where does Rt estimation actually begin at each site?
sapply(list(FL = fit_rsv, MA = fit_ma, CA = fit_ca, ZH = fit_zurich),
       function(f) as.character(min(f$summary$R[seeding == FALSE, date])))

# when does each site's concentration signal start rising?
# (compare with the PPC figure — look for the first sustained non-zero stretch)
#flow robust check
# rebuild MA measurements (same filter as the original fit)
meas_ma <- one_ma[, .(concentration = mean(pcr_target_avg_conc, na.rm = TRUE)),
                  by = .(date = as.Date(sample_collect_date))][order(date)]

# constant-flow version of the same data
ww_ma_const <- sewer_data(measurements = meas_ma,
                          flows = meas_ma[, .(date, flow = 1e9)])

options(mc.cores = 4)
fit_ma_const <- EpiSewer(data = ww_ma_const, assumptions = assume_rsv,
                         fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                           iter_warmup = 1000, iter_sampling = 1000, chains = 4, seed = 42)))

fit_ma_const$fitted$diagnostic_summary()
#flow overlay
flowcheck <- rbind(
  fit_ma$summary$R[seeding == FALSE][,       .(date, median, lower_0.95, upper_0.95, flow = "measured")],
  fit_ma_const$summary$R[seeding == FALSE][, .(date, median, lower_0.95, upper_0.95, flow = "constant")]
)

p_flow <- ggplot(flowcheck, aes(date, median, colour = flow, fill = flow)) +
  geom_ribbon(aes(ymin = lower_0.95, ymax = upper_0.95), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) + geom_hline(yintercept = 1, linetype = 2) +
  labs(x = NULL, y = expression(R[t]), colour = NULL, fill = NULL) + theme_minimal()
p_flow

m <- merge(fit_ma$summary$R[seeding == FALSE, .(date, real = median)],
           fit_ma_const$summary$R[seeding == FALSE, .(date, const = median)], by = "date")
m[, .(max_abs_diff = round(max(abs(real - const)), 3),
      median_abs_diff = round(median(abs(real - const)), 3))]
#check all 4 sites again
# refit each site at the two extremes and overlay full Rt curves
extremes <- c(4.5, 13.5)
curve_data <- rbindlist(lapply(names(site_data), function(s) {
  rbindlist(lapply(extremes, function(m) {
    a <- sewer_assumptions(generation_dist = gen_rsv,
                           shedding_dist = get_discrete_gamma(gamma_shape = 1/0.8^2, gamma_scale = m*0.8^2),
                           shedding_reference = "infection")
    f <- EpiSewer(data = site_data[[s]], assumptions = a,
                  fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                    iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42)))
    f$summary$R[seeding == FALSE][, .(t = as.integer(date - min(date)), median,
                                      site = s, shed = paste0(m, " d"))]
  }))
}))

ggplot(curve_data, aes(t, median, colour = factor(shed))) +
  geom_line(linewidth = 0.8) + geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ site, scales = "free_x") +
  labs(x = "Days since season start", y = expression(R[t]), colour = "Shedding mean") +
  theme_minimal()
#same plot with ribbons
extremes <- c(4.5, 13.5)
curve_data <- rbindlist(lapply(names(site_data), function(s) {
  rbindlist(lapply(extremes, function(m) {
    a <- sewer_assumptions(generation_dist = gen_rsv,
                           shedding_dist = get_discrete_gamma(gamma_shape = 1/0.8^2, gamma_scale = m*0.8^2),
                           shedding_reference = "infection")
    f <- EpiSewer(data = site_data[[s]], assumptions = a,
                  fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                    iter_warmup = 500, iter_sampling = 500, chains = 4, seed = 42)))
    f$summary$R[seeding == FALSE][, .(t = as.integer(date - min(date)),
                                      median, lower_0.95, upper_0.95,
                                      site = s, shed = paste0(m, " d"))]
  }))
}))

ggplot(curve_data, aes(t, median, colour = factor(shed), fill = factor(shed))) +
  geom_ribbon(aes(ymin = lower_0.95, ymax = upper_0.95), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(~ site, scales = "free_x") +
  labs(x = "Days since season start", y = expression(R[t]),
       colour = "Shedding mean", fill = "Shedding mean") +
  theme_minimal()
#cross site consitency
prep <- function(fit, lab) fit$summary$R[seeding == FALSE][, .(t = as.integer(date - min(date)), median, lower_0.95, upper_0.95, series = lab)]
all4 <- rbind(
  prep(fit_rsv,    "Florida"),
  prep(fit_ma,     "Massachusetts"),
  prep(fit_ca,     "California"),
  prep(fit_zurich, "Zurich")
)

p_crosssite <- ggplot(all4, aes(t, median, colour = series, fill = series)) +
  geom_ribbon(aes(ymin = lower_0.95, ymax = upper_0.95), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = 2) +
  coord_cartesian(xlim = c(0, 220)) +
  labs(x = "Days since season start (~Aug 1)", y = expression(R[t]),
       colour = NULL, fill = NULL) +
  theme_minimal()
p_crosssite

ggsave("~/Desktop/thesis_figs/fig_crosssite_overlay.png", plot = p_crosssite, width = 9, height = 6, dpi = 300)

#check to see how massachusetts would look with different estimable window
p_crosssite <- ggplot(all4, aes(t, median, colour = series, fill = series)) +
  geom_ribbon(aes(ymin = lower_0.95, ymax = upper_0.95), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = 2) +
  coord_cartesian(xlim = c(0, 300)) +
  labs(x = "Days since season start (Aug 1)", y = expression(R[t]),
       colour = NULL, fill = NULL) +
  theme_minimal()

#checking mass data
meas_ma[, .(first = min(date), last = max(date), n = .N)]
meas_ma[concentration > 0][order(date)][1:5]
meas_fl[concentration > 0][order(date)][1]
meas_ca[concentration > 0][order(date)][1]
#check flow max
m <- merge(fit_ma$summary$R[seeding == FALSE, .(date, real = median)],
           fit_ma_const$summary$R[seeding == FALSE, .(date, const = median)], by = "date")
m[, .(max_abs_diff = round(max(abs(real - const)), 3),
      median_abs_diff = round(median(abs(real - const)), 3))]
exists("ww_data"); exists("ww_data_ma"); exists("ww_ca"); exists("ww_zurich")
ls()
list.files("~/Desktop/thesis_figs")
# 1. what survived
ls()

# 2. Zurich CSV still on disk?
file.exists("/mnt/user-data/uploads/1784318285347_Fig_2c.csv")

# 3. current baseline mean
p <- assume_rsv$shedding_dist; sum(p*(seq_along(p)-1))/sum(p)
list.files("/mnt/user-data/uploads", full.names = TRUE)
ls()
list.files("~/Desktop", pattern = "\\.RData$", full.names = TRUE)
file.info(list.files("~/Desktop", pattern="\\.RData$", full.names=TRUE))[, "mtime", drop=FALSE]
load("/Users/cyruskohlmetz/Desktop/rsv_thesis_session.RData")
ls()
exists("fit_zurich"); exists("fit_ca"); exists("fit_ma")
fit_zurich$fitted$diagnostic_summary()$ebfmi   # should be ~0.78-0.86
fit_ca$fitted$diagnostic_summary()$ebfmi        # should be ~0.66-0.70
ls()  # find the sweep object name — likely "sweep"
str(sweep, max.level = 1)   # see its structure
head(sweep)                  # see what columns it has
ls()
fit_zurich$fitted$diagnostic_summary()
fit_ca$fitted$diagnostic_summary()
fit_rsv$fitted$diagnostic_summary()
fit_ma$fitted$diagnostic_summary()
#appendix building
# n and non-detects per site — adjust column names to your data
meas_fl[, .(n = .N, non_detect = sum(is.na(concentration) | concentration == 0))]
meas_ma[, .(n = .N, non_detect = sum(is.na(concentration) | concentration == 0))]
meas_ca[, .(n = .N, non_detect = sum(is.na(concentration) | concentration == 0))]
zmeas[, .(n = .N, non_detect = sum(is.na(concentration) | concentration == 0))]
save.image("/Users/cyruskohlmetz/Desktop/rsv_thesis_session.RData")
fo
fo_sweep <- function(d, a) EpiSewer(data = d, assumptions = a,
                                    fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(
                                      iter_warmup = 1000, iter_sampling = 1000, chains = 4, seed = 42)))
str(ww_ca, max.level = 1)
max(f_test$summary$R[seeding == FALSE]$median)
a_test <- sewer_assumptions(generation_dist = gen_rsv,
                            shedding_dist = get_discrete_gamma(gamma_shape = 1/0.8^2, gamma_scale = 4.6*0.8^2),
                            shedding_reference = "infection")
f_test <- fo_sweep(ww_ca, a_test)
max(f_test$summary$R[seeding == FALSE]$median)
library(EpiSewer)
library(data.table)
library(ggplot2)
a_test <- sewer_assumptions(generation_dist = gen_rsv,
                            shedding_dist = get_discrete_gamma(gamma_shape = 1/0.8^2, gamma_scale = 4.6*0.8^2),
                            shedding_reference = "infection")
f_test <- fo_sweep(ww_ca, a_test)
max(f_test$summary$R[seeding == FALSE]$median)
save.image
str(ww_data, max.level = 1)   # should show ~218 measurements for Florida
cv <- 0.8
shed_means <- c(4.5, 6, 8, 10, 12, 13.5)
sites <- list(Florida = ww_data, Massachusetts = ww_data_ma,
              California = ww_ca, Zurich = ww_zurich)

sweep_fits <- list()
sweep_results <- rbindlist(lapply(names(sites), function(sn) {
  rbindlist(lapply(shed_means, function(m) {
    a <- sewer_assumptions(
      generation_dist = gen_rsv,
      shedding_dist   = get_discrete_gamma(gamma_shape = 1/cv^2, gamma_scale = m*cv^2),
      shedding_reference = "infection")
    f <- fo_sweep(sites[[sn]], a)
    sweep_fits[[paste(sn, m, sep="_")]] <<- f
    r <- f$summary$R[seeding == FALSE]
    ebfmi <- f$fitted$diagnostic_summary()$ebfmi
    data.table(site = sn, shed_mean = m,
               peak_Rt   = round(max(r$median), 2),
               peak_date = r$date[which.max(r$median)],
               min_ebfmi = round(min(ebfmi), 2),
               divergent = sum(f$fitted$diagnostic_summary()$num_divergent))
  }))
}))
print(sweep_results)
save.image("/Users/cyruskohlmetz/Desktop/rsv_thesis_session.RData")
#appendix diagnostics
# richer diagnostics from the stored fits
sweep_diag_full <- rbindlist(lapply(names(sweep_fits), function(nm) {
  d <- sweep_fits[[nm]]$fitted$diagnostic_summary()
  data.table(
    fit = nm,
    divergent = sum(d$num_divergent),
    max_treedepth = sum(d$num_max_treedepth),
    min_ebfmi = round(min(d$ebfmi), 3),
    mean_ebfmi = round(mean(d$ebfmi), 3),
    n_chains_low_ebfmi = sum(d$ebfmi < 0.3)
  )
}))
print(sweep_diag_full)
# write to CSV to paste/import into your thesis
fwrite(sweep_diag_full, "/Users/cyruskohlmetz/Desktop/sweep_diagnostics.csv")
save.image("/Users/cyruskohlmetz/Desktop/rsv_thesis_session.RData")
load("/Users/cyruskohlmetz/Desktop/rsv_thesis_session.RData")
library(EpiSewer); library(data.table); library(ggplot2)
exists("sweep_fits")   # confirm the stored sweep fits are there
length(sweep_fits)     # should be 24
names(sweep_fits)      # e.g. "Florida_4.5", "Florida_6", ...
exists("sweep_fits")
length(sweep_fits)
names(sweep_fits[[1]]$summary$R)

#run appendix sweep
# Build trajectory table from all 24 sweep fits
sweep_traj <- rbindlist(lapply(names(sweep_fits), function(nm) {
  parts <- strsplit(nm, "_")[[1]]
  s <- parts[1]
  m <- as.numeric(parts[2])
  r <- sweep_fits[[nm]]$summary$R[seeding == FALSE]
  data.table(site = s, shed_mean = m, date = r$date, median = r$median)
}))

# One overlay plot per site (six shedding-mean trajectories)
dir.create("~/Desktop/thesis_figs", showWarnings = FALSE)
for (s in unique(sweep_traj$site)) {
  p <- ggplot(sweep_traj[site == s],
              aes(date, median, colour = factor(shed_mean), group = shed_mean)) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    labs(title = paste("Shedding-mean sensitivity —", s),
         y = expression(R[t]), x = NULL, colour = "Shedding\nmean (days)") +
    theme_minimal(base_size = 12)
  ggsave(paste0("~/Desktop/thesis_figs/sweep_", tolower(s), ".png"),
         p, width = 8, height = 5, dpi = 300)
}

list.files("~/Desktop/thesis_figs", pattern = "sweep_")
sweep_traj[site == "California" & date >= as.Date("2024-10-01") & date <= as.Date("2025-01-15"),
           .(interior_peak = round(max(median), 3)), by = shed_mean][order(shed_mean)]
sweep_traj[site == "California" & date >= as.Date("2024-10-01") & date <= as.Date("2025-01-15"),
           .(peak_date = date[which.max(median)]), by = shed_mean][order(shed_mean)]
# example — adjust column name and non-detect coding
nrow(meas_fl); nrow(meas_ma); nrow(meas_ca); nrow(zmeas)  # totals (you have these)
# non-detects: depends on how they're stored
names(meas_fl)
head(meas_fl)
# non-detect counts and proportions per site
rbind(
  data.table(site = "Florida",       n = nrow(meas_fl), non_detect = meas_fl[concentration == 0, .N]),
  data.table(site = "Massachusetts", n = nrow(meas_ma), non_detect = meas_ma[concentration == 0, .N]),
  data.table(site = "California",    n = nrow(meas_ca), non_detect = meas_ca[concentration == 0, .N]),
  data.table(site = "Zurich",        n = nrow(zmeas),   non_detect = zmeas[concentration == 0, .N])
)[, pct := round(100 * non_detect / n, 1)][]
table_a1 <- data.table(
  Site = c("Florida", "Massachusetts", "California", "Zurich"),
  Season = c("2024-25", "2024-25", "2024-25", "2023-24"),
  n = c(218, 84, 286, 334),
  Non_detects = c(108, 46, 148, 38),
  Non_detect_pct = c(49.5, 54.8, 51.7, 11.4),
  Conc_method = c("none", "ceres nanotrap", "none", "Eawag dPCR"),
  Units = c("copies/L", "copies/L", "copies/L", "gc/mL"),
  Flow = c("constant", "measured", "constant", "constant")
)
print(table_a1)

# markdown table (paste-friendly)
knitr::kable(table_a1, format = "pipe")

fwrite(table_a1, "/Users/cyruskohlmetz/Desktop/table_a1.csv")

diag_baseline <- rbindlist(lapply(
  list(Florida = fit_rsv, Massachusetts = fit_ma, California = fit_ca, Zurich = fit_zurich),
  function(f) {
    d <- f$fitted$diagnostic_summary()
    data.table(divergent = sum(d$num_divergent),
               max_treedepth = sum(d$num_max_treedepth),
               min_ebfmi = round(min(d$ebfmi), 3),
               max_ebfmi = round(max(d$ebfmi), 3))
  }), idcol = "Site")
print(diag_baseline)
fwrite(diag_baseline, "/Users/cyruskohlmetz/Desktop/diag_baseline.csv")

names(fit_zurich$summary)

names(fit_zurich$summary$concentration)
head(fit_zurich$summary$concentration)

coverage_fn <- function(fit, obs_data) {
  pred <- fit$summary$concentration[type == "estimate", .(date, lower_0.95, upper_0.95)]
  obs  <- obs_data[, .(date, observed = concentration)]
  m <- merge(obs, pred, by = "date")
  # detected samples only (observed > 0)
  det <- m[observed > 0]
  round(mean(det$observed >= det$lower_0.95 & det$observed <= det$upper_0.95), 3)
}

coverage_table <- data.table(
  Site = c("Florida","Massachusetts","California","Zurich"),
  Coverage = c(
    coverage_fn(fit_rsv,    meas_fl),
    coverage_fn(fit_ma,     meas_ma),
    coverage_fn(fit_ca,     meas_ca),
    coverage_fn(fit_zurich, zmeas)
  )
)
print(coverage_table)

str(ww_zurich, max.level = 1)
head(ww_zurich$flows)
# check if flow is constant (all the same value) or varying
length(unique(ww_zurich$flows$flow))   # or whatever the flow column is called

save.image("/Users/cyruskohlmetz/Desktop/rsv_thesis_session.RData")