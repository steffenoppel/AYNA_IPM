# LIBRARIES #####
library(nimble)
library(coda)
library(tidyverse)
library(tidybayes)
library(strex)
library(beepr)
library(postpack)
library(here)
library(ggstream)
library(wesanderson)

# PLOT THEME ###

theme_murres <- function(){ 
  font <- "Helvetica"   #assign font family up front
  
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font, color = "black",           #set font family
        size = 14,                #set font size
        #face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font, color = "black",           #font family
        size = 14),               #font size
      
      axis.title = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 14),               #font size
      
      axis.ticks.x.bottom = element_blank(), 
      axis.ticks.y.left = element_blank(), 
      
      axis.text = element_text(              #axis text
        family = font, color = "black",           #axis famuly
        size = 12),                #font size
      
      legend.title = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 14),               #font size
      
      legend.text = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 12),              #font size
      
      strip.text = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 12),               #font size
      
      strip.background = element_blank()
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}

# LOAD SAMPLES #####
load("~/Documents/AYNA_IPM/samples_statespace_marginal_loaf_reduced_incolony_informed_COVARIATES_chain1.Rdata")
chain1 <- out1
load("~/Documents/AYNA_IPM/samples_statespace_marginal_loaf_reduced_incolony_informed_COVARIATES_chain2.Rdata")
chain2 <- out1
load("~/Documents/AYNA_IPM/samples_statespace_marginal_loaf_reduced_incolony_informed_COVARIATES_chain3.Rdata")
chain3 <- out1

out1 <- list(chain1 = chain1, chain2 = chain2, chain3 = chain3) %>% 
  as.mcmc.list()
nb <- dim(chain1)[1]/2
out1_wburnin <- lapply(out1, function(x) x[(nrow(x)-nb+1):nrow(x), ] %>% as.mcmc()) %>% 
  as.mcmc.list()
out1_wburnin_thinned <- post_thin(out1_wburnin, keep_iters = nb/10)
out1_wburnin_thinned <- post_subset(out1_wburnin_thinned, get_params(out1_wburnin_thinned, type = "base_index"),
                                    matrix = T) %>% 
  as.data.frame() %>% 
  filter(!is.nan(`IM[13, 1, 1]`) & !is.nan(`IM[13, 1, 3]`) & !is.nan(`Ntot[13]`)) %>% 
  slice(1:12000) %>% 
  bind_cols(ITER = rep(1:4000, times = 3), CHAIN = rep(1:3, each = 4000))
out1_wburnin_thinned <- post_convert(out1_wburnin_thinned %>% as.matrix())

summ <- t(post_summ(out1_wburnin_thinned, get_params(out1_wburnin_thinned, type = "base_index"), 
                    neff = TRUE, Rhat = TRUE , probs = c(0.025, 0.5, 0.975)
                    )) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name")

summ_noIM <- summ %>% 
  filter(!str_detect(name, "IM")) %>% 
  select(-neff)
beep(sound = 1)

write_csv(summ_noIM, file = here("Figures", "posterior_summary.csv"))

# Figures

## Methods ####

# study area map
  # DONE - in figures folder

# life cycle diagram
  # DONE - in figures folder

# mitigation effort plot

longline <- readRDS(here("longline.RDS"))
longline <- longline %>% 
  filter(Year >= 1985) %>% 
  select(-n_hooks) %>% 
  pivot_longer(-Year) %>% 
  mutate(name = case_when(
    str_detect(name, "ICCAT") ~ "ICCAT",
    str_detect(name, "NAM") ~ "Namibia",
    str_detect(name, "RSA") ~ "South Africa",
    TRUE ~ "Uruguay"
  )) %>% 
  arrange(Year, name)

pal <- wes_palette("Darjeeling2", 5, type = "continuous")

plot_mit <- ggplot(longline, aes(x = Year, y = value, color = name)) +
  geom_line(size = 1.25, alpha = 0.9) +
  theme_murres() + 
  ylab("Proportion of fleet using mitigation") + xlab("Year") +
  theme(legend.position = "top") + 
  guides(color=guide_legend(title="")) +
  scale_color_manual(values= pal)

png(here("Figures", "mitigation.png"), width = 6, height = 4, units = "in", res = 300)
plot_mit
dev.off()

## Results ####

# timeseries of survival by age

longline <- readRDS(here("longline.RDS")) %>% 
  filter(Year >= 1985)

w <- post_subset(out1_wburnin_thinned, "w\\[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() 

beta.mit <- post_subset(out1_wburnin_thinned, "beta.mit", matrix = T, iters = F, chains = F) %>% 
  as.data.frame()

phi.ad <- post_subset(out1_wburnin_thinned, "mean.phi.ad|sigma.phi|eps.phi", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  bind_cols(w, beta.mit = beta.mit$`beta.mit[2]`) %>% 
  mutate(mu.phi.ad = log(mean.phi.ad / (1-mean.phi.ad))) %>% 
  pivot_longer(1:36) %>% 
  mutate(year = str_first_number(name)) %>% 
  bind_cols(
    mit.ICCAT = rep(longline$mit.ICCAT[1:36], times = 15000),
    mit.NAM = rep(longline$mit.NAM[1:36], times = 15000),
    mit.RSA = rep(longline$mit.RSA[1:36], times = 15000),
    mit.URU = rep(longline$mit.URU[1:36], times = 15000)
  ) %>% 
  mutate(mit.index = `w[1]`*mit.ICCAT + `w[2]`*mit.NAM + `w[3]`*mit.RSA + `w[4]`*mit.URU) %>% 
  mutate(phi.ad = plogis(mu.phi.ad + value + beta.mit * mit.index)) %>% 
  mutate(year = year + 1985 - 1) %>% 
  select(year, phi.ad) %>% 
  rename(phi = phi.ad) %>% 
  mutate(Age = "Adult")

phi.im <- post_subset(out1_wburnin_thinned, "mean.phi.im|sigma.phi|eps.phi", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  bind_cols(w, beta.mit = beta.mit$`beta.mit[1]`) %>% 
  mutate(mu.phi.im = log(mean.phi.im / (1-mean.phi.im))) %>% 
  pivot_longer(1:36) %>% 
  mutate(year = str_first_number(name)) %>% 
  bind_cols(
    mit.ICCAT = rep(longline$mit.ICCAT[1:36], times = 15000),
    mit.NAM = rep(longline$mit.NAM[1:36], times = 15000),
    mit.RSA = rep(longline$mit.RSA[1:36], times = 15000),
    mit.URU = rep(longline$mit.URU[1:36], times = 15000)
  ) %>% 
  mutate(mit.index = `w[1]`*mit.ICCAT + `w[2]`*mit.NAM + `w[3]`*mit.RSA + `w[4]`*mit.URU) %>% 
  mutate(phi.im = plogis(mu.phi.im + value + beta.mit * mit.index)) %>% 
  mutate(year = year + 1985 - 1) %>% 
  select(year, phi.im) %>% 
  rename(phi = phi.im) %>% 
  mutate(Age = "Immature")

phi.juv <- post_subset(out1_wburnin_thinned, "mean.phi.juv|sigma.phi|eps.phi", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  bind_cols(w, beta.mit = beta.mit$`beta.mit[1]`) %>% 
  mutate(mu.phi.juv = log(mean.phi.juv / (1-mean.phi.juv))) %>% 
  pivot_longer(1:36) %>% 
  mutate(year = str_first_number(name)) %>% 
  bind_cols(
    mit.ICCAT = rep(longline$mit.ICCAT[1:36], times = 15000),
    mit.NAM = rep(longline$mit.NAM[1:36], times = 15000),
    mit.RSA = rep(longline$mit.RSA[1:36], times = 15000),
    mit.URU = rep(longline$mit.URU[1:36], times = 15000)
  ) %>% 
  mutate(mit.index = `w[1]`*mit.ICCAT + `w[2]`*mit.NAM + `w[3]`*mit.RSA + `w[4]`*mit.URU) %>% 
  mutate(phi.juv = plogis(mu.phi.juv + value + beta.mit * mit.index)) %>% 
  mutate(year = year + 1985 - 1) %>% 
  select(year, phi.juv) %>% 
  rename(phi = phi.juv) %>% 
  mutate(Age = "Juvenile")

toplot_survival <- bind_rows(phi.ad, phi.im, phi.juv) %>% 
  filter(Age != "Immature")
toplot_surv_means <- toplot_survival %>% 
  group_by(Age) %>% 
  summarise(mean = mean(phi))

pal <- wes_palette("Zissou1", 3, type = "continuous") 

plot_surv <- ggplot(toplot_survival, aes(x = year, y = phi, color = Age, fill = Age)) +
  stat_eye(alpha = 0.5) + 
  theme_murres() + 
  theme(legend.position = 'none') +
  ylab("Annual Survival Probability") + xlab("Year") +
  ylim(c(0.5,1)) +
  scale_fill_manual(values= pal) +
  scale_color_manual(values= pal) +
  facet_wrap(~ Age, nrow = 2) + 
  geom_hline(data = toplot_surv_means, aes(yintercept=mean, col = Age), alpha = 0.5) 
  
png(here("Figures", "annual_surv.png"), width = 6, height = 4, units = "in", res = 300)
plot_surv
dev.off()

# fecundity timeseries

fec_samps <- post_subset(out1_wburnin_thinned, "mean.fec|eps.fec", 
                              matrix = T) %>% 
  as.data.frame() %>% 
  mutate(mu.fec = log(mean.fec / (1-mean.fec))) %>% 
  mutate(across(1:13, ~ plogis(mu.fec + .))) %>% 
  select(1:13) %>% 
  pivot_longer(everything()) %>% 
  mutate(Year = str_first_number(name) + 2008 - 1)

toplot_fec <- fec_samps[, -1]

pal <- wes_palette("Darjeeling2", 5, type = "continuous") 

plot_fec <- ggplot(toplot_fec, aes(x = Year, y = value)) +
  stat_eye(color = pal[2], fill = pal[4]) + 
  geom_hline(yintercept=mean(toplot_fec$value), alpha = 0.5) + 
  theme_murres() + 
  ylab("Fecundity") + xlab("Year") +
  ylim(c(0,1))

png(here("Figures", "annual_fec.png"), width = 6, height = 4, units = "in", res = 300)
plot_fec
dev.off()

# recruitment curve

ages <- matrix(1:15, nrow = 15000, ncol = length(1:15), byrow = T) %>% 
  as.data.frame()
agebeta <- post_subset(out1_wburnin_thinned, "mu.p.juv|agebeta", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  bind_cols(ages) %>% 
  mutate(across(3:17, ~ plogis(mu.p.juv + agebeta * .)))

toplot_recruit <- agebeta %>% 
  pivot_longer(3:17) %>% 
  mutate(Age = str_first_number(name))

pal <- wes_palette("Chevalier1", 5, type = "continuous") 

plot_recruit <- ggplot(toplot_recruit, aes(x = Age, y = value)) +
  stat_eye(color = pal[1], fill = pal[1], alpha = 0.75) + 
  theme_murres() + 
  ylab("Recruitment Probability") + xlab("Age") +
  ylim(c(0,1)) 

png(here("Figures", "recruitment_prob.png"), width = 6, height = 4, units = "in", res = 300)
plot_recruit
dev.off()

# weights of each fleet on mitigation

toplot_w <- w %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(name = case_when(
    str_first_number(name) == 1 ~ "ICCAT",
    str_first_number(name) == 2 ~ "Namibia",
    str_first_number(name) == 3 ~ "South Africa",
    TRUE ~ "Uruguay"
  )) %>% 
  group_by(name) %>% 
  mutate(iter = row_number()) %>% 
  ungroup() %>% 
  arrange(iter, name)

pal <- wes_palette("Darjeeling2", 5, type = "continuous")

plot_w <- ggplot(toplot_w[1:1000, ], aes(x = iter, y = value, fill = name)) +
  geom_stream(type = "proportional", alpha = 0.75) +
  theme_murres() + 
  ylab("Weight") + xlab("") +
  theme(axis.text.x = element_blank(),, 
        legend.position = "top") + 
  guides(fill=guide_legend(title="")) +
  scale_fill_manual(values= pal)

png(here("Figures", "mitigation_weights.png"), width = 6, height = 4, units = "in", res = 300)
plot_w
dev.off()

# abundance timeseries, with population projections by scenario

load(here("PVA_results.RData"))

Ntot_samps <- post_subset(out1_wburnin_thinned, "Ntot\\[", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Year = str_first_number(name) + 2008 - 1)

PVA_Ntot <- apply(Ntot, 2, rbind) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(scenario = case_when(
    as.numeric(rowname) %% 3 == 1 ~ 1,
    as.numeric(rowname) %% 3 == 2 ~ 2,
    as.numeric(rowname) %% 3 == 0 ~ 3
  )) %>% 
  select(-rowname) %>% 
  pivot_longer(-scenario) %>% 
  mutate(name = str_first_number(name) + 2020) %>% 
  rename(Year = name) 

toplot_N <- Ntot_samps %>% 
  mutate(scenario = 1) %>% 
  select(scenario, Year, value) %>% 
  bind_rows(PVA_Ntot) %>% 
  group_by(Year, scenario) %>% 
  summarize(mean = mean(value), 
            lower = quantile(value, 0.025), 
            upper = quantile(value, 0.975)) 
  
pal <- wes_palette("Moonrise3", 5, type = "continuous")

plot_abund_total <- ggplot(toplot_N) +
  geom_ribbon(aes(x=Year,
                  ymin=lower,
                  ymax=upper,
                  fill=as.character(scenario)
                  ),
              alpha=0.35) +
  geom_line(aes(x=Year,
                y=mean,
                color = as.character(scenario)
                ),
            size = 0.75) +
  theme_murres() + 
  theme(legend.position = 'none') +
  ylab("Females in Study Area") + xlab("Year") +
  scale_fill_manual(values= pal) +
  scale_color_manual(values= pal)

png(here("Figures", "pop_size_projections.png"), width = 6, height = 4, units = "in", res = 300)
plot_abund_total
dev.off()

Ntot.breed_samps <- post_subset(out1_wburnin_thinned, "Ntot.breed", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Year = str_first_number(name) + 2008 - 1)

PVA_Ntot.breed <- apply(Ntot.breed, 2, rbind) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(scenario = case_when(
    as.numeric(rowname) %% 3 == 1 ~ 1,
    as.numeric(rowname) %% 3 == 2 ~ 2,
    as.numeric(rowname) %% 3 == 0 ~ 3
  )) %>% 
  select(-rowname) %>% 
  pivot_longer(-scenario) %>% 
  mutate(name = str_first_number(name) + 2020) %>% 
  rename(Year = name) 

toplot_N <- Ntot.breed_samps %>% 
  mutate(scenario = 1) %>% 
  select(scenario, Year, value) %>% 
  bind_rows(PVA_Ntot.breed) %>% 
  group_by(Year, scenario) %>% 
  summarize(mean = mean(value), 
            lower = quantile(value, 0.025), 
            upper = quantile(value, 0.975)) 

pal <- wes_palette("Moonrise3", 5, type = "continuous")

plot_abund_breed <- ggplot(toplot_N) +
  geom_ribbon(aes(x=Year,
                  ymin=lower,
                  ymax=upper,
                  fill=as.character(scenario)
  ),
  alpha=0.35) +
  geom_line(aes(x=Year,
                y=mean,
                color = as.character(scenario)
  ),
  size = 0.75) +
  theme_murres() + 
  theme(legend.position = 'none') +
  ylab("Breeding Pairs in Study Area") + xlab("Year") +
  scale_fill_manual(values= pal) +
  scale_color_manual(values= pal)

png(here("Figures", "breed_size_projections.png"), width = 6, height = 4, units = "in", res = 300)
plot_abund_breed
dev.off()

N.atsea_samps <- post_subset(out1_wburnin_thinned, "N.atsea", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Year = str_first_number(name) + 2008 - 1)

PVA_N.atsea <- apply(N.atsea, 2, rbind) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(scenario = case_when(
    as.numeric(rowname) %% 3 == 1 ~ 1,
    as.numeric(rowname) %% 3 == 2 ~ 2,
    as.numeric(rowname) %% 3 == 0 ~ 3
  )) %>% 
  select(-rowname) %>% 
  pivot_longer(-scenario) %>% 
  mutate(name = str_first_number(name) + 2020) %>% 
  rename(Year = name) 

toplot_N <- N.atsea_samps %>% 
  mutate(scenario = 1) %>% 
  select(scenario, Year, value) %>% 
  bind_rows(PVA_N.atsea) %>% 
  group_by(Year, scenario) %>% 
  summarize(mean = mean(value), 
            lower = quantile(value, 0.025), 
            upper = quantile(value, 0.975)) 

pal <- wes_palette("Moonrise3", 5, type = "continuous")

plot_abund_atsea <- ggplot(toplot_N) +
  geom_ribbon(aes(x=Year,
                  ymin=lower,
                  ymax=upper,
                  fill=as.character(scenario)
  ),
  alpha=0.35) +
  geom_line(aes(x=Year,
                y=mean,
                color = as.character(scenario)
  ),
  size = 0.75) +
  theme_murres() + 
  theme(legend.position = 'none') +
  ylab("Females at Sea") + xlab("Year") +
  scale_fill_manual(values= pal) +
  scale_color_manual(values= pal)

png(here("Figures", "sea_size_projections.png"), width = 6, height = 4, units = "in", res = 300)
plot_abund_atsea
dev.off()

N.loaf_samps <- post_subset(out1_wburnin_thinned, "N.loaf", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Year = str_first_number(name) + 2008 - 1)

PVA_N.loaf <- apply(N.loaf, 2, rbind) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(scenario = case_when(
    as.numeric(rowname) %% 3 == 1 ~ 1,
    as.numeric(rowname) %% 3 == 2 ~ 2,
    as.numeric(rowname) %% 3 == 0 ~ 3
  ))  %>% 
  select(-rowname) %>% 
  pivot_longer(-scenario) %>% 
  mutate(name = str_first_number(name) + 2020) %>% 
  rename(Year = name) 

toplot_N <- N.loaf_samps %>% 
  mutate(scenario = 1) %>% 
  select(scenario, Year, value) %>% 
  bind_rows(PVA_N.loaf) %>% 
  group_by(Year, scenario) %>% 
  summarize(mean = mean(value), 
            lower = quantile(value, 0.025), 
            upper = quantile(value, 0.975)) 

pal <- wes_palette("Moonrise3", 5, type = "continuous")

plot_abund_loaf <- ggplot(toplot_N) +
  geom_ribbon(aes(x=Year,
                  ymin=lower,
                  ymax=upper,
                  fill=as.character(scenario)
  ),
  alpha=0.35) +
  geom_line(aes(x=Year,
                y=mean,
                color = as.character(scenario)
  ),
  size = 0.75) +
  theme_murres() + 
  theme(legend.position = 'none') +
  ylab("Females Loafing") + xlab("Year") +
  scale_fill_manual(values= pal) +
  scale_color_manual(values= pal)

png(here("Figures", "loaf_size_projections.png"), width = 6, height = 4, units = "in", res = 300)
plot_abund_loaf
dev.off()

JUV_samps <- post_subset(out1_wburnin_thinned, "JUV", matrix = T, iters = F, chains = F) %>% 
  as.data.frame() %>% 
  pivot_longer(everything()) %>% 
  mutate(Year = str_first_number(name) + 2008 - 1)

PVA_JUV <- apply(JUV, 2, rbind) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(scenario = case_when(
    as.numeric(rowname) %% 3 == 1 ~ 1,
    as.numeric(rowname) %% 3 == 2 ~ 2,
    as.numeric(rowname) %% 3 == 0 ~ 3
  ))  %>% 
  select(-rowname) %>% 
  pivot_longer(-scenario) %>% 
  mutate(name = str_first_number(name) + 2020) %>% 
  rename(Year = name) 

toplot_N <- JUV_samps %>% 
  mutate(scenario = 1) %>% 
  select(scenario, Year, value) %>% 
  bind_rows(PVA_JUV) %>% 
  group_by(Year, scenario) %>% 
  summarize(mean = mean(value), 
            lower = quantile(value, 0.025), 
            upper = quantile(value, 0.975)) 

pal <- wes_palette("Moonrise3", 5, type = "continuous")

plot_abund_juv <- ggplot(toplot_N) +
  geom_ribbon(aes(x=Year,
                  ymin=lower,
                  ymax=upper,
                  fill=as.character(scenario)
  ),
  alpha=0.35) +
  geom_line(aes(x=Year,
                y=mean,
                color = as.character(scenario)
  ),
  size = 0.75) +
  theme_murres() + 
  theme(legend.position = 'none') +
  ylab("Juvenile Females") + xlab("Year") +
  scale_fill_manual(values= pal) +
  scale_color_manual(values= pal)

png(here("Figures", "juv_size_projections.png"), width = 6, height = 4, units = "in", res = 300)
plot_abund_juv
dev.off()



