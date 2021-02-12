
### Load packages ####

library(tidyverse)
library(lubridate)
library(broom)
library(Metrics)
library(MuMIn)
library(viridis)
library(scales)

#### Load and explore data ----
list.files()
cav <- read.csv("Monitoring_data.csv", head=TRUE, sep=";", dec = ".")
cav$Date <- mdy(cav$Date)

head(cav)
str(cav)
names(cav)
cav %>% dplyr::select(Date, Month, Year) %>% summary()
cav %>% distinct(Species) %>% nrow()
cav %>% distinct(Cave) %>% nrow()

#### Prepare data ----
## Exclude selected species - Presence/absence records or taxonomically problematic species
exclude <- c("aff Xyccarph sp.1", "Ochyrocera sp.1", "Pipa arrabali", "Scydmaenus sp.1", "Zelurus sp.", 
  "Carajas paraua", "Circoniscus buckupi", "Coarazuphium amazonicus", "Copelatus cessaima", 
  "Pseudonannolene sp.", "Turbellaria")

## Exclude species and make cave area histogram
cav_area <- cav %>% dplyr::select(Cave, Area, Species) %>% 
  dplyr::filter(! Species %in% exclude) %>% group_by(Cave) %>%
  summarise(Area = mean(Area))

area.plot <- cav_area %>% ggplot(aes(Area)) + geom_histogram(bins=50) + 
  geom_vline(xintercept=25, col="red") +
  ylab("Number of caves") + xlab("Area (m²)") + theme_bw() + 
  theme(
    axis.title.x = element_text(size=15, face="bold"),
    axis.title.y = element_text(size=15, face="bold"))

area.plot
ggsave("area.plot.png", plot = area.plot, device = "png", width = 6, height = 6, dpi = 300)

## Exclude species and nest data by cave
cav_nested <- cav %>% dplyr::select(Cave, Date, Season, Abundance, Species) %>% 
  dplyr::filter(! Species %in% exclude) %>%
  group_by(Cave) %>% nest()

## Calculate number of field surveys per cave and plot histogram
cav_nested <- cav_nested %>% mutate(Ncamp = map_dbl(data, function(.x){ 
    .x %>% distinct(Date, .keep_all= TRUE) %>% nrow()})) 

cav_nested %>% unnest(data)
table(cav_nested$Ncamp)
cav_nested %>%  ggplot(aes(Ncamp)) + geom_histogram(binwidth=1)

## Save file map file
Caves_map <- cav %>% dplyr::filter(! Species %in% exclude) %>%
  group_by(Cave) %>% nest() %>%
  mutate(Ncamp = map_dbl(data, function(.x){ 
    .x %>% distinct(Date, .keep_all= TRUE) %>% nrow()})) %>%
  unnest(data) %>%
  dplyr::select(c(Cave, Long, Lat, Ncamp)) %>%
  distinct()

write.csv(Caves_map, file="Caves_map.csv")

#### Select species 
species <- cav_nested %>% unnest(data) %>%
  group_by(Species) %>% summarize(n = n()) %>% 
  arrange(desc(n)) %>% print(n=33)

sp.names <- as.character(species$Species)

##### Plot abundance ----
cav_nested$data[[2]] %>% 
  ggplot(aes(x=Date, y=Abundance, color = Species)) + geom_line() 

cav_nested$data[[9]] %>% filter(Species=="Phalangopsis sp.1") %>% 
  #slice_sample(.data, n=2, replace = FALSE) %>%
  ggplot(aes(x=Date, y=Abundance)) + geom_point() + 
  geom_smooth(method="lm", se=FALSE, size=1, col="black") # glm.nb

cav_nested$data[[7]] %>% filter(Species=="Phalangopsis sp.1") %>%
  ggplot(aes(x=Date, y=Abundance, col=Season)) + geom_point() + 
  geom_smooth(method="lm", se=FALSE, size=1) # glm.nb

#### Test the effect of seasonality ----
#### Define seasonality.lm function
## data = data with a data list column, species = species name, dry.min = min sample size for the dry season, rain.min = min sample size for rainy season
seasonality.lm <- function(data, species, dry.min, rain.min) {
  
  ## Filter by speciesn split by season and filter by minimum sample size
  cav_nested_seas <- data %>% mutate(Species = map(data, function(.x){
    .x %>% filter(Species==species)})) %>%
    mutate(Dry = map(Species, function(.x){
      .x %>% filter(Season == "Seca") %>%
        group_by(Date) %>% summarize(Ntot = sum(Abundance), .groups = "drop")})) %>%
    mutate(N.Dry = map_dbl(Dry, ~nrow(.x))) %>% 
    mutate(Rain = map(Species, function(.x){
      .x %>% filter(Season == "Chuvosa") %>%
        group_by(Date) %>% summarize(Ntot = sum(Abundance), .groups = "drop")})) %>%
    mutate(N.Rain = map_dbl(Rain, ~nrow(.x))) %>%
    filter(N.Dry>=dry.min, N.Rain>=rain.min) %>%
    mutate(data.cc = map(Species, function(.x){ .x %>%
        group_by(Date, Season) %>% summarize(Ntot = sum(Abundance), .groups = "drop")}))
  
  ## Run lm if nrow!=0 and retrieve p-value for interaction effect
  if(nrow(cav_nested_seas)==0){
    return(NA)} else{
      
  res <- cav_nested_seas %>% mutate(LM = map(data.cc, ~ lm(formula = Ntot ~ Date*Season, data = .x))) %>%
  mutate(res = map(LM, ~tidy(.x))) %>%
  unnest(res) %>% filter(term=="Date:SeasonSeca") %>% 
  dplyr::select(c(Cave, p.value)) %>%
  mutate(species = species)

  return(res)
  }
}

#### Run loop for all species
res.seas <- list()
for(i in 1:length(sp.names)){
  res.seas[[i]] <- seasonality.lm(cav_nested, sp.names[i], 5, 5)
}

res.seas.full <- bind_rows(res.seas[is.na(res.seas) == FALSE]) %>% na.omit()
res.seas.full %>% nrow()
res.seas.full$BH.p.value <- p.adjust(res.seas.full$p.value, "BH")
#res.seas.full %>% arrange(BH.p.value)

Season_plot <- res.seas.full %>% 
  ggplot(aes(x=species, y=BH.p.value)) + geom_point(alpha=0.4, size=2) + 
  geom_hline(yintercept = 0.05, col="red", linetype = "dashed", size=0.5) +
  ylab("Adjusted p-value for interaction effect \n(Date:Season)") + 
  xlab("Taxa") + 
  theme_bw() + theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(size=15, face="bold"),
    axis.title.y = element_text(size=10, face="bold"))

Season_plot
ggsave("Season_plot.png", plot = Season_plot, device = "png", width = 6, height = 6, dpi = 300)

#### Test the effect of subsampling ----
#### Define functions to filter by sample size and fit lm for each season
## glm.nb is problematic with small sample sizes so using lm

#### coef.resample function
## data = dataset containing a data list column, n = number of resamples, Season = character indicating season
coef.resample <- function(data, n=2, Season){
  require(Metrics)
  slist <- list()
  for(i in 1:10){
    slist[[i]] <- data %>%
      mutate(Sp.sub = map(Date, function(.x){
        .x %>% slice_sample(n=n, replace=FALSE)})) %>%
      mutate(LM.sub = map(Sp.sub, ~ lm(formula = Ntot ~ Date, data = .x))) %>%
      mutate(res.sub = map(LM.sub, ~tidy(.x))) %>%
      unnest(res.sub) %>% filter(term=="Date") %>% 
      rename(estimate.sub = estimate) %>% 
      dplyr::select(Cave, estimate.full, estimate.sub)
  }
  
  res <- bind_rows(slist) %>% mutate(resample = n) %>%
    mutate(season = Season)
  return(res)
}

#### sp.filter.lm function
## data = dataset containing a Species list column, n.min = minimum number of surveys to be included, n.max= max number of surveys to be included
sp.filter.lm <- function(data, species){

  ## Filter dataset by one species and split by seasons
  cav_nested_sp <- data %>% mutate(Species = map(data, function(.x){
    .x %>% filter(Species==species)})) %>%
    mutate(Dry = map(Species, function(.x){
      .x %>% filter(Season == "Seca") %>%
        group_by(Date) %>% summarize(Ntot = sum(Abundance), .groups = "drop")})) %>%
    mutate(N.Dry = map_dbl(Dry, ~nrow(.x))) %>% 
    mutate(Rain = map(Species, function(.x){
      .x %>% filter(Season == "Chuvosa") %>%
        group_by(Date) %>% summarize(Ntot = sum(Abundance), .groups = "drop")})) %>%
    mutate(N.Rain = map_dbl(Rain, ~nrow(.x))) %>% 
    dplyr::select(-c(data, Ncamp))

  ## Filter data by at least 3 surveys and create dry subset  
  lm_sp_dry <- cav_nested_sp %>% 
    filter(N.Dry>=3) %>% dplyr::select(-contains("rain")) %>%
    mutate(Species2 = map(Species, function(.x){ .x %>%
    group_by(Date) %>% summarize(Ntot = sum(Abundance), .groups = "drop")})) %>%
    mutate(LM = map(Species2, ~ lm(formula = Ntot ~ Date, data = .x))) %>%
    mutate(res = map(LM, ~tidy(.x))) %>%
    unnest(res) %>% filter(term=="Date") %>% 
    dplyr::select(-c(Species, Species2, N.Dry, LM, term, std.error, statistic, p.value)) %>%
    rename(estimate.full = estimate, Date = Dry)

  ## Filter data by at least 3 surveys and create rain subset  
  lm_sp_rain <- cav_nested_sp %>% 
    filter(N.Rain>=3) %>% dplyr::select(-contains("dry")) %>%
    mutate(Species2 = map(Species, function(.x){ .x %>%
    group_by(Date) %>% summarize(Ntot = sum(Abundance), .groups = "drop")})) %>%
    mutate(LM = map(Species2, ~ lm(formula = Ntot ~ Date, data = .x))) %>%
    mutate(res = map(LM, ~tidy(.x))) %>%
    unnest(res) %>% filter(term=="Date") %>% 
    dplyr::select(-c(Species, Species2, N.Rain, LM, term, std.error, statistic, p.value)) %>%
    rename(estimate.full = estimate, Date = Rain)

  ## Define n.max
  n.max.dry <- max(unique(cav_nested_sp$N.Dry))
  n.max.rain <- max(unique(cav_nested_sp$N.Rain))
  
  ## Run coef.resample function in loops
  res.dry <- list()
  for(i in 2:n.max.dry){ 
    res.dry[[i]] <- coef.resample(lm_sp_dry, i, "Dry")
  }

  res.rain <- list()
  for(i in 2:n.max.rain){
    res.rain[[i]] <- coef.resample(lm_sp_rain, i, "Rain")
  }

  ## Merge dataframes and compute rmse
  sp_final <- bind_rows(res.dry, res.rain) %>% 
    group_by(resample, season) %>%
    summarise(rmse = rmse(estimate.full, estimate.sub), .groups = "drop") %>%
    mutate(species = species)

  ## Return results
  # print(table(cav_nested_sp$N.Dry))
  # print(table(cav_nested_sp$N.Rain))
  return(sp_final)
}

#### Run sp.filter.lm function for each species
length(sp.names)
test <- sp.filter.lm(data = cav_nested, species = sp.names[33])

rmse.results <- list()
for(i in c(1:33)){ 
  rmse.results[[sp.names[i]]] <- sp.filter.lm(data = cav_nested, species = sp.names[i])
}

length(rmse.results)
rmse.results[33]
names(rmse.results[33])
save(rmse.results, file = "rmse.results.RData")

#### Plot results
## Bind results
load("rmse.results.RData")
rmse_bind <- bind_rows(rmse.results)

show_col(viridis_pal()(20))
show_col(rainbow(10))
pal <- c("#FF0000FF", "#3300FFFF")

##
RMSE_plot <- rmse_bind  %>% ggplot(aes(x=resample, y = rmse, color=season)) + geom_line(size=2) +
  labs(x = "Number of surveys", y = "Root Mean Square Error", color = "Season") +
  facet_wrap(~ species,  scales = "free") +
  scale_color_manual(values=pal) +
  theme_bw() + theme(
    legend.direction = "horizontal",
    legend.position=c(0.8,0.05), 
    legend.title = element_text(size=20, face="bold"),
    legend.text = element_text(size=15, face="bold"),
    strip.text.x = element_text(size = 10, face="bold"),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"))

RMSE_plot

rmse_bind %>% filter(species == sp.names[10]) %>% ggplot(aes(x=resample, y = rmse, color=season)) + geom_line(size=2) +
  xlab("Number of surveys") + ylab("Root Mean Square Error") + 
  theme(legend.text = element_text(size=10, face="bold"))

ggsave("RMSE_plot.png", plot = RMSE_plot, device = "png", width = 15, height = 10, dpi = 300)

## Check min and max years
cav_nested %>% unnest(data) %>% filter(Season == "Chuvosa") %>% 
  dplyr::select(Date) %>% summary()
  
cav_nested %>% unnest(data) %>% filter(Season == "Seca") %>%
  dplyr::select(Date) %>% summary()

###### Identify indicator species ----
## Filter by species and by Nyears and run lm Abundance ~ Date
lm.per.sp <- function(data, species){
  data_lm <- data %>% 
    mutate(Species = map(data, function(.x){
      .x %>% filter(Species==species)})) %>%
    mutate(Nyears = map_dbl(Species, function(.x){ 
      .x %>% distinct(year(Date), .keep_all= TRUE) %>% nrow()})) %>%
    filter(Nyears >= 3) %>%
    mutate(data.cc = map(Species, function(.x){ .x %>%
        group_by(Date, Season) %>% summarize(Ntot = sum(Abundance), .groups = "drop")})) %>%
    mutate(LM = map(data.cc, ~ lm(formula = Ntot ~ Date, data = .x))) %>%
    mutate(res = map(LM, ~tidy(.x))) %>%
    unnest(res) %>% filter(term=="Date") %>% 
    dplyr::select(c(Cave, estimate)) %>%
    mutate(sp.name = species)
}

lm.res <- list()
for(i in 1:length(sp.names)){
  lm.res[[sp.names[i]]] <- lm.per.sp(cav_nested, sp.names[i])
}

lm_sp_cave <- bind_rows(lm.res)
hist(lm_sp_cave$estimate, breaks=50)

## Load landscape variables and append to table
list.files()
Lvars <- as_tibble(read.csv("Landscape_metrics.csv", sep = ",")) %>%
  dplyr::select(Cave, contains(c("Minning", "tdm")))
Lvars[, 2:ncol(Lvars)] <- scale(Lvars[, 2:ncol(Lvars)])
head(Lvars)
# Lvars_2019 <- Lvars %>% dplyr::select(c(cave1, contains('2019'), -contains('B500')))
# names(Lvars_2019)

lm_sp_cave_lan <- lm_sp_cave %>% inner_join(Lvars, by = "Cave") %>% 
  group_by(sp.name) %>% nest() %>%
  mutate(N = map_dbl(data, function(.x){ 
    .x %>% nrow()})) %>%
  filter(N >= 10)

lm_sp_cave_lan$data[[1]]
vars <- names(Lvars)[-1]

#### AICc
lm_sp_cave_lan_res <- lm_sp_cave_lan %>% mutate(sel.var = map(data, function(.x){
  res.lm <- data.frame(AICc = NA, var = NA, estimate = NA, p.val = NA)
  for(i in 1:length(vars)){
    formula <- paste("estimate ~", vars[i])
    lm <- lm(formula, data = .x)
    res.lm[i, 1] <- AICc(lm)
    res.lm[i, 2] <- vars[i]
    res.lm[i, 3] <- summary(lm)$coefficients[2,1]
    res.lm[i, 4] <- summary(lm)$coefficients[2,4]
  }
  return(res.lm[which.min(res.lm$AICc), 2:4])
})) %>% unnest(sel.var) %>% 
  filter(p.val <= 0.05) %>%
  dplyr::select(-c(data, N)) 

lm_sp_cave_lan_res
write.csv(lm_sp_cave_lan_res, file="lm_sp_cave_lan_res.csv", row.names=F)

#### LRT
lm_sp_cave_lan_res <- lm_sp_cave_lan %>% mutate(sel.var = map(data, function(.x){
  res.lrt <- data.frame(var = NA, lrt.p = NA, estimate = NA, p.val = NA)
  for(i in 1:length(vars)){
    formula <- paste("estimate ~", vars[i])
    lm <- lm(formula, data = .x)
    lrt.p <- drop1(lm, test="Chisq")[2, 5]
    res.lrt[i, 1] <- vars[i]
    res.lrt[i, 2] <- lrt.p
    res.lrt[i, 3] <- summary(lm)$coefficients[2,1]
    res.lrt[i, 4] <- summary(lm)$coefficients[2,4]
  }
  return(res.lrt %>% filter(lrt.p<=0.05))
})) %>% mutate(n.sel.var = map_dbl(sel.var, function(.x){.x %>% nrow()})) %>%
  dplyr::filter(n.sel.var>0) %>%
  dplyr::select(-c(N, sel.var, n.sel.var))

lm_sp_cave_lan_res

lm_sp_cave_lan_res <- lm_sp_cave_lan_res %>% mutate(LM = map(data, function(.x){
  res.tidy <- list()
  for(i in 1:length(vars)){
    formula <- paste("estimate ~", vars[i])
    res.tidy[[i]] <- tidy(lm(formula, data = .x)) 
  }
  return(bind_rows(res.tidy))}))

lm_sp_cave_lan_res$LM[[1]] 
lm_sp_cave_lan_res_plot <- lm_sp_cave_lan_res %>% unnest(LM) %>%
  dplyr::filter(term != "(Intercept)") %>% filter(p.value<=0.05) %>%
  mutate(metric = if_else(str_detect(term, "(Minning_)"), "Mining cover", "Distance to mine")) %>%
  mutate(year = if_else(str_detect(term, "(2015)"), "2015", 
                        if_else(str_detect(term, "(2016)"), "2016", 
                                if_else(str_detect(term, "(2017)"), "2017", 
                                        if_else(str_detect(term, "(2018)"), "2018", "2019"))))) %>%
  mutate(scale = if_else(str_detect(term, "B500"), "500", "1000")) %>%
  mutate(year = as.factor(year)) %>% mutate(scale = as.factor(scale))
  # mutate(term = recode(term, B500_Minning_pland_2015 = "Mining cover in 2015 (500m)",
  #                      B500_Minning_pland_2016 = "Mining cover in 2016 (500m)",
  #                      B500_Minning_pland_2017 = "Mining cover in 2017 (500m)",
  #                      B500_Minning_pland_2018 = "Mining cover in 2018 (500m)", 
  #                      B500_Minning_pland_2019 = "Mining cover in 2019 (500m)",
  #                      B1000_Minning_pland_2015 = "Mining cover in 2015 (1000m)", 
  #                      B1000_Minning_pland_2016 = "Mining cover in 2016 (1000m)", 
  #                      B1000_Minning_pland_2017 = "Mining cover in 2017 (1000m)", 
  #                      B1000_Minning_pland_2018 = "Mining cover in 2018 (1000m)", 
  #                      B1000_Minning_pland_2019 = "Mining cover in 2019 (1000m)", 
  #                      tdm_2015 = "Distance to mine in 2015", tdm_2016 = "Distance to mine in 2016",               
  #                      tdm_2017 = "Distance to mine in 2017", tdm_2018 = "Distance to mine in 2018", 
  #                      tdm_2019 = "Distance to mine in 2019"))
                           
lm_sp_cave_lan_res_plot

Coef_plot <- lm_sp_cave_lan_res_plot %>% 
  ggplot(aes(x = estimate, y = fct_reorder(sp.name, desc(sp.name)),
                          xmin = estimate - 2*std.error,
                          xmax = estimate + 2*std.error,
                          group = year, color = year, shape=scale)) +
  geom_pointrange(position = position_dodge(width = 0.5)) + 
  facet_wrap( ~ metric, scales = "free") +
  geom_vline(xintercept=0, col="black", linetype = "dashed", size=0.5) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = "Estimate", y = "Taxa", color = "Year", shape="Mining cover scale") +
  theme_bw() + theme(
    strip.text = element_text(size = 15, face="bold"),
    legend.position="bottom",
    legend.title = element_text(size=10, face="bold"),
    legend.text = element_text(size=10),
    axis.title.x = element_text(size=15, face="bold"),
    axis.title.y = element_text(size=15, face="bold"))
                          
Coef_plot
ggsave("Coeficients_plot.png", plot = Coef_plot, device = "png", width = 8, height = 5, dpi = 300)
