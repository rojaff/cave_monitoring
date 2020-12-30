### Load packages ####

library(tidyverse)
library(lubridate)
library(vegan)
library(ade4)
library(mvabund)
library(usdm)

### Load data ####

# CKS monitoring dataset
list.files()
cav <- read.csv("Monitoring_data.csv", head=TRUE, sep=",", dec = ".")
cav$Date <- mdy(cav$Date)
cav<-as_tibble(cav)


head(cav)
str(cav)
cav$Location <- as.factor(cav$Location)
cav$Year <- as.factor(cav$Year)
cav$Abundance <- as.numeric(as.character(cav$Abundance))
cav$Area <- as.numeric(as.character(cav$Area))
names(cav)


## Exclude selected species - Presence/absence records or taxonomically problematic species
exclude <- c("aff Xyccarph sp.1", "Ochyrocera sp.1", "Pipa arrabali", "Scydmaenus sp.1", "Zelurus sp.", 
             "Carajas paraua", "Circoniscus buckupi", "Coarazuphium amazonicus", "Copelatus cessaima", 
             "Pseudonannolene sp.", "Turbellaria")

cav <- cav %>%  
  dplyr::filter(! Species %in% exclude)

# Landscape data

landscape_metrics <- read.csv("Landscape_metrics.csv", head=TRUE, sep=",", dec = ".")
landscape_metrics <- as_tibble(landscape_metrics)

landscape_metrics2 <- landscape_metrics %>%
  dplyr::select("Cave","B500_Forest_pland_2015":"B1000_Allvegetation_pland_2019") %>%
  pivot_longer(!Cave, names_to = "metric", values_to = "value") %>%
  separate(metric, c("Scale", "Class", "metric", "Year"))

topodist_mine2 <- landscape_metrics %>%
  dplyr::select("Cave","tdm_2015":"tdm_2019") %>%
  pivot_longer(!Cave, names_to = "metric", values_to = "value") %>%
  separate(metric, c("metric", "Year"))



### Inspecting data ####

## Evaluate sampling tables
cav_sampling <- cav %>%
  dplyr::select(Cave, Date, Year, Season) %>%
  group_by(Cave, Year, Season) %>%
  distinct() %>%
  summarise(Ncamp=n()) %>%
  pivot_wider(names_from = c(Year,Season), values_from = Ncamp)%>%
  dplyr::select( "Cave","2015_Chuvosa","2015_Seca","2016_Chuvosa":"2019_Seca")

# Note there is no sampling in whole Serra Sul during 2015, only SN

# We can also check individual species records by year/season

cav_sampling2 <- cav %>%
  dplyr::select(Cave, Month, Year, Season, Species) %>%
  group_by(Cave, Year, Season, Month) %>%
  count(N=Species)


## Nest data by cave
cav_nested <- cav %>% dplyr::select(Cave, Date, Month, Year, Season, Abundance, Species) %>% 
  group_by(Cave) %>% 
  nest()


## Calculate number of field surveys per cave and plot histogram
cav_nested <- cav_nested %>% mutate(Ncamp = map_dbl(data, function(.x){ 
  .x %>% distinct(Date, .keep_all= TRUE) %>% nrow()})) 

cav_nested %>%  ggplot(aes(Ncamp)) + geom_histogram(binwidth=1)

# Note by comparing this histogram with the cave_sampling table that cave with 5 or less
# survey occasions were those from SSul not sampled in 2015 and 2019
# SSul has no surveys in 2015

### Arrange data matrices for analysis (Decisions taken) ####

# Considered samples from 2015, but then filtered only caves sampled in all remaining years

cav_comm <- as.data.frame(cav %>%
                            dplyr::select(Cave, Location, Date,Year, Season, Abundance, Species)%>%
                            group_by(Cave) %>% 
                            nest() %>% 
                            mutate(Ncamp = map_dbl(data, function(.x){ .x %>% 
                                distinct(Date, .keep_all= TRUE) %>% 
                                nrow()})) %>%
                            filter(Ncamp >= 8) %>%
                            unnest(data) %>% 
  dplyr::select(-Ncamp, -Date) %>% 
  arrange(Cave, Year, Season) %>%
  dplyr::group_by(Cave,Year, Season) %>%
  tidyr::pivot_wider(names_from = Species, values_from =Abundance, values_fn = list(Abundance = mean)) %>%
  mutate(across(c("Mesabolivar spp.":"Chelodesmidae sp."), as.integer)) %>%
  mutate_all(replace_na, 0))

str(cav_comm)

# Generate Taxonomic abundance-matrix 

cav_comm_df <- cav_comm[,-c(1:4)]
str(cav_comm_df)


cav_comm_df %>%
  rowSums() # All caves have records in all years


# Landscape predictors matrix

cav_land <- as.data.frame(cav_comm %>%
                            dplyr::select(Cave, Location, Year, Season) %>%
                            group_by(Cave, Location, Year, Season) %>%
                            distinct() %>%
                            inner_join(y=landscape_metrics2) %>%
                            pivot_wider(names_from = c("Scale", "Class", "metric"), values_from = "value") %>%
                            dplyr::select(Cave, Location, Year, Season, B500_Forest_pland,  B1000_Forest_pland, B500_Minning_pland, B1000_Minning_pland,
                                          B500_Canga_pland, B1000_Canga_pland) %>%
                            inner_join(y=topodist_mine2) %>%
                            pivot_wider(names_from = "metric", values_from = "value") %>%
                            arrange(Cave, Location, Year, Season))


cav_land_df <- cav_land[,-c(1:4)]


# Environmental predictors matrix

names(cav)

cav_env <- cav %>%
  dplyr::select(Cave, Location, Year, Season, Precipitation, Temperature, Area) %>%
  dplyr::group_by(Cave,Location, Year, Season) %>%
  summarise(Precipitation=mean(Precipitation),
            Temperature=mean(Temperature), Area = mean(Area)) %>%
  arrange(Cave, Location, Year, Season) %>%
  drop_na()



# Some environmental data is missing in some caves, which where dropped. We have to 
# adjust the other data matrix

# First let´s arrange the environmental matrix to match community 
cav_DF<- as.data.frame(cav_env %>%
                         inner_join(y=cav_comm) %>%
                         arrange(Cave, Location, Year, Season)) 



# And now add and correct landscape data 
cav_DF <- as.data.frame(cav_land %>%
                          inner_join(y=cav_DF) %>%
                          arrange(Cave, Location, Year, Season))

# Finally, order Year data and transform it, together with season, in a factor

cav_DF$Year <- factor(cav_DF$Year, order = TRUE, levels = c("2015", "2016", "2017", "2018", "2019"))

cav_DF$Season <- factor(cav_DF$Season, order = F, levels = c("Chuvosa", "Seca"))

str(cav_DF)


# Remove species with < 10 records

cav_DF_pa <- cav_DF[,c(15:47)]

cav_DF_pa[cav_DF_pa > 0] <- 1

cav_DF_pa <- bind_cols(tibble(Cave=cav_DF$Cave),tibble(cav_DF_pa))


Cave_registers <- cav_DF_pa %>%
  group_by(Cave) %>% 
  summarise_all(funs(sum))

Cave_registers_pa <- Cave_registers[,c(2:34)]

Cave_registers_pa[Cave_registers_pa > 0] <- 1

Cave_registers_pa <- bind_cols(tibble(Cave=Cave_registers$Cave),tibble(Cave_registers_pa))

caves <- as.data.frame(Cave_registers_pa %>%
                         dplyr::select("Mesabolivar spp.":"Chelodesmidae sp.") %>%
                         colSums) 


retirar <- c("Dytiscidae sp.1", "Glomeridesmus spelaeus", "Circoniscus sp.", "Chelodesmidae sp.",
             "Pseudonannolenidae sp.")
            


cav_DF2 <- cav_DF %>%
  dplyr::select(-(retirar))




### Overall constrained analysis (RDA)- Redundancy Analysis ####
# What contributes most to the patterns in species composition found: space, time or environment?

# Generate time variable to be used in the following analysis

# First, creating a variable that nests Season in Years
cav_DF3 <- as.data.frame(cav_DF2 %>%
                           unite("time", c("Year","Season"), sep = "_"))

cav_DF3$time <- factor(cav_DF3$time, order = TRUE, 
                       levels = c("2015_Chuvosa","2015_Seca","2016_Chuvosa","2016_Seca","2017_Chuvosa","2017_Seca","2018_Chuvosa","2018_Seca","2019_Chuvosa","2019_Seca"))


Season <- dudi.mix(cav_DF2$Season, scannf=F, nf=2)$tab
Year <- dudi.mix(cav_DF2$Year, scannf=F, nf=2)$tab
time <- dudi.mix(cav_DF3$time, scannf=F, nf=1)$tab #interaction of Year and Season




## Pre-RDA

# Hellinger transformation for abundance data (some outliers in abundance)
cav_comm_df.hel <- decostand(cav_DF3[14:41], "hellinger") 

## VIFS - checking for colinearity

cor(cav_DF3[,c(4:13)])

vif(cav_DF3[,c(4:10)]) # We´ll retain only the 1000 scale variables
cav_land_df <- cav_DF3[,c(5,7,9,10)]
vif(cav_land_df) #ok


vif(cav_DF3[,c(11:13)]) # ok
cav_env_df <- cav_DF3[,c(11:13)]


## Final preparation of predictor datasets to run RDA

Location <- cav_DF3$Location

alldata <- -cbind(Year, Season,cav_land_df,cav_env_df)
#or 
alldata2 <- -cbind(time,cav_land_df,cav_env_df) #One variable  nesting Year:Season


# Full RDA

# check between sets of variables
rda.comp.full<-rda(cav_comm_df.hel~.,data=alldata2,add=TRUE)

var<-varpart(cav_comm_df.hel,time, cav_land_df,cav_env_df)

plot(var, 
     digits = 2, 
     bg = c("red", "blue", "green"), 
     Xnames = c("Temporal", "Landscape", "Environmental"), 
     id.size = 0.7
)


# pRDA - Partial Redundancy analysis

rda.comp.p <- rda(cav_comm_df.hel~time$df.L+cav_land_df$B1000_Canga_pland+cav_land_df$B1000_Forest_pland+cav_land_df$B1000_Minning_pland+cav_land_df$tdm+
                    cav_env_df$Area+cav_env_df$Temperature+cav_env_df$Precipitation+Condition(Location),add=TRUE)

anova(rda.comp.p,permutations = how(nperm = 999), by="margin") # acess p-values

RsquareAdj(rda.comp.p) # acess adjusted R-squared


### Modelling multivariate composition data using GLMs (mvAbund approach) ####


# generate multivariate response variable data for mvAbund package

cav_spp2 <- mvabund(cav_DF3[,14:41])

meanvar.plot(cav_spp2) # clearly not a good sign for traditional ordination techniques  


# Does disturbance variables influence multi-species abundance patterns? 

# Multi-model selection using AIC, evaluating between time and area variables in
# a model in which disturbance variables are fixed (selected from the results of the RDA) 

N.Model_ta <- manyglm(cav_spp2~cav_DF3$time+cav_DF3$Area+cav_DF3$B1000_Minning_pland+cav_DF3$tdm, family="negative_binomial")

N.Model_a <-manyglm(cav_spp2~cav_DF3$Area+cav_DF3$B1000_Minning_pland+cav_DF3$tdm, family="negative_binomial")

N.Model_t <- manyglm(cav_spp2~cav_DF3$time+cav_DF3$B1000_Minning_pland+cav_DF3$tdm, family="negative_binomial")

N.Model_ta$AICsum # best model, including time and area
N.Model_a$AICsum
N.Model_t$AICsum


# extract values
output_coef <- N.Model_ta$coefficients

output_coef2 <- N.Model_ta$stderr.coefficients

term  <- dimnames(N.Model_ta$coefficients)[[1]]


# wrangle them 
# coefficient estimate
output_coef <- cbind(term,output_coef)

output_coef <- as_tibble(output_coef)

output_coef <- output_coef %>%
  pivot_longer(-term,names_to = "sp.name", values_to = "estimate") 

# standard error
output_coef2 <- cbind(term,output_coef2)

output_coef2 <- as_tibble(output_coef2)

output_coef2 <- output_coef2 %>%
  pivot_longer(-term,names_to = "sp.name", values_to = "std.error") 

# join both

output_coefficients <- as.data.frame(output_coef %>%
                                       inner_join(y=output_coef2))

str(output_coefficients)

output_coefficients$estimate <- as.numeric(output_coefficients$estimate)
output_coefficients$std.error <- as.numeric(output_coefficients$std.error)


output_coefficients_f <- output_coefficients %>%
  filter(term=="cav_DF3$tdm" | term=="cav_DF3$B1000_Minning_pland") %>%
  mutate(term = case_when(term == "cav_DF3$B1000_Minning_pland" ~ "Mining cover",
                          term == "cav_DF3$tdm" ~ "Distance to mine",
                          TRUE ~ "Other")) %>%
  filter(sp.name=="Charinus.ferreos" & term == "Mining cover"| sp.name=="Pyrgodesmidae.sp.1" & term == "Mining cover" |
           sp.name=="Escadabiidae.sp.1" & term == "Mining cover" | sp.name=="Sphendononema.guildingii"  & term == "Mining cover" |
           sp.name=="Stygnidae.sp.1" & term == "Mining cover" | sp.name=="Thecadactylus.rapicauda" & term == "Mining cover" |
           sp.name=="Pristimantis.fenestratus" & term == "Mining cover"| 
           sp.name=="Prodidomidae.sp."  & term == "Mining cover" | sp.name=="Phalangopsis.sp.1"  & term == "Mining cover" |
           sp.name=="Theraphosidae" & term == "Mining cover" |
           sp.name=="Astieae.sp.1" & term == "Distance to mine" |
           sp.name=="Spirostreptida.sp." & term == "Distance to mine"| sp.name=="Pristimantis.fenestratus" & term == "Distance to mine"|
           sp.name=="Plato.spp." & term == "Distance to mine"| sp.name=="Latebraria.sp." & term == "Distance to mine"|
           sp.name=="Leptodactylus.pentadactylus" & term == "Distance to mine"| sp.name=="Prodidomidae.sp."  & term == "Distance to mine" |
           sp.name=="Prionostemma.sp." & term == "Distance to mine"| sp.name=="Protimesius.sp." & term == "Distance to mine"|
           sp.name=="Theraphosidae" & term == "Distance to mine" | sp.name=="Cydninae.sp.1"  & term == "Distance to mine" |
           sp.name=="Escadabiidae.sp.1" & term == "Distance to mine" | sp.name=="Escadabiidae.sp.2" & term == "Distance to mine" |
           sp.name=="Uvaroviella.sp." & term == "Distance to mine") %>%
  mutate(sp.name = case_when(sp.name=="Charinus.ferreos" ~ "Charinus ferreus",
                             sp.name=="Astieae.sp.1" ~ "Astieae sp.1",
                             sp.name=="Spirostreptida.sp." ~ "Spirostreptida sp.",
                             sp.name=="Escadabiidae.sp.1" ~ "Escadabiidae sp.1",
                             sp.name=="Escadabiidae.sp.2" ~ "Escadabiidae sp.2",
                             sp.name=="Cydninae.sp.1" ~ "Cydninae sp.1",
                             sp.name=="Stygnidae.sp.1" ~ "Stygnidae sp.1",
                             sp.name=="Uvaroviella.sp." ~  "Uvaroviella sp.",
                             sp.name=="Pristimantis.fenestratus" ~  "Pristimantis fenestratus",
                             sp.name=="Prionostemma.sp." ~ "Prionostemma sp.",
                             sp.name=="Protimesius.sp." ~ "Protimesius sp.",
                             sp.name=="Prodidomidae.sp." ~ "Prodidomidae sp.",
                             sp.name=="Phalangopsis.sp.1" ~ "Phalangopsis sp.1",
                             sp.name=="Latebraria.sp." ~ "Latebraria sp.",
                             sp.name=="Plato.spp." ~ "Plato sp.",
                             sp.name=="Sphendononema.guildingii" ~ "Sphendononema guildingii",
                             sp.name=="Heterophrynus.longicornis" ~ "Heterophrynus longicornis", sp.name=="Theraphosidae" ~ "Theraphosidae",
                             sp.name=="Leptodactylus.pentadactylus" ~ "Leptodactylus pentadactylus",
                             sp.name=="Pyrgodesmidae.sp.1" ~ "Pyrgodesmidae sp.1", sp.name=="Spirostreptida.sp." ~ "Spirostreptida sp.",
                             sp.name=="Thecadactylus.rapicauda" ~ "Thecadactylus rapicauda")) %>%
  mutate(sp.name = fct_reorder(sp.name,desc(sp.name)))



Coef_plot <- ggplot(output_coefficients_f,aes(x = estimate, y = sp.name,
                                              xmin = estimate - 1.96*std.error,
                                              xmax = estimate + 1.96*std.error)) +
  geom_pointrange() +
  facet_wrap( ~ term, scales = "free") +
  geom_vline(xintercept=0, col="black", linetype = "dashed", size=0.5) +
  labs(x = "Estimate", y = "Taxa") +
  theme_bw() + theme(
    strip.text = element_text(size = 15, face="bold"),
    legend.position="bottom",
    legend.title = element_text(size=10, face="bold"),
    legend.text = element_text(size=10),
    axis.title.x = element_text(size=15, face="bold"),
    axis.title.y = element_text(size=15, face="bold"))


png("Fig4_Final.png",width =1800, height = 900, res=200, bg = "white")
Coef_plot
dev.off()



