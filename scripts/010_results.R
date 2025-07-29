# Kyle's dissertation paper on land use and river nutrients
# Script edited by Christian John 18 July 2024 thru most recent update

# Kyle Neumann, Christian John, and Jordan Hollarsmith
# Chapter 2 of Kyle's Dissertation


#### Import packages ----
library(sf)
library(vegan)
library(lmerTest)
library(visreg)
library(tidyverse)
library(ggfortify)
library(ggpattern)
library(gginnards)
source("scripts/000_colorPalettes.R")


#### Define parameters ----
set.seed(07291954) ## Fellowship of the Ring published

#### Load data and some wrangling ----
## .. Watershed data ----
watershed_parameters <- read_sf("data/WaterSample_watersheds.shp") %>%
  st_drop_geometry() %>%
  rename("Watershed" = Watrshd,
         "Area" = Area_km,
         "StreamOrder" = StrmOrd,
         "StreamLength" = Strm_km) %>%
  select(Watershed,
         Area,
         StreamOrder,
         StreamLength) %>%
  mutate(Watershed = case_when(Watershed == "Haapiti" ~ "Ha'apiti",
                               Watershed == "Maatea" ~ "Ma'atea",
                               TRUE ~ Watershed))

## .. River data ----
## Dataset with Mg/L calculated in addition to uM
chp2rawdata <- read_csv("data/RiverData.csv",
                        na = c("","NA", "#VALUE!", "#DIV/0!")) %>% 
  # add id column for every sample event
  # make id column numeric
  mutate(id=rownames(.),
         id = as.numeric(id)) %>%
  # change month and monthday values in the singular missing date row to NA
  # Note that that date will fail to parse
  mutate(Sample_Month = ifelse(Sample_Month == -1, 
                               NA, 
                               Sample_Month),
         Sample_MonthDay = ifelse(Sample_MonthDay == -1,
                                  NA,
                                  Sample_MonthDay),
         Sample_Date = lubridate::ymd(paste(Sample_Year,
                                            Sample_Month,
                                            Sample_MonthDay)))

## Modify Paopao (Pao Pao), Ha'apiti (Haapiti), and Ma'atea (Maatea) syntax 
## to match other datasets
## Remove Afareaitu and Vaiare from every analysis
## TSS less than 0 are an artefact of clear water. Raise these values to 0.
river.data = chp2rawdata %>%
  dplyr::select("Sample_Date", "Sample_Year", "Season",
                "Site", "Shore", 
                "NH4_uM", "NH4_mgL", 
                "NO2_uM", "NO2_mgL", 
                "PO4_uM", "PO4_mgL",
                "NO3_uM", "NO3_mgL", 
                "DIN_uM", "DIN_mgL", 
                "N.P", 
                "TSS_mgL") %>%
  mutate(Site = gsub("Pao Pao", "Paopao", Site),
         Site = ifelse(Site == "Haapiti", "Ha'apiti", Site),
         Site = ifelse(Site == "Maatea", "Ma'atea", Site)) %>%
  filter(!( Site %in% c("Afareaitu", "Vaiare"))) %>%
  mutate(TSS_mgL = ifelse(TSS_mgL < 0, 0, TSS_mgL))

## .. Census data ----
river.pops = read_csv("data/census_by_watershed.csv") %>%
  dplyr::select(-IDDIST17) %>%
  group_by(Site) %>%
  summarize(Population = sum(Ensemble)) %>%
  mutate(Site = ifelse(Site == "Vaiane", "Vaianae", Site)) 

## .. Land use data ----
water.shed <- read.csv("data/landuse_wide2018.csv") %>%
  mutate(Watershed = case_when(Watershed == "Haapiti" ~ "Ha'apiti",
                               Watershed == "Maatea" ~ "Ma'atea",
                               TRUE ~ Watershed)) %>%
  ## Create a column with "cleared" as the sum of non-forested, non-water area
  ## And then calculate actual area of classes by watershed
  mutate(Cleared = BufferAg + Buildings + Dirt + Monoculture + Paved,
         BufferAg_area = BufferAg*WatershedArea,
         Buildings_area = Buildings*WatershedArea,
         Dirt_area = Dirt*WatershedArea,
         Forest_area = Forest*WatershedArea,
         Monoculture_area = Monoculture*WatershedArea,
         Paved_area = Paved*WatershedArea,
         Water_area = Water*WatershedArea,
         Cleared_area = Cleared*WatershedArea) %>%
  rename(Site = Watershed) 

## .. Precipitation data ----
## These were generated in 002_meteoFrance.R
precip = read_csv("data/rolling_precip.csv") %>%
  mutate(Sample_Date = Date) %>%
  rename(Precipitation = rollPrecip) %>%
  select(Sample_Date, Shore, Precipitation)

precipMonthly = read_csv("data/monthly_precip.csv")


#### Methods summary values ----
## .. Sampling effort ----
## In total, 228 samples were collected over the four sampling seasons (37 and 56 samples in the 2018 dry season and rainy season, respectively, and 53 and 82 samples in the 2019 dry season and rainy season, respectively).
## Summarize sampling effort
sampleSummary = chp2rawdata %>%
  filter(!(is.na(NH4_uM) & is.na(NO2_uM) &
             is.na(PO4_uM) & is.na(NO3_uM) &
             is.na(DIN_uM) & is.na(N.P) & 
             is.na(TSS_mgL)))
## total number of samples
sampleSummary %>%
  count()
## replicates per sampling season
sampleSummary %>%
  count(Sample_Year, Season)

## .. Detection limits ----
## x = micromolarity
## y = molar mass
## solving for z concentration in mg/L
## x uM =  x * 1e-6 * mol / L
## y g/mol * 1000 mg/g = y * 1e3 mg/mol
## x * 1e-6 mol/L * y * 1e3 mg/mol = z mg/L
## x * y * 1e-3 = z mg/L

detectionNH4nM = 50 # micromol/L
(detectionNH4nM*(18.04))/1000 # mg/L

detectionPO4nM = 3
(detectionPO4nM*(94.97))/1000

detectionNO3nM = 2
(detectionNO3nM*(62.01))/1000

detectionNO2nM = 2
(detectionNO2nM*(46.01))/1000

#### Data wrangling ----
## .. Combine datasets ----
mydata = river.data %>%
  mutate(Site = case_when(Site == "Haapiti" ~ "Ha'apiti",
                          Site == "Maatea" ~ "Ma'atea",
                          TRUE ~ Site)) %>%
  select(Sample_Date, Sample_Year, Season, Site, Shore, 
         NH4_uM, NH4_mgL, 
         NO2_uM, NO2_mgL, 
         NO3_uM, NO3_mgL, 
         PO4_uM, PO4_mgL, 
         DIN_uM, DIN_mgL, 
         N.P, 
         TSS_mgL) %>%
  left_join(water.shed) %>%
  left_join(river.pops) %>%
  left_join(precip)


#### Manuscript Tables 1 and S1 ----
## .. > Table 1 ----
## Summarizing watershed parameters
table1 = water.shed %>%
  select(Site, WatershedArea, Cleared) %>%
  left_join(river.pops) %>%
  rename("Watershed" = Site) %>%
  filter(Watershed != "Hauru") %>%
  left_join(watershed_parameters) %>%
  arrange(Watershed) %>%
  mutate(WatershedArea = round(WatershedArea/1000000, 1),
         Cleared = round(Cleared*100,1),
         Shore = c("W","W","E","E","N","N","N","N",
                   "N","N","N","N","N","N","E","W"),
         StreamLength = round(StreamLength, 1)) %>%
  select(Watershed, Shore, WatershedArea, StreamOrder, StreamLength,
         Population, Cleared) %>%
  rename("Watershed area (km2)" = WatershedArea,
         "Stream order" = StreamOrder,
         "Stream length (km)" = StreamLength,
         "Cleared (%)" = Cleared)
table1
sjPlot::tab_df(table1,file = "tables/table1.doc",
               col.header = names(table1))
write_csv(table1,file = "tables/table1.csv")

## .. > Supplementary Materials S1 ----
## Have to make an intermediate object in order to point toward its col names
tableS1 = river.data %>%
  select(Site, DIN_mgL, PO4_mgL, N.P, TSS_mgL) %>%
  pivot_longer(-Site, names_to = "variable", values_to = "value") %>%
  group_by(Site, variable) %>%
  filter(sum(!is.na(value)) > 1) %>%
  summarize(sd = round(sd(value, na.rm = T),2),
            minVal = round(min(value, na.rm=T),2),
            maxVal = round(max(value, na.rm=T),2),
            mean = round(mean(value, na.rm=T),2)) %>%
  ungroup() %>%
  mutate(range = paste0(minVal,"-",maxVal)) %>%
  select(Site, variable, mean, sd, range) %>%
  pivot_wider(names_from = variable, 
              values_from = c(mean, range, sd),
              names_sep = "_") %>%
  select(Site, 
         names(.)[grep("DIN_mgL",names(.))], 
         names(.)[grep("PO4",names(.))], 
         names(.)[grep("N.P",names(.))], 
         names(.)[grep("TSS",names(.))])
tableS1
sjPlot::tab_df(tableS1,file = "tables/S1.doc",
               col.header = names(tableS1))
write.csv(tableS1,file = "tables/S1.csv")


#### Manuscript text ----
## .. Results paragraph 1 ----
## .. > Figure 2 ----
## Precipitation in Moorea was seasonal but occurred during every month of the study period (Figure 2). 
## Note, Fig 2 is produced in 002_methoFrance.R.

## .... Precipitation ----
samplingDates = river.data %>%
  pull(Sample_Date) %>%
  unique()
samplingMonths = paste0(lubridate::year(samplingDates),
                        "-",
                        lubridate::month(samplingDates),"-01")
precipMonthlySummary = precipMonthly %>%
  mutate(sm = paste0(lubridate::year(Date),
                     "-",lubridate::month(Date),"-01")) %>%
  filter(sm %in% samplingMonths)

## Within our sampling dates, February of 2018 had the highest cumulative monthly rainfall with 866 mm of precipitation on the North shore, 659 mm on the East shore, and 772 mm on the West shore. 
precipMonthlySummary %>%
  group_by(Shore) %>%
  summarize(maxPrec = max(Precip_mm,na.rm=T),
            maxDate = Date[which(Precip_mm == max(Precip_mm, na.rm=T))])


## The rainy season in 2019 had lower monthly precipitation than in 2018 with a peak of 234 mm in February on the West shore, and 465 and 274 mm in March on the North shore and East shore, respectively.
precipMonthlySummary %>%
  filter(lubridate::year(Date) == 2019) %>%
  group_by(Shore) %>%
  summarize(maxPrec = max(Precip_mm,na.rm=T),
            maxDate = Date[which(Precip_mm == max(Precip_mm, na.rm=T))])

## August was the driest month of both years and precipitation totaled 103 mm on the North shore, 56 mm on the East shore, and 39 mm on the West shore in 2018; and 86 mm of on the North shore, 23 mm on the East shore, and 28 mm on the West shore in 2019. 
precipMonthlySummary %>%
  filter(lubridate::year(Date) == 2018) %>%
  group_by(Shore) %>%
  summarize(minPrec = min(Precip_mm,na.rm=T),
            minDate = Date[which(Precip_mm == min(Precip_mm, na.rm=T))])
precipMonthlySummary %>%
  filter(lubridate::year(Date) == 2019) %>%
  group_by(Shore) %>%
  summarize(minPrec = min(Precip_mm,na.rm=T),
            minDate = Date[which(Precip_mm == min(Precip_mm, na.rm=T))])


## .. Results paragraph 2 ----
## Water chemistry varied considerably by nutrient, season, and watershed across the sampling period (Figure 3, Supplementary materials S1). 
## .. > Figure 3 ----
## Prep data for Fig 3, removing Pihaena since there are no dry season obs and sorting site and season the way we want
fig3Data = mydata %>%
  mutate(Season = ifelse(Season == "d", "dry","rainy")) %>%
  filter(Site != "Pihaena") %>%
  mutate(Site = factor(Site, levels = rev(sort(unique(Site)))),
         Season = factor(Season, levels = c("rainy","dry"))) 

## ...... Generate river chemistry boxplots ----
fig3a = fig3Data %>%
  ggplot(aes(x=Site, y=DIN_mgL, fill=Site)) + 
  geom_hline(yintercept = 0.1,
             col = "grey60",
             lwd = 1,
             lty = 1) +
  geom_hline(yintercept = 0.18,
             col = "grey60",
             lwd = 1,
             lty = 2) +
  ggpattern::geom_boxplot_pattern(aes(pattern = Season),
                                  colour          = 'black',
                                  pattern_fill = "black",
                                  pattern_colour  = 'black',
                                  pattern_size = 0.01,
                                  pattern_spacing = 0.015,
                                  pattern_angle = 45,
                                  outlier.size = 0.35,
                                  size = 0.15,
                                  width = 1) +
  scale_fill_manual(values = simpsonCols) +
  scale_pattern_manual(name="Season", 
                       values = c("rainy"="stripe","dry"="none")) +
  ylab("DIN (mg/L)") +
  xlab("Watershed") +
  coord_flip() +
  CJsBasics::BasicTheme +
  theme(legend.position = "none",
        axis.title.y = element_blank())
fig3b = fig3Data %>%
  ggplot(aes(x=Site, y=NO3_mgL, fill=Site)) + 
  ggpattern::geom_boxplot_pattern(aes(pattern = Season),
                                  colour          = 'black',
                                  pattern_fill = "black",
                                  pattern_colour  = 'black',
                                  pattern_size = 0.01,
                                  pattern_spacing = 0.015,
                                  pattern_angle = 45,
                                  outlier.size = 0.35,
                                  size = 0.15,
                                  width = 1) +
  scale_fill_manual(values = simpsonCols) +
  scale_pattern_manual(name="Season", 
                       values = c("rainy"="stripe","dry"="none")) +
  ylab(expression(NO[3]^" -"*" (mg/L)")) +
  xlab("Watershed") +
  coord_flip() +
  CJsBasics::BasicTheme +
  theme(legend.position = "none",
        axis.title.y = element_blank())
fig3c = fig3Data %>%
  ggplot(aes(x=Site, y=NH4_mgL, fill=Site)) + 
  ggpattern::geom_boxplot_pattern(aes(pattern = Season),
                                  colour          = 'black',
                                  pattern_fill = "black",
                                  pattern_colour  = 'black',
                                  pattern_size = 0.01,
                                  pattern_spacing = 0.015,
                                  pattern_angle = 45,
                                  outlier.size = 0.35,
                                  size = 0.15,
                                  width = 1) +
  scale_fill_manual(values = simpsonCols) +
  scale_pattern_manual(name="Season", 
                       values = c("rainy"="stripe","dry"="none")) +
  ylab(expression(NH[4]^" +"*" (mg/L)")) +
  xlab("Watershed") +
  coord_flip() +
  CJsBasics::BasicTheme +
  theme(legend.position = "none",
        axis.title.y = element_blank())
fig3d = fig3Data %>%
  ggplot(aes(x=Site, y=NO2_mgL, fill=Site)) + 
  ggpattern::geom_boxplot_pattern(aes(pattern = Season),
                                  colour          = 'black',
                                  pattern_fill = "black",
                                  pattern_colour  = 'black',
                                  pattern_size = 0.01,
                                  pattern_spacing = 0.015,
                                  pattern_angle = 45,
                                  outlier.size = 0.35,
                                  size = 0.15,
                                  width = 1) +
  scale_fill_manual(values = simpsonCols) +
  scale_pattern_manual(name="Season", 
                       values = c("rainy"="stripe","dry"="none")) +
  ylab(expression(NO[2]^" -"*" (mg/L)")) +
  xlab("Watershed") +
  coord_flip() +
  CJsBasics::BasicTheme +
  theme(legend.position = "none",
        axis.title.y = element_blank())
fig3e = fig3Data %>%
  ggplot(aes(x=Site, y=PO4_mgL, fill=Site)) + 
  ggpattern::geom_boxplot_pattern(aes(pattern = Season),
                                  colour          = 'black',
                                  pattern_fill = "black",
                                  pattern_colour  = 'black',
                                  pattern_size = 0.01,
                                  pattern_spacing = 0.015,
                                  pattern_angle = 45,
                                  outlier.size = 0.35,
                                  size = 0.15,
                                  width = 1) +
  scale_fill_manual(values = simpsonCols) +
  scale_pattern_manual(name="Season", 
                       values = c("rainy"="stripe","dry"="none")) +
  ylab(expression(PO[4]^" 3-"*" (mg/L)")) +
  xlab("Watershed") +
  coord_flip() +
  CJsBasics::BasicTheme +
  theme(legend.position = "none",
        axis.title.y = element_blank())
fig3f = fig3Data %>%
  ggplot(aes(x=Site, y=N.P, fill=Site)) + 
  ggpattern::geom_boxplot_pattern(aes(pattern = Season),
                                  colour          = 'black',
                                  pattern_fill = "black",
                                  pattern_colour  = 'black',
                                  pattern_size = 0.01,
                                  pattern_spacing = 0.015,
                                  pattern_angle = 45,
                                  outlier.size = 0.35,
                                  size = 0.15,
                                  width = 1) +
  scale_fill_manual(values = simpsonCols) +
  scale_pattern_manual(name="Season", 
                       values = c("rainy"="stripe","dry"="none")) +
  ylab("N:P") +
  xlab("Watershed") +
  coord_flip() +
  CJsBasics::BasicTheme +
  theme(legend.position = "none",
        axis.title.y = element_blank())
fig3g = fig3Data %>%
  ## Raise less-than-1 mg/L values to 1 to facilitate log scale
  ## These values are so small this transformation isn't significant for results
  ## but fixes the x-intercept/log scale issue
  mutate(TSS_mgL = ifelse(TSS_mgL < 1, 1, TSS_mgL)) %>%
  ggplot(aes(x=Site, y=TSS_mgL, fill=Site)) + 
  geom_hline(yintercept = 5,
             col = "grey60",
             lwd = 1,
             lty = 1) +
  geom_hline(yintercept = 50,
             col = "grey60",
             lwd = 1,
             lty = 2) +
  ggpattern::geom_boxplot_pattern(aes(pattern = Season),
                                  colour          = 'black',
                                  pattern_fill = "black",
                                  pattern_colour  = 'black',
                                  pattern_size = 0.01,
                                  pattern_spacing = 0.015,
                                  pattern_angle = 45,
                                  outlier.size = 0.35,
                                  size = 0.15,
                                  width = 1) +
  scale_fill_manual(values = simpsonCols) +
  scale_pattern_manual(name="Season", 
                       values = c("rainy"="stripe","dry"="none")) +
  scale_y_continuous(trans = "log",
                     breaks = c(1,10,100,1000),
                     limits = c(1,1000)) +
  ylab("TSS (mg/L)") +
  xlab("Watershed") +
  coord_flip() +
  CJsBasics::BasicTheme +
  theme(legend.position = "none",
        axis.title.y = element_blank())
# fig3g
## ...... Compile Figure 3 ----
fig3Compiled = cowplot::plot_grid(
  cowplot::plot_grid(
    fig3a+theme(plot.margin = margin(0,0,0,0,"in")),
    fig3b+theme(plot.margin = margin(0,0,0,0,"in")),
    fig3c+theme(plot.margin = margin(0,0,0,0,"in")),
    fig3d+theme(plot.margin = margin(0,0,0,0,"in")),
    fig3e+theme(plot.margin = margin(0,0,0,0,"in")),
    fig3f+theme(plot.margin = margin(0,0,0,0,"in")),
    nrow = 3,
    ncol = 2,
    labels = c("a","b","c","d","e","f")),
  fig3g+theme(plot.margin = margin(0,0,0,0,"in")),
  labels = c("","g"),
  nrow = 2,
  rel_heights = c(3,1))
## ...... Save Figure 3 ----
ggsave("plots/figure3.jpg",fig3Compiled,
       width = 8, height = 10, units = "in", dpi = 300)


## .... DIN ----
# DIN concentrations ranged from a minimum of 0.008 mg/L measured in a sample collected from Ma’atea in 2019, to 0.528 mg/L measured in a sample collected in Paopao in 2018.
mydata %>%
  filter(Site != "Pihaena") %>%
  select(Sample_Date, DIN_mgL, Site) %>%
  filter(DIN_mgL %in% c(min(DIN_mgL, na.rm = T),
                        max(DIN_mgL, na.rm = T)))
# DIN differed significantly across watersheds (p < 0.001) and seasons (p=0.04; Figure 3a, Supplementary materials S2). 
dinsummary = mydata %>%
  filter(Site != "Pihaena") %>%
  mutate(logN_mgL = log(DIN_mgL + 1))
shedDIN = lm(logN_mgL ~ Season*Site, data = dinsummary)
anova(shedDIN)
# We found no significant differences between seasons within each watershed in the post-hoc Tukey test, though the significant interaction between watershed and season in the overall model indicates that rainy season DIN concentrations were higher than dry season in some watersheds (e.g., Papetoai and Maharepa) while the opposite pattern was observed in other watersheds (e.g., Ha’apiti and Haumi). 
TukeyHSD(aov(shedDIN), "Season:Site")$`Season:Site` %>%
  as.data.frame() %>%
  mutate(comp = rownames(.)) %>%
  separate(comp, into = c("s1","s2"),sep="-") %>%
  separate(s1, into = c("sea1","sit1"),sep=":") %>%
  separate(s2, into = c("sea2","sit2"),sep=":") %>%
  filter(sit1 == sit2) %>%
  filter(`p adj` < 0.05)
visreg::visreg(shedDIN, "Season","Site")
# The watershed with the highest mean±s.e. DIN concentration was Paopao (0.256±0.018 mg/L), and the lowest mean concentration was found in Vaianae (0.057±0.004 mg/L).
dinsummary %>% 
  group_by(Site) %>%
  summarize(sedin = CJsBasics::se(DIN_mgL, na.rm = T),
            meandin = mean(DIN_mgL, na.rm = T)) %>%
  arrange(meandin)

dinsummary  %>% 
  summarize(meandinMGL = mean(DIN_mgL, na.rm = T),
            meandinuM = mean(DIN_uM, na.rm = T))


## .. Results paragraph 3 ----
## .... Nitrate ----
# The N species comprising the majority of DIN was NO3-, followed by NH4+ and finally NO2-. 
mydata %>% 
  ## Check this using molarity rather than full weight since we care 
  ## about the N part and not the O, H parts
  mutate(compNO3 = NO3_uM/DIN_uM,
         compNO2 = NO2_uM/DIN_uM,
         compNH4 = NH4_uM/DIN_uM) %>%
  group_by(Site) %>%
  summarize(NO3 = mean(compNO3, na.rm = T),
            NO2 = mean(compNO2, na.rm = T),
            NH4 = mean(compNH4, na.rm = T)) %>%
  pivot_longer(-Site) %>%
  ggplot() + 
  geom_density(aes(x = value, 
                   group = name, 
                   fill = name), 
               alpha = 0.5) +
  scale_fill_discrete("Species") +
  xlab("Proportion contribution") +
  ylab("Density") +
  CJsBasics::BasicTheme +
  theme(legend.position = c(0.85,0.85))

## NO3- concentrations differed significantly by watershed (p < 0.001) and season (p = 0.02), and there was a significant interaction between watershed and season (p < 0.001); like DIN, the post-hoc Tukey test revealed no significant differences in NO3- concentrations between seasons within watersheds (Figure 3b).
no3summary = mydata %>%
  filter(Site != "Pihaena") %>%
  mutate(logNO3mgL = log(NO3_mgL + 1))
shedNO3 = lm(logNO3mgL ~ Season*Site, data = no3summary)
anova(shedNO3)
TukeyHSD(aov(shedNO3), "Season:Site")$`Season:Site` %>%
  as.data.frame() %>%
  mutate(comp = rownames(.)) %>%
  separate(comp, into = c("s1","s2"),sep="-") %>%
  separate(s1, into = c("sea1","sit1"),sep=":") %>%
  separate(s2, into = c("sea2","sit2"),sep=":") %>%
  filter(sit1 == sit2) %>%
  filter(`p adj` < 0.05)

## .... Ammonium ----
# NH4+ concentrations differed significantly by watershed (p < 0.001), while the interaction with season was not significant. 
nh4summary = mydata %>%
  filter(Site != "Pihaena") %>%
  mutate(logNH4mgL = log(NH4_mgL + 1))
shedNH4 = lm(logNH4mgL ~ Season*Site, data = nh4summary)
anova(shedNH4)
TukeyHSD(aov(shedNH4),"Site")$Site %>% as.data.frame() %>% 
  filter(`p adj` < 0.05) %>%
  arrange(diff)

# Teavaro and Haumi had higher NH4+ levels than the other watersheds (Figure 3c; the mean±s.e. NH4+ concentration in these watersheds was 0.030±0.006 mg/L and 0.023±0.003 mg/L, respectively, whereas across watersheds the mean NH4+ concentration was 0.012±0.002 mg/L). 
nh4summary %>% 
  group_by(Site) %>%
  summarize(seNH4=CJsBasics::se(NH4_mgL, na.rm = T),
            meanNH4 = mean(NH4_mgL, na.rm = T)) %>%
  arrange(meanNH4)
nh4summary %>% 
  group_by(Site) %>%
  summarize(meanNH4 = mean(NH4_mgL, na.rm = T)) %>%
  ungroup() %>%
  summarize(seNH4=CJsBasics::se(meanNH4, na.rm = T),
            meanNH4 = mean(meanNH4, na.rm = T)) 

## .... Nitrite ----
# NO2- concentrations differed significantly by watershed (p < 0.001) but not by season although there was a significant interaction of watershed and season (p = 0.001). 
no2summary = mydata %>%
  filter(Site != "Pihaena") %>%
  mutate(logNO2mgL = log(NO2_mgL + 1))
shedNO2 = lm(logNO2mgL ~ Season*Site, data = no2summary)
anova(shedNO2)
# The rainy season NO2- concentrations were significantly lower than the dry season concentrations in Haumi (difference±s.e. = 0.030±0.011 mg/L, p < 0.001). 
TukeyHSD(aov(shedNO2), "Season:Site")$`Season:Site` %>%
  as.data.frame() %>%
  mutate(comp = rownames(.)) %>%
  separate(comp, into = c("s1","s2"),sep="-") %>%
  separate(s1, into = c("sea1","sit1"),sep=":") %>%
  separate(s2, into = c("sea2","sit2"),sep=":") %>%
  filter(sit1 == sit2) %>%
  filter(`p adj` < 0.05) %>%
  mutate(se = (diff-lwr)/1.96)
TukeyHSD(aov(shedNO2),"Site")$Site %>% as.data.frame() %>% 
  filter(`p adj` < 0.05) %>%
  arrange(diff)
# Like NH4+, Teavaro and Haumi had higher NO2- values relative to the other watersheds (Figure 3d; NO2- concentration in these watersheds was 0.035±0.007 mg/L and 0.033±0.008 mg/L, respectively, whereas across watersheds the mean NH4+ concentration was 0.009±0.003 mg/L). 
no2summary %>% 
  # filter(Site %in% c("Teavaro","Haumi")) %>%
  group_by(Site) %>%
  summarize(seNO2=CJsBasics::se(NO2_mgL, na.rm = T),
            meanNO2 = mean(NO2_mgL, na.rm = T)) %>%
  arrange(meanNO2)
no2summary %>% 
  group_by(Site) %>%
  summarize(meanNO2 = mean(NO2_mgL, na.rm = T)) %>%
  ungroup() %>%
  summarize(seNO2=CJsBasics::se(meanNO2, na.rm = T),
            meanNO2 = mean(meanNO2, na.rm = T)) 

## .. Results paragraph 4 ----
## .... Phosphate ----
# PO43- concentrations ranged from a minimum of 0.018 mg/L from a sample collected in Teavaro in the dry season to a maximum of 0.72 mg/L in a sample collected from Paopao 3 in the rainy season. 
po4summary = mydata %>%
  filter(Site != "Pihaena") %>%
  mutate(logPO4mgL = log(PO4_mgL+1))
po4summary %>%
  group_by(Site) %>%
  summarize(minPO4 = min(PO4_mgL,na.rm=T),
            minDate = Sample_Date[which(PO4_mgL == min(minPO4, na.rm=T))]) %>%
  arrange(minPO4)
po4summary %>%
  group_by(Site) %>%
  summarize(maxPO4 = max(PO4_mgL,na.rm=T),
            maxDate = Sample_Date[which(PO4_mgL == max(maxPO4, na.rm=T))]) %>%
  arrange(maxPO4)
# PO43- differed significantly among watersheds (p < 0.001) and season (p = 0.03) and had a nearly significant interaction of watershed and season (p = 0.053) (Figure 3e). 
shedPO4 = lm(logPO4mgL ~ Season*Site, data = po4summary)
anova(shedPO4)
# The highest mean PO43- concentration from rivers that flowed year-round was observed in Paopao 3 (0.257±0.053 mg/L) while the lowest was in Paopao 1 (0.093±0.024 mg/L). 
mydata %>%
  filter(Site != "Pihaena") %>%
  group_by(Site) %>%
  summarize(sePO4=CJsBasics::se(PO4_mgL, na.rm = T),
            meanPO4 = mean(PO4_mgL, na.rm=T)) %>%
  arrange(meanPO4)

#  PO43- concentrations were significantly associated with TSS concentrations (β = 0.04, p < 0.001, and r2 = 0.33). 
summary(lm(PO4_mgL ~ log(TSS_mgL+1), data = po4summary))

## .... N:P Ratio ----
## The N:P ratio ranged between a minimum of 0.29 in Opunohu 1, and a maximum of 54.7 in Paopao 1. 
mydata %>%
  filter(Site != "Pihaena") %>%
  group_by(Site) %>%
  summarize(minNP = min(N.P,na.rm=T)) %>%
  arrange(minNP)
mydata %>%
  filter(Site != "Pihaena") %>%
  group_by(Site) %>%
  summarize(maxNP = max(N.P,na.rm=T)) %>%
  arrange(maxNP)

## The mean N:P ratio was highest in Paopao 1 (24.2±5.97) and lowest in Maharepa (2.29±0.64).
mydata %>%
  filter(Site != "Pihaena") %>%
  group_by(Site) %>%
  summarize(seNP=CJsBasics::se(N.P, na.rm = T),
            meanNP = mean(N.P,na.rm=T)) %>%
  arrange(meanNP)

# N:P differed significantly across watersheds (p < 0.001) and had a significant interaction with watershed and season (p < 0.001), with Paopao 1 and Opunohu 2 having significantly higher N:P molar ratio in dry season than the rainy season (Figure 3).
npsummary = mydata %>%
  filter(Site != "Pihaena") 
shedNP = lm(N.P ~ Season*Site, data = npsummary)
anova(shedNP)
TukeyHSD(aov(shedNP), "Season:Site")$`Season:Site` %>%
  as.data.frame() %>%
  mutate(comp = rownames(.)) %>%
  separate(comp, into = c("s1","s2"),sep="-") %>%
  separate(s1, into = c("sea1","sit1"),sep=":") %>%
  separate(s2, into = c("sea2","sit2"),sep=":") %>%
  filter(sit1 == sit2) %>%
  filter(`p adj` < 0.05)


## .. Results paragraph 5 ----
## .... Total suspended solids ----
# TSS concentrations ranged from below the detection limit in some samples from Atiha, Ha’apiti, and Maharepa to a maximum of 902 mg/L found in a sample collected in Paopao 1 during the rainy season. 
tsssummary = mydata %>%
  filter(!(Site %in% c("Pihaena",
                       "Paopao",
                       "Paopao 3")))  %>%
  mutate(logTSS_mgL = log(TSS_mgL + 1))
tsssummary %>%
  group_by(Site) %>%
  summarize(minTSS = min(TSS_mgL,na.rm=T),
            minDate = Sample_Date[which(TSS_mgL == min(minTSS, na.rm=T))]) %>%
  arrange(minTSS)
tsssummary %>%
  group_by(Site) %>%
  summarize(maxTSS = max(TSS_mgL,na.rm=T),
            maxDate = Sample_Date[which(TSS_mgL == max(maxTSS, na.rm=T))]) %>%
  arrange(maxTSS)
# Mean TSS was higher during the rainy season than the dry season in 10 of the 12 watersheds for which year-round TSS data were collected (i.e. all sites except Pihaena, Paopao, and Paopao 3). Differences in TSS among seasons were significant (p = 0.002) but the interaction between watershed and season was not significant (p > 0.05) (Figure 3). 
shedTSS = lm(logTSS_mgL ~ Season*Site, data = tsssummary)
anova(shedTSS)
visreg::visreg(shedTSS, "Season", "Site")

# Although differences in TSS among watersheds were not significant after accounting for season, variability in TSS concentration across watersheds was impressive: the lowest mean TSS concentration was recorded in Ma’atea at 2.5±0.44 mg/L while the highest was recorded in Paopao 1 at 137±128 mg/L. The next highest mean TSS concentration was in Paopao 2 at 117±111 mg/L, followed by Papetoai with 46±36 mg/L, less than half the concentration of the two largest Paopao sub-watersheds.
tsssummary %>%
  group_by(Site) %>%
  filter(sum(!is.na(TSS_mgL)) > 1) %>%
  summarize(seTSS=CJsBasics::se(TSS_mgL, na.rm = T),
            meanTSS = mean(TSS_mgL,na.rm=T)) %>%
  arrange(meanTSS)


## .. > Table S2 ----
## .... ANOVA table outputs ----
myANOVA_output = function(x, response){
  myANOVA = anova(x)
  as.data.frame(myANOVA) %>%
    mutate("Response" = response,
           "Predictor" = rownames(.),
           `Sum Sq` = round(`Sum Sq`, 2),
           `Mean Sq` = round(`Mean Sq`, 2),
           `F value` = round(`F value`, 2),
           `Pr(>F)` = round(`Pr(>F)`, 3)) %>%
    select(Response, Predictor, Df, `Sum Sq`, `Mean Sq`, `F value`, `Pr(>F)`)
}

myANOVA_output(shedDIN, "DIN") %>%
  rbind(myANOVA_output(shedNO2, "NO2")) %>%
  rbind(myANOVA_output(shedNO3, "NO3")) %>%
  rbind(myANOVA_output(shedNH4, "NH4")) %>%
  rbind(myANOVA_output(shedPO4, "PO4")) %>%
  rbind(myANOVA_output(shedTSS, "TSS")) %>%
  mutate(pStar = CJsBasics::pStars(`Pr(>F)`)) %>%
  write_csv(file = "tables/S2.csv")




## .. Results paragraph 6 ----
## .... Land use ----
# Our analysis of land use focused on the percent cleared land – and specifically percent of land being used for monoculture –  using the consensus land cover classification from a series of Worldview-3 images collected in 2018 (Figure 1). Island-wide model accuracy was 98% based on a spatially random sample of validation points. Using a hierarchical sampling scheme to balance across training classes, classification accuracy was 82%. Most of the classification discrepancies came from errors of omission of exposed soil which was often co-occurring with intensive monoculture and captured by that class, and errors of omission of small buildings which were lost as noise to their surrounding vegetation buffer. The watershed with the highest proportion of cleared land in 2018 was Paopao 2 (23%), followed by Opunohu 1 (10.4%) while the lowest was Ha’apiti (0.8%) (Table 1).  
table1 %>%
  arrange(`Cleared (%)`) %>%
  select(Watershed, `Cleared (%)`)

## The human population on Moorea was 17,463 residents in 2017, with the largest population in the Paopao basin (1281 residents) and the smallest in Opunohu 2 (10 residents) (Table 1).
table1 %>%
  arrange(Population) %>%
  select(Watershed, Population)

## Population was not significantly related to watershed size or percentage of the watershed that was classified as cleared land (p > 0.05 but β > 0 in all cases), but these null associations were driven by the large area and broad-scale land clearing in the Opunohu basin (which has a very low population density): among all watersheds except those in the Opunohu basin, population was positively associated with watershed size and percentage of watershed that was cleared (R2 = 0.60 and 0.32, and p = 0.002 and 0.04, respectively).
shedPops = water.shed %>%
  left_join(river.pops) %>%
  mutate(WatershedAreakm2 = WatershedArea/1000000)
shedPops %>% ggplot(aes(x = WatershedAreakm2, y = Population)) + 
  geom_point(aes(col = Site)) +
  stat_smooth(method = "lm") + 
  scale_color_manual(values = CJsBasics::KellyCols[4:20])
summary(lm(Population ~ WatershedAreakm2, data = shedPops))
summary(lm(Population ~ Cleared, data = shedPops))

## Removing Opunohu basin reveals positive associations between population and land use in other watersheds
shedPops2 = shedPops %>%
  filter(!(Site %in% c("Opunohu","Opunohu 1","Opunohu 2")))
shedPops2 %>% ggplot(aes(x = WatershedAreakm2, y = Population)) +
  geom_point(aes(col = Site)) + 
  stat_smooth(method = "lm") + 
  scale_color_manual(values = CJsBasics::KellyCols[4:20])
summary(lm(Population ~ WatershedAreakm2, data = shedPops2))
summary(lm(Population ~ Cleared, data = shedPops2))


## .. Results paragraph 7 ----
## .... PCA data prep ----
mydata_pca_prep <- mydata %>%
  select(Site, Season, 
         DIN_mgL, PO4_mgL, N.P, TSS_mgL, 
         Precipitation,
         Cleared, Population, WatershedArea) %>%
  filter(complete.cases(.)) %>%
  mutate(Season = ifelse(Season == "d", "dry","rainy")) %>%
  mutate(DIN_mgL = log(DIN_mgL + 1),
         PO4_mgL = log(PO4_mgL + 1),
         TSS_mgL = log(TSS_mgL + 1),
         Cleared = log(Cleared + 1),
         WatershedArea = log(WatershedArea + 1),
         Population = log(Population + 1),
         Precipitation = log(Precipitation + 1))

## Separate into river chemistry data (chem) and environmental (meta) for the PCA
pca_meta <- mydata_pca_prep %>%
  select(Site, Season, Precipitation, 
         Cleared, Population, WatershedArea)
pca_chem <- mydata_pca_prep %>% 
  select(-c(Site, Season, Precipitation, 
            Cleared, Population, WatershedArea))

## .... Run PCA ----
## Calculate principal components using water chemistry
pc <- prcomp(pca_chem, center = T, scale = T)


## .. > Figure 4 ----
## Sample points across the watersheds and seasons were distributed across the ordination space in the PCA (Figure 4). 
## Compile loadings from PCA
chem.pc <- as.data.frame(
  pc$rotation[1:4,1:2]) %>%
  mutate(species = rownames(.)) %>%
  mutate(specieslab = case_when(
    species == "N.P" ~ "N:P ratio",
    species == "DIN_mgL" ~ "DIN (mg/L)",
    species == "TSS_mgL" ~ "TSS (mg/L)",
    species == "PO4_mgL" ~ "PO<sub>4</sub><sup>3-</sup> (mg/L)"),
    PCarrow1 = 0.5*PC1,
    PCarrow2 = 0.5*PC2,
    PClab1 = PCarrow1 + 0.025*sign(PC1),
    PClab2 = PCarrow2 + 0.025*sign(PC2))


## ...... Compile Figure 4 ----
figure4 = autoplot(pc, 
                   data = mydata_pca_prep, 
                   colour = "Site", 
                   shape = "Season",
                   size=3,
                   # alpha = 0.75,
                   loadings = TRUE, 
                   loadings.colour = 'transparent',
                   loadings.label = TRUE, 
                   loadings.label.color = "transparent") +
  scale_color_manual(values = simpsonCols) +
  geom_segment(data=chem.pc,
               aes(x=0,xend=PCarrow1,y=0,yend=PCarrow2),
               arrow = arrow(length = unit(0.3, "cm")),
               lwd = 1,
               col = "grey70") +
  ggtext::geom_richtext(data=chem.pc,
                        aes(x=PClab1,y=PClab2,
                            label=specieslab),
                        size=4,
                        nudge_x = c(0.05,0,0,0),
                        col = "grey40",
                        fill = alpha(c("white"),0.9)) +
  xlim(-0.4,0.4) +
  coord_equal() +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.spacing.y = unit(0.0005,"in"),
        legend.key.height = unit(0.0005,"in"))
figure4_arrowsUnder = gginnards::move_layers(
  figure4, idx = c(4,5), position = "bottom")
## ...... Save Figure 4 ----
ggsave("plots/figure4.jpg",
       figure4_arrowsUnder, 
       width = 9, height = 4.5, units = "in", dpi = 600)

## .... PCA interpretation ----
## The following interpretation comes from Figure 4 and Table 3
## The first PC axis explained 46.6% of the variability in river chemistry, and was negatively associated with PO43- and TSS (PC1 loadings of -0.53 and -0.26, respectively), and positively correlated with N:P and DIN (PC1 loadings of 0.64 and 0.49, respectively). The second PC axis explained 37.2% of the variability and was negatively correlated with all of the watershed chemistry variables (PC2 loadings of -0.45, -0.68, -0.52, and -0.25, for PO43-, TSS, DIN, and N:P, respectively) (Table 2).
pc

## .. > Table 2 ----
## Extract rotations to a data.frame and save the table
pcdf = as.data.frame(pc$rotation) %>%
  mutate(Parameter = rownames(.)) %>%
  select(Parameter, PC1, PC2, PC3, PC4) %>%
  mutate(PC1 = round(PC1, 2),
         PC2 = round(PC2, 2),
         PC3 = round(PC3, 2),
         PC4 = round(PC4, 2))
write_csv(pcdf, file = "tables/table2.csv")


## The third and fourth PC axes explained 9.5% and 6.6% of the variability, respectively.
(pc$sdev^2)/sum(pc$sdev^2)

## Watershed identity and season were significantly related to the first two axes of the ordination space (r2 = 0.50 and 0.07, respectively; p = 0.001 in both cases). 
m2 <- envfit(pc ~ 
               Site +
               Season,
             data = mydata_pca_prep, perm = 999)
m2

## .. Results paragraph 8 ----
## .... Connecting watershed characteristics with river chemistry ----
## Prep daily data with log-transformed variables
sheddaily_summaries = mydata %>%
  mutate(logDIN = log(DIN_mgL + 1),
         logPO4 = log(PO4_mgL + 1),
         logTSS = log(TSS_mgL + 1),
         logCleared = log(Cleared + 1),
         logPopulation = log(Population + 1),
         logPrecipitation = log(Precipitation  + 1))

## Linear mixed effects models comparing river chemistry with environmental variables
lmerDIN = lmer(
  logDIN ~ 
    logCleared*Season +
    logPopulation +
    logPrecipitation +
    (1|Site),
  data = sheddaily_summaries) 
lmerTSS = lmer(
  logTSS ~ 
    logCleared*Season +
    logPopulation +
    logPrecipitation +
    (1|Site),
  data = sheddaily_summaries)
lmerPO4 = lmer(
  logPO4 ~ 
    logCleared*Season +
    logPopulation +
    logPrecipitation +
    (1|Site),
  data = sheddaily_summaries)

# Linear mixed effects models including population and precipitation, and cleared land percentage, season, and their interaction as fixed effects and watershed identity as a random effect revealed differential impacts of precipitation and cleared land on different components of surface water chemistry (Figure 5, Supplementary materials S3). DIN was significantly positively associated with recent precipitation (p = 0.001), and a significant positive interaction between land clearing and season revealed that DIN was significantly higher under intense land clearing during the rainy season (p = 0.002). Similarly, PO43- was significantly higher following recent precipitation (p < 0.001), and higher under increased land clearing during the rainy season (p = 0.046). TSS was not significantly associated with land clearing, but was significantly higher following recent precipitation (p < 0.001). 
lmerDIN %>%
  summary()
lmerPO4 %>%
  summary()
lmerTSS %>%
  summary()

MuMIn::r.squaredGLMM(lmerDIN)
MuMIn::r.squaredGLMM(lmerPO4)
MuMIn::r.squaredGLMM(lmerTSS)

## .. > Supplementary table S3 ----
summary(lmerDIN)
summary(lmerPO4)
summary(lmerTSS)

sjPlot::tab_model(lmerDIN,        
                  file = "tables/supplementary_DIN.doc")
sjPlot::tab_model(lmerTSS,        
                  file = "tables/supplementary_TSS.doc")
sjPlot::tab_model(lmerPO4,        
                  file = "tables/supplementary_PO4.doc")


dinPredsLogscale = predict(lmerDIN,
        newdata = data.frame(cleared = c(0,20),
                             population = c(100,100),
                             Season = c("r","r"),
                             precip = c(50,50),
                             Site = c("Pihaena",
                                      "Pihaena")) %>%
          mutate(logCleared = log(cleared + 1),
                 logPopulation = log(population + 1),
                 logPrecipitation = log(precip + 1))) 
exp(dinPredsLogscale - 1)[2]/exp(dinPredsLogscale - 1)[1]

tssPredsLogscale = predict(lmerTSS,
                   newdata = data.frame(cleared = c(20,20),
                                        population = c(100,100),
                                        Season = c("d","d"),
                                        precip = c(10,50),
                                        Site = c("Pihaena",
                                                 "Pihaena")) %>%
                     mutate(logCleared = log(cleared + 1),
                            logPopulation = log(population + 1),
                            logPrecipitation = log(precip + 1))) 
exp(tssPredsLogscale - 1)[2]/exp(tssPredsLogscale - 1)[1]

## Check for multicollinearity
## vif > 5 is cause for concern
car::vif(lmerDIN)
car::vif(lmerPO4)
car::vif(lmerTSS)
## all good.

## .. > Figure 5 ----
## A function to extract visreg outputs
myvisregFxCleared = function(mylmer){
  visregDat = visreg(mylmer, "logCleared", by="Season",
                     type = "contrast", plot = F) 
  fitData = visregDat$fit %>%
    mutate(cleared = (exp(logCleared) - 1)*100)
  resData = visregDat$res[c("logCleared", "Season", "visregRes")] %>%
    mutate(cleared = (exp(logCleared) - 1)*100)
  return(list(fitData, resData))
}

myvisregFxPrecip = function(mylmer){
  visregDat = visreg(mylmer, "logPrecipitation", by="Season",
                     type = "contrast", plot = F) 
  fitData = visregDat$fit %>%
    mutate(precip = (exp(logPrecipitation) - 1))
  resData = visregDat$res[c("logPrecipitation", "visregRes")] %>%
    mutate(precip = (exp(logPrecipitation) - 1))
  return(list(fitData, resData))
}

## Pull out visreg info for the two lmer objects
dinFitCleared = myvisregFxCleared(lmerDIN)[[1]] %>%
  as_tibble() %>%
  mutate(response = "DIN") %>%
  rename("Concentration" = logDIN)
dinResCleared =  myvisregFxCleared(lmerDIN)[[2]] %>%
  as_tibble() %>%
  mutate(response = "DIN")
tssFitCleared = myvisregFxCleared(lmerTSS)[[1]] %>%
  as_tibble() %>%
  mutate(response = "TSS") %>%
  rename("Concentration" = logTSS)
tssResCleared = myvisregFxCleared(lmerTSS)[[2]] %>%
  as_tibble() %>%
  mutate(response = "TSS")
po4FitCleared = myvisregFxCleared(lmerPO4)[[1]] %>%
  as_tibble() %>%
  mutate(response = "PO4") %>%
  rename("Concentration" = logPO4)
po4ResCleared = myvisregFxCleared(lmerPO4)[[2]] %>%
  as_tibble() %>%
  mutate(response = "PO4")


dinFitPrecip = myvisregFxPrecip(lmerDIN)[[1]] %>%
  as_tibble() %>%
  mutate(response = "DIN") %>%
  rename("Concentration" = logDIN)
dinResPrecip = myvisregFxPrecip(lmerDIN)[[2]] %>%
  as_tibble() %>%
  mutate(response = "DIN")
tssFitPrecip = myvisregFxPrecip(lmerTSS)[[1]] %>%
  as_tibble() %>%
  mutate(response = "TSS") %>%
  rename("Concentration" = logTSS)
tssResPrecip = myvisregFxPrecip(lmerTSS)[[2]] %>%
  as_tibble() %>%
  mutate(response = "TSS")
po4FitPrecip = myvisregFxPrecip(lmerPO4)[[1]] %>%
  as_tibble() %>%
  mutate(response = "PO4") %>%
  rename("Concentration" = logPO4)
po4ResPrecip = myvisregFxPrecip(lmerPO4)[[2]] %>%
  as_tibble() %>%
  mutate(response = "PO4")


## Now join into a single data.frame()
visreg_fitCleared = do.call("rbind", 
                            list(dinFitCleared, 
                                 tssFitCleared, 
                                 po4FitCleared)) %>%
  mutate(sea = as.character(Season),
         sea = ifelse(sea == "r", "rainy", "dry"))
visreg_resCleared = do.call("rbind", 
                            list(dinResCleared,
                                 tssResCleared,
                                 po4ResCleared)) %>%
  mutate(sea = as.character(Season),
         sea = ifelse(sea == "r", "rainy", "dry"))

visreg_fitPrecip = do.call("rbind", 
                           list(dinFitPrecip, 
                                tssFitPrecip, 
                                po4FitPrecip)) 
visreg_resPrecip = do.call("rbind", 
                           list(dinResPrecip,
                                tssResPrecip, 
                                po4ResPrecip)) 

## Create a label key to format fig 5 facet labs
labelkey = c(
  "DIN" = "DIN *' (mg/L)'",
  "TSS" = "TSS *' (mg/L)'",
  "PO4" = "PO[4]^ 3^-{} *' (mg/L)'"
)


## Plot the contrasts for DIN and TSS along a cleared land gradient
fig5a = ggplot() +
  geom_ribbon(data = visreg_fitPrecip,
              aes(x = precip, 
                  ymin = visregLwr, ymax = visregUpr,
                  group = paste0(response)),
              col = "#50505040",
              fill = "#99999940",
              alpha = 0.25) +
  geom_point(data = visreg_resPrecip,
             aes(x = precip,
                 y = visregRes,
                 # col = sea,
                 group = paste0(response)),
             col = "grey40",
             shape = 1,
             size = 2,
             stroke = 0.75) +
  geom_line(data = visreg_fitPrecip,
            aes(x = precip, 
                y = visregFit,
                # col = sea,
                group = paste0(response)),
            lwd = 1.5) +
  facet_wrap(~response,
             scales = "free_y",
             labeller = as_labeller(labelkey, label_parsed),
             strip.position = "right",
             nrow = 3) +
  ylab("Contrast") +
  xlab("Precipitation (mm)") +
  scale_color_manual("Season",
                     values = c("#d4b895", "#2297e6")) +
  CJsBasics::BasicTheme +
  theme(strip.text = element_text(size = 12))
fig5b = ggplot() +
  geom_ribbon(data = visreg_fitCleared,
              aes(x = cleared, 
                  ymin = visregLwr, ymax = visregUpr,
                  group = paste0(sea,response)),
              col = "#50505040",
              fill = "#99999940",
              alpha = 0.25) +
  geom_point(data = visreg_resCleared,
             aes(x = cleared,
                 y = visregRes,
                 col = sea,
                 group = paste0(response)),
             col = "grey40",
             shape = 1,
             size = 2,
             stroke = 0.75) +
  geom_line(data = visreg_fitCleared,
            aes(x = cleared, 
                y = visregFit,
                col = sea,
                group = paste0(sea,response)),
            lwd = 1.5) +
  facet_wrap(~response,
             scales = "free_y",
             labeller = as_labeller(labelkey, label_parsed),
             strip.position = "right",
             nrow = 3) +
  ylab("Contrast") +
  xlab("Watershed cleared (%)") +
  scale_color_manual("Season",
                     values = c("#d4b895", "#2297e6")) +
  CJsBasics::BasicTheme +
  theme(strip.text = element_text(size = 12),
        legend.position = "none")


fig5 = cowplot::plot_grid(fig5a,
                          fig5b,
                          nrow = 1,
                          labels = c("a","b"))
ggsave(filename = "plots/figure5.jpg",
       fig5,
       height = 6, width = 6, units = "in", dpi = 300)


## .. Results paragraph 9 ----
## Threshold exceedance rates
## .. > Table 3 ----
exceedanceData = mydata %>%
  select(Site,
         TSS_mgL,
         DIN_mgL)

# Water samples from all watersheds, except for Ma’atea, exceeded the lower TSS water quality threshold of 5 mg/L at least once. On average, watersheds on Moorea exceeded this threshold in 39% of sampling events (Table 3). The maximum exceedance rate was found in Pihaena at 78% of the time. Eight of the watersheds included in this study exceeded the upper TSS threshold of 50 mg/L at least once. Across Moorea, this threshold was exceeded on average 8% of the time. Six watersheds exceeded Hawaiian water quality standards with values of 50 mg/L at least 10% of the time (Haumi: 11%, Opunohu 1: 10%, Paopao 1 and 2: 14%, Papetoai: 13%, and Pihaena: 22%). All watersheds except Atiha and Vaianae exceeded both the lower and the upper DIN thresholds (Table 3). On average, stream samples across Moorea exceeded 0.1 mg/L DIN 46% of the time and 0.18 mg/L 20% of the time. Opunohu 2, Paopao 1, and Paopao 2 exceeded the lower standard in every Sample_ 
exceedances = exceedanceData %>%
  mutate(Site = case_when(Site == "Haapiti" ~ "Ha'apiti",
                          Site == "Maatea" ~ "Ma'atea",
                          Site == "Pao Pao" ~ "Paopao",
                          Site == "Pao Pao 1" ~ "Paopao 1",
                          Site == "Pao Pao 2" ~ "Paopao 2",
                          Site == "Pao Pao 3" ~ "Paopao 3",
                          TRUE ~ Site)) %>%
  group_by(Site) %>%
  filter(sum(!is.na(TSS_mgL))>1) %>%
  summarize(
    exc_tss_05 = 100*sum(
      TSS_mgL[!is.na(TSS_mgL)] > 5)/
      sum(!is.na(TSS_mgL)),
    exc_tss_50 = 100*sum(
      TSS_mgL[!is.na(TSS_mgL)] > 50)/
      sum(!is.na(TSS_mgL)),
    exc_din_1 = 100*sum(
      DIN_mgL[!is.na(DIN_mgL)] > 0.1)/
      sum(!is.na(DIN_mgL)),
    exc_din_18 = 100*sum(
      DIN_mgL[!is.na(DIN_mgL)] > 0.18)/
      sum(!is.na(DIN_mgL))) %>%
  rbind(data.frame(Site = "Overall",
                   exc_tss_05 = mean(.$exc_tss_05),
                   exc_tss_50 = mean(.$exc_tss_50),
                   exc_din_1 = mean(.$exc_din_1),
                   exc_din_18 = mean(.$exc_din_18)))
sjPlot::tab_df(exceedances,
               file = "tables/table3.doc",
               col.header = c("Watershed",
                              "> 5 mg/L TSS",
                              "> 50 mg/L TSS",
                              "> 0.1 mg/L DIN",
                              "> 0.18 mg/L DIN"),
               digits = 0)
write_csv(exceedances %>% mutate_if(is.numeric, function(x){round(x,0)}),
          file = "tables/table3.csv")

exceedanceStats = exceedances %>%
  filter(Site != "Overall") %>%
  left_join(water.shed) %>%
  left_join(river.pops)

# Although the frequency by which rivers exceeded thresholds was not significantly associated with the proportion of the watershed cleared of forest, the slope of the relationship between threshold exceedance frequencies and watershed clearing was almost always positive (i.e. β > 0, p > 0.05, except for the 0.18 mg/L DIN threshold). 

exceedanceStats %>%
  select(Site,
         exc_tss_05,
         exc_tss_50,
         exc_din_1,
         exc_din_18,
         Cleared,
         Population) %>%
  pivot_longer(-c(Site, Cleared, Population), 
               names_to = "Threshold", values_to = "Rate") %>%
  ggplot() +
  geom_point(aes(x = Cleared,
                 y = Rate,
                 col = Threshold)) +
  stat_smooth(aes(x = Cleared,
                  y = Rate,
                  col = Threshold),
              method = "lm", se = F) +
  xlab("Cleared land (%)") +
  ylab("Exceedance rate (%)") +
  CJsBasics::BasicTheme

summary(lm(exc_tss_05 ~ Population + Cleared,
           data = exceedanceStats))
summary(lm(exc_tss_50 ~ Population + Cleared,
           data = exceedanceStats))
summary(lm(exc_din_1 ~ Population + Cleared,
           data = exceedanceStats))
summary(lm(exc_din_18 ~ Population + Cleared,
           data = exceedanceStats))


#### End ----
writeLines(capture.output(sessionInfo()), "sessionInfo_results.txt")
