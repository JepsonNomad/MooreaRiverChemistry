#### This script plots relevant census districts and summarizes population
#### according to overlapping watersheds. Note that census data are private
#### and the summary output for this script is provided in the data subdirectory

#### Load packages ----
library(tidyverse)
library(sf)
library(terra)
source("scripts/000_colorPalettes.R")

#### Import data ----
## > Raster ----
## Island hillshade
hillshade = rast("data/hillshade.tif") %>%
  project("epsg:32706") %>%
  as.data.frame(xy = T) %>%
  as_tibble()

## > Vector ----
## Moorea outline
moorea = read_sf("data/Moorea_LOC_ILE.shp") %>%
  st_transform(32706)

## Watershed boundaries
watersheds = read_sf("data/WaterSample_watersheds.shp") %>%
  rename("Watershed" = Watrshd,
         "Shape_Area" = Shap_Ar,
         "StreamOrder" = StrmOrd,
         "StreamLength" = StrmLng) %>%
  st_transform(32706) %>%
  filter(!(Watershed %in% c("Afareaitu", "Hauru", "Vaiare"))) %>%
  mutate(Watershed = case_when(Watershed == "Haapiti" ~ "Ha'apiti",
                               Watershed == "Maatea" ~ "Ma'atea",
                               TRUE ~ Watershed))

## Moorea rivers
rivers = read_sf("data/Moorea_HYD_TRONCON_EAU_L.shp") %>%
  st_transform(st_crs(moorea))

## Reef crest
reefcrest = read_sf("data/crest.shp") %>%
  st_transform(st_crs(moorea))

## Census boundaries
census_raw = read_sf("data/PRIVATE/ISPFdistrict_shapefile/Districts_MOO_RP17.shp") %>%
  st_transform(st_crs(moorea)) 

#### Data wrangling ----
## Calculate overlap between census tracts and watershed layers
census_sf = census_raw %>%
  st_intersection(watersheds) %>% 
  mutate(intersect_area = st_area(.))
census_area = census_raw %>%
  group_by(IDDIST17) %>%
  summarize(census_ar = st_area(geometry)) %>%
  st_drop_geometry()
census_sel = census_sf %>%
  left_join(census_area) %>%
  select(IDDIST17, Watershed, intersect_area, census_ar) %>%
  mutate(census_rep = as.numeric(intersect_area) / as.numeric(census_ar)) %>%
  filter(census_rep > 0.1)
census_watershed_IDs = census_sel %>%
  select(IDDIST17, Watershed) %>%
  st_drop_geometry() %>%
  filter(!(IDDIST17=="294001650" & Watershed == "Paopao 2"),
         !(IDDIST17=="294001690" & Watershed == "Paopao 1"))
census_watershed_plot = census_raw %>%
  left_join(census_watershed_IDs) %>%
  select(IDDIST17, Watershed) %>%
  filter(!is.na(Watershed))

## Import actual population data
## Have to do this after the wrangling above to join properly
census_numbers = readxl::read_xlsx("data/PRIVATE/ISPFdata/Districts RP2017_Moorea_Univ Californie.xlsx",skip = 2) %>%
  mutate("IDDIST17" = as.character(`...1`)) %>%
  select(IDDIST17,Ensemble) %>%
  right_join(census_watershed_IDs) %>%
  mutate(Site = Watershed)

#### Generate plots ----
censusWatersheds = ggplot() +
  geom_sf(data = moorea) +
  geom_raster(data = hillshade,
              aes(x = x, y = y, fill = hillshade),
              show.legend = F) +
  scale_fill_gradient(low = "black", high = "white") +
  ggnewscale::new_scale_fill() +
  geom_sf(data = moorea,
          fill = "white",
          alpha = 0.2,
          col = "grey40",
          lwd = 1) +
  geom_sf(data = census_watershed_plot,
          aes(fill = Watershed),
          alpha = 0.5) +
  geom_sf(data = watersheds,
          col = "black",
          fill = "transparent",
          lwd = 2)  +
  geom_sf(data = watersheds,
          aes(col = Watershed),
          fill = "transparent",
          lwd = 1) +
  scale_color_manual("Watershed",
                     values = simpsonCols) +
  scale_fill_manual("Census districts",
                    values = simpsonCols) +
  ggtitle("Focal watershed boundaries\nand census districts on Moorea")  +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(0.25,0.25,0.25,0.25,"in"),
        legend.box = "horizontal")
censusWatersheds
ggsave("plots/censusWatersheds.jpg",
       censusWatersheds,
       width = 10, height = 6, units = "in", dpi = 300)

#### Export results ----
write_csv(census_numbers,
          "data/census_by_watershed.csv")

census_numbers %>%
  group_by(Watershed) %>%
  summarise(sum(Ensemble))
