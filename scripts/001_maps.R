#### This script generates a map of Moorea using hillshaded relief, landcover,
#### relevant watershed boundaries, and a reef crest layer. This is Figure 1
#### in Neumann et al 2025 Scientific Reports.


#### Load packages ----
library(tidyverse)
library(sf)
library(terra)
source("scripts/000_colorPalettes.R")

#### Import data ----
## Raster 
## Note, consensusLandUse.tif is on EDI and should be downloaded from there.
lulc = rast("data/consensusLandUse.tif") %>%
  aggregate(4, fun="modal")
hillshade = rast("data/hillshade.tif") %>%
  project(lulc)

lulc
hillshade

## Vector
moorea = read_sf("data/Moorea_LOC_ILE.shp") %>%
  st_transform(st_crs(lulc))
watersheds = read_sf("data/WaterSample_watersheds.shp") %>%
  rename("Watershed" = Watrshd,
         "Shape_Area" = Shap_Ar,
         "StreamOrder" = StrmOrd,
         "StreamLength" = StrmLng) %>%
  st_transform(st_crs(lulc)) %>%
  filter(!(Watershed %in% c("Afareaitu", "Hauru", "Vaiare"))) %>%
  mutate(Watershed = case_when(Watershed == "Haapiti" ~ "Ha'apiti",
                               Watershed == "Maatea" ~ "Ma'atea",
                               TRUE ~ Watershed))
rivers = read_sf("data/Moorea_HYD_TRONCON_EAU_L.shp") %>%
  st_transform(st_crs(lulc))
reefcrest = read_sf("data/crest.shp") %>%
  st_transform(st_crs(lulc))

## Metadata
## attribute codes from EDI
## Download from EDI.
lulc_attr = readxl::read_xlsx(path = "data/metadata_landuse.xlsx",
                              sheet = "Attribute Codes") %>%
  select(Code, Definition)

#### Data wrangling ----
## > Vectors ----
# combine subdivided watersheds
watersheds_aggregated = watersheds %>%
  mutate(grandshed = case_when(Watershed %in% c("Paopao 1",
                                                "Paopao 2",
                                                "Paopao 3") ~ "Paopao",
                               Watershed %in% c("Opunohu 1",
                                                "Opunohu 2") ~ "Opunohu",
                               TRUE ~ Watershed)) %>%
  group_by(grandshed) %>%
  summarize(geometry = st_union(geometry)) %>%
  ungroup()


## > Rasters ----
## Prune raster layers to ROI and aggregate
lulc_mask = mask(lulc, vect(moorea))
hill_mask = mask(hillshade, vect(moorea))
## Convert rasters to data.frames that are ggplot-friendly
lulc_df = as.data.frame(lulc_mask, xy=T) %>%
  tibble() %>%
  left_join(lulc_attr %>%
              rename(classnm = Code)) %>%
  mutate(Definition = factor(Definition,
                             levels = 
                               sort(lulc_attr$Definition)))
hill_df = as.data.frame(hill_mask, xy=T) %>%
  tibble()



#### Plot results ----
## Draw plot
watershedMap = ggplot() +
  geom_sf(data = reefcrest,
          col = "grey80",
          lwd = 0.35) +
  geom_raster(data = lulc_df,
              aes(x = x, y = y, fill = Definition),
              show.legend = FALSE) +
  scale_fill_manual(values = c("#b3fc5cff",
                               "black",
                               "orange4",
                               "#b3fc5cff",
                               "darkgreen",
                               "orange2",
                               "#b3fc5cff",
                               "#b3fc5cff",
                               "beige",
                               "lightblue")) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = hill_df,
              aes(x = x, y = y, fill = hillshade),
              alpha = 0.35,
              show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = "black") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = hill_df,
              aes(x = x, y = y),
              fill = "white",
              alpha = 0.25,
              show.legend = FALSE) +
  geom_sf(data = moorea,
          fill = "transparent",
          col = "grey40",
          lwd = 0.35) +
  geom_sf(data = rivers,
          col = "dodgerblue4",
          lwd = 0.2) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = watersheds,
          col = "black",
          fill = "transparent",
          lwd = 1) +
  geom_sf(data = watersheds,
          aes(col = Watershed,
              fill = Watershed),
          alpha = 0.2,
          lwd = 0.5) +
  scale_color_manual("Watershed", values = simpsonCols) +
  scale_fill_manual("Watershed", values = simpsonCols) +
  ggnewscale::new_scale_colour() +
  geom_segment(aes(x = 190500, xend = 193500,
                   y = 8052200, yend = 8052200),
               col = "black",
               lwd = 1.25) +
  geom_text(label = "3 km",
            aes(x = 192000,
                y = 8052700)) +
  theme_void() +
  theme(legend.key.height = unit(0.05,"in"))


#### Save results ----
ggsave("plots/figure1.jpg",
       watershedMap,
       width = 6, height = 4, units = "in", dpi = 600)
