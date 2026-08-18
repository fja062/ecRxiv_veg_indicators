# load packages
library(tidyverse)
library(pwr)
library(lme4)
library(simr)
library(readxl)
library(formattable)
library(sf)
library(WorldFlora)
library(zen4R)


data("WFO.example")
wfo.example.prepared <- WFO.prepare(WFO.example$scientificName)
WFO.match(wfo.example.prepared)

# load in ano data
ano_sp <- st_read("P:/41201785_okologisk_tilstand_2022_2023/data/ANO/naturovervaking_eksport.gdb",
                  layer="ANO_Art", quiet = T)
ano_geo <- st_read("P:/41201785_okologisk_tilstand_2022_2023/data/ANO/naturovervaking_eksport.gdb",
                   layer="ANO_SurveyPoint", quiet = T)

# extract species data
ano_species_distinct <- distinct(ano_sp, art_navn) |> 
  mutate(scientific_name_original = str_replace_all(art_navn, "_", " ")) |> 
  separate_wider_delim(scientific_name_original, delim = " ", names = c("genus", "species", "subspecies"), too_few = "align_start", cols_remove = FALSE) |> 
  mutate(scientific_name_subsp = if_else(
    # select subspecies
    !is.na(subspecies) & !grepl("agg.", scientific_name_original), paste(genus, species, "subsp.", subspecies), 
    # select aggregates
    if_else(!is.na(subspecies) & grepl("agg.", scientific_name_original), paste(genus, species, subspecies),
            # and species identified to species level
            paste(genus, species))
  )) |> 
  mutate(scientific_name_for_matching = paste(genus, species))


# check for incongruencies
ano_sp_prepared <- WFO.prepare(ano_species_distinct$scientific_name_for_matching)


# download WFO plant list from zenodo
download_zenodo("https://doi.org/10.5281/zenodo.20782718", files = "_DwC_backbone_R.zip", path = "/Users/francesca.jaroszynsk/OneDrive - NINA/nina_projects/ANO/lowlands/ecRxiv_veg_indicators/data")

wfo_backbone <- read_delim("C:/Users/francesca.jaroszynsk/OneDrive - NINA/nina_projects/ANO/lowlands/ecRxiv_veg_indicators/data/_DwC_backbone_R/classification.csv", delim = "\t")


# standardise names to the WFO backbone
ano_sp_matched <- WFO.match.fuzzyjoin(spec.data = ano_species_distinct$scientific_name_for_matching, WFO.data = wfo_backbone)

# filter for species with matches in WFO

# 1. new.accepted == TRUE
# 2. old.name is populated
# 3. otherwise scientificName is accepted.



