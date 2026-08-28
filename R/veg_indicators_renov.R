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
library(tidylog)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# download WFO plant list from zenodo
download_zenodo("https://doi.org/10.5281/zenodo.20782718", files = "_DwC_backbone_R.zip", path = "/Users/francesca.jaroszynsk/OneDrive - NINA/nina_projects/ANO/lowlands/ecRxiv_veg_indicators/data")

wfo_backbone <- read_delim("C:/Users/francesca.jaroszynsk/OneDrive - NINA/nina_projects/ANO/lowlands/ecRxiv_veg_indicators/data/_DwC_backbone_R/classification.csv", delim = "\t")



### 7.3. Dataset C
Import plant indicator values from Tyler et al. (2021):
  ```{r}
ind.Tyler <- readRDS("P:/41201785_okologisk_tilstand_2022_2023/data/functional plant indicators/ind.Tyler.RDS")
```

The Swedish plant indicator values dataset published by Tyler et al. (2021) contains a large collection of plant indicators based on the Swedish flora, which is well representative of the Norwegian flora as well. From this set, we use indicator data for moisture and Moisture as these are thought to be subject to potential change due to ongoing pressures in the respective ecosystems (see details above under 3.4 'Impact factors').


## 8. Spatial units
The spatial units for functional plant community indicators are governed by the ANO data, which consist of 1000 randomly chosen sites in Norway. Each site is a 500 x 500 m grid cell with 18 monitoring points of 250sqm and a central 1 x 1 m vegetation plot. The basic unit for which observational community data exist, and for which every functional plant community indicator is computed, is this central 1 x 1 m vegetation plot. Aggregation of these basic spatial units to higher ones - like sites, municipalities, counties, regions, or the national level - should be done with appropriate consideration of the spatial structure of the data, i.e. imbalances and clustering (highly unequal number of points of an ecosystem type between sites).




## 9. Analyses

<!--# 
  
  Use this header for documenting the analyses. Put code in separate code chunks, and annotate the code in between using normal text (i.e. between the chunks, and try to avoid too many hashed out comments inside the code chunks). Add subheaders as needed. 

Code folding is activated, meaning the code will be hidden by default in the html (one can click to expand it).

Caching is also activated (from the top YAML), meaning that rendering to html will be quicker the second time you do it. This will create a folder inside you project folder (called INDICATORID_cache). Sometimes caching created problems because some operations are not rerun when they should be rerun. Try deleting the cash folder and try again.

-->
  
  #### Data handling
  - Checking for errors
- Checking species nomenclature in the different species lists to make species and indicator data possible to merge
- Merging indicator data with monitoring data and indicator data with reference data
(not shown here, but documented in the code)

```{r dataHandling, results='hide', warning=F, message=F}

#### Plant indicator data
ind_tyler <- ind.Tyler |>
  rename(scientific_name = Scientific_name) |> 
  mutate(scientific_name_original = scientific_name,
         scientific_name = str_replace_all(scientific_name, "ssp.", "subsp."),
         scientific_name = str_replace_all(scientific_name, " x ", " \u00D7 ")) |> 
  filter(!is.na(scientific_name), !scientific_name == "") |> 
  distinct(scientific_name_original, scientific_name, Moisture, Nitrogen) |> 
  tibble()

tyler_prepared_wfo <- WFO.prepare(ind_tyler$scientific_name)

# 
tyler_prepared <- tyler_prepared_wfo |>
  tibble() |> 
  mutate(
    # first word of Authorship
    first_author = if_else(!is.na(Authorship), word(Authorship, 1, 1), ""),
    
    # 1) If first_author looks like a missed subspecies epithet (lowercase),
    #    AND the name is NOT a hybrid, append "subsp. <first_author>"
    spec.name = if_else(
      first_author != "" &
        str_detect(first_author, "^[a-z]") &               # starts lowercase
        !str_detect(spec.name, " x$") &                    # not ending with " x"
        !str_detect(spec.name, "\u00D7$") &                # not ending with "×"
        !str_detect(spec.name, "^\\s*[x\u00D7]\\b"),       # not starting with x/× hybrid marker
      paste(spec.name, "subsp.", first_author),
      spec.name
    ),
    
    # 2) Hybrids: append first_author after trailing " x" or "×"
    spec.name = if_else(
      str_detect(spec.name, " x$"),
      paste0(spec.name, " ", first_author),
      spec.name
    ),
    spec.name = if_else(
      str_detect(spec.name, "\u00D7$"),
      paste0(spec.name, " ", first_author),
      spec.name
    ),
    spec.name = if_else(
      str_detect(Authorship, "^\u00D7"),
      paste0(spec.name, Authorship),
      spec.name
    ),
    
    # 3) Clean sect. and whitespace
    spec.name = spec.name |>
      str_replace_all("\\bsect\\b\\.?", " ") |>
      str_squish()
  ) |>
  rename(clean_string = spec.name) |> 
  distinct(spec.full, clean_string)



# standardise names to the WFO backbone
tyler_sp_matched <- WFO.match(spec.data = tyler_prepared$clean_string,
                             WFO.data = wfo_backbone,
                             Fuzzy = 0.15,
                             Fuzzy.max = 50,
                             Fuzzy.one = FALSE)



# create accepted name column according to latest taxonomical nomenclature
tyler_sp_clean <- tyler_sp_matched |>
  # First, store the original string clearly
  rename(clean_string = spec.name.ORIG) |>
  group_by(clean_string) |>
  summarise(
    # 1) Flag multiple scientificName suggestions
    flag_multiple_suggestions = n_distinct(scientificName) > 1,
    
    # 2) Candidate accepted_name from Old.name when possible
    accepted_name = case_when(
      any(New.accepted == TRUE & Old.name != "") ~ 
        # take one Old.name where New.accepted == TRUE and Old.name non-empty
        Old.name[New.accepted == TRUE & Old.name != ""][1],
      TRUE ~ 
        # otherwise fall back to (one) scientificName
        scientificName[1]
    ),
    
    # 3) Where did accepted_name come from?
    accepted_from = case_when(
      any(New.accepted == TRUE & Old.name != "") ~ "Old.name",
      TRUE ~ "scientificName"
    ),
    .groups = "drop"
  )


# check the species that change name where many options were available
tyler_sp_clean |> filter(clean_string != accepted_name) |> view()
tyler_sp_clean |> filter(flag_multiple_suggestions == TRUE, clean_string != accepted_name) |> view()


# correct incorrect corrections. haha
#tyler_sp_clean <- tyler_sp_clean |> 
#  mutate(flag_species_revert = case_when(
#    flag_multiple_suggestions == TRUE & clean_string != accepted_name ~ paste("reverted to #original from", accepted_name),
#    TRUE ~ ""
#  ),
#  accepted_name = case_when(
#    flag_multiple_suggestions == TRUE & clean_string != accepted_name ~ clean_string,
#    TRUE ~ accepted_name
#  ))



# bind new species names onto indicator dataset
ind_tyler <- ind_tyler |> 
  left_join(tyler_sp_matched_for_join, by = "clean_string") |> 
  # filter out sect. species and subspecies
  #filter(!grepl("subsp.", scientific_name_original)) |> 
  distinct()


ind_tyler |> filter(is.na(accepted_name))


ind_tyler |> 
  group_by(accepted_name) |> 
  filter(n() > 1)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# GRUK

# read in data
GRUK_species <- read_excel("P:/41201785_okologisk_tilstand_2022_2023/data/GRUK/GRUK_alle_artsdata_2020-24.xlsx", sheet="Arter i ruter")
GRUK_ruter <- read_excel("P:/41201785_okologisk_tilstand_2022_2023/data/GRUK/GRUK_alle_artsdata_2020-24.xlsx", sheet="Ruter")
GRUK_sirkler <- read_excel("P:/41201785_okologisk_tilstand_2022_2023/data/GRUK/GRUK_alle_artsdata_2020-24.xlsx", sheet="Sirkler")
GRUK_polygoner <- read_excel("P:/41201785_okologisk_tilstand_2022_2023/data/GRUK/GRUK_alle_artsdata_2020-24.xlsx", sheet="Polygoner")


## GRUK species data handling
GRUK_species <- GRUK_species |> 
  rename(norsk_navn = `Norsk navn`,
         scientific_name = `Latinsk navn`,
         art_dekning = `Dekning %`)


# fix species names
GRUK_species <- GRUK_species |>
  mutate(
    scientific_name_original = scientific_name,
    scientific_name = str_replace_all(scientific_name, "ssp\\.", "subsp."),
    scientific_name = str_replace_all(scientific_name, " x ", " \u00D7 ")
  ) |>
  filter(!is.na(scientific_name), scientific_name != "") |>
  tibble()

GRUK_prepared_wfo <- WFO.prepare(GRUK_species$scientific_name)

GRUK_prepared <- GRUK_prepared_wfo |>
  mutate(
    subspecies_epithet = if_else(!is.na(Authorship), word(Authorship, 1, 1), ""),
    
    # Attach subspecies epithet only if:
    #   - we have an epithet
    #   - AND Authorship does not contain "agg."
    spec.name = if_else(
      subspecies_epithet != "" & !str_detect(Authorship, "\\bagg\\.?"),
      paste(spec.name, "subsp.", subspecies_epithet),
      spec.name
    ),
   # # attach agg. info
   # spec.name = if_else(
   #   str_detect(Authorship, "\\bagg\\.?") & !str_detect(spec.name, "\\bagg\\.?"),
   #   paste(spec.name, subspecies_epithet),  # just add "agg." as a trailing word
   #   spec.name
   # ),
    
    # Handle hybrids / ×
    spec.name = if_else(grepl(" x$", spec.name), paste(spec.name, subspecies_epithet), spec.name),
    spec.name = if_else(grepl("\u00D7$", spec.name), paste(spec.name, subspecies_epithet), spec.name),
    spec.name = if_else(grepl("^\u00D7", Authorship), paste(spec.name, subspecies_epithet), spec.name),
    
    # Clean out 'sect.'
    spec.name = spec.name |>
      str_replace_all("\\bsect\\b\\.?", " ") |>
      str_squish() |> 
     str_to_sentence()
  ) |>
  rename(clean_string = spec.name) |>
  distinct(spec.full, clean_string)


# standardise names to the WFO backbone
GRUK_sp_matched <- WFO.match(spec.data = GRUK_prepared$clean_string,
                             WFO.data = wfo_backbone,
                             Fuzzy = 0.15,
                             Fuzzy.max = 50,
                             Fuzzy.one = FALSE)


# create accepted name column according to latest taxonomical nomenclature
GRUK_sp_clean <- GRUK_sp_matched |>
  # First, store the original string clearly
  rename(clean_string = spec.name.ORIG) |>
  group_by(clean_string) |>
  summarise(
    # 1) Flag multiple scientificName suggestions
    flag_multiple_suggestions = n_distinct(scientificName) > 1,
    
    # 2) Candidate accepted_name from Old.name when possible
    accepted_name = case_when(
      any(New.accepted == TRUE & Old.name != "") ~ 
        # take one Old.name where New.accepted == TRUE and Old.name non-empty
        Old.name[New.accepted == TRUE & Old.name != ""][1],
      TRUE ~ 
        # otherwise fall back to (one) scientificName
        scientificName[1]
    ),
    
    # 3) Where did accepted_name come from?
    accepted_from = case_when(
      any(New.accepted == TRUE & Old.name != "") ~ "Old.name",
      TRUE ~ "scientificName"
    ),
    .groups = "drop"
  )


# check the species that change name where many options were available
GRUK_sp_clean |> filter(clean_string != accepted_name)
GRUK_sp_clean |> filter(flag_multiple_suggestions == TRUE, clean_string != accepted_name)


# correct incorrect corrections. haha
GRUK_sp_clean <- GRUK_sp_clean |> 
  mutate(  flag_species_revert = case_when(
    flag_multiple_suggestions == TRUE & clean_string != accepted_name ~ paste("reverted to original from", accepted_name),
    TRUE ~ ""
  ),
         accepted_name = case_when(
    flag_multiple_suggestions == TRUE & clean_string != accepted_name ~ clean_string,
    TRUE ~ accepted_name
  ))


# bind new species names onto original dataset
GRUK_species_clean <- left_join(GRUK_prepared, GRUK_sp_clean, by = "clean_string") |> 
  full_join(GRUK_species, by = join_by(spec.full == scientific_name_original)) |> 
  # filter out sect. species and subspecies
  #filter(!grepl("subsp.", scientific_name_original)) |> 
  distinct()


  
GRUK_species_clean |> 
  filter(is.na(accepted_name)) |> 
  tibble()





# merge species data with indicators
GRUK.species.ind <- merge(x=GRUK.species[,c("Species", "art_dekning", "ParentGlobalID","PolygonID","RuteID")], 
                          y= ind.dat[,c("species","CC", "SS", "RR","Light", "Nitrogen", "Soil_disturbance")],
                          by.x="Species", by.y="species", all.x=T)
summary(GRUK.species.ind)

# checking which species didn't find a match
unique(GRUK.species.ind[is.na(GRUK.species.ind$Light & 
                                is.na(GRUK.species.ind$RR)),'Species'])

# fix species name issues
ind.dat <- ind.dat %>% 
  #  mutate(species=str_replace(species,"Aconitum lycoctonum", "Aconitum septentrionale")) %>% 
  #  mutate(species=str_replace(species,"Carex simpliciuscula", "Kobresia simpliciuscula")) %>%
  #  mutate(species=str_replace(species,"Carex myosuroides", "Kobresia myosuroides")) %>%
  #  mutate(species=str_replace(species,"Artemisia rupestris", "Artemisia norvegica")) %>%
  mutate(species=str_replace(species,"Cotoneaster simonsii", "Cotoneaster symondsii")) %>%
  mutate(species=str_replace(species,"Rosa vosagica", "Rosa vosagiaca"))

GRUK.species <- GRUK.species %>%
  mutate(species=str_replace(Species,"Acinos arvensis", "Clinopodium acinos")) %>%
  mutate(Species=str_replace(Species,"Arabis wahlenbergii", "Arabis hirsuta")) %>%
  #  mutate(Species=str_replace(Species,"Arctous alpinus", "Arctous alpina")) %>%
  #  mutate(Species=str_replace(Species,"Betula tortuosa", "Betula pubescens")) %>%
  #  mutate(Species=str_replace(Species,"Blysmopsis rufa", "Blysmus rufus")) %>%
  #  mutate(Species=str_replace(Species,"Cardamine nymanii", "Cardamine pratensis")) %>%
  #  mutate(Species=str_replace(Species,"Carex adelostoma", "Carex buxbaumii")) %>%
  #  mutate(Species=str_replace(Species,"Carex leersii", "Carex echinata")) %>%
  mutate(Species=str_replace(Species,"Carex paupercula", "Carex magellanica")) %>%
  #  mutate(Species=str_replace(Species,"Carex simpliciuscula", "Kobresia simpliciuscula")) %>%
  mutate(Species=str_replace(Species,"Carex viridula", "Carex flava")) %>%
  #  mutate(Species=str_replace(Species,"Chamaepericlymenum suecicum", "Cornus suecia")) %>%
  #  mutate(Species=str_replace(Species,"Cicerbita alpina", "Lactuca alpina")) %>%
  mutate(Species=str_replace(Species,"Cotoneaster scandinavicus", "Cotoneaster integerrimus")) %>%
  # mutate(Species=str_replace(Species,"Cotoneaster symondsii", "Cotoneaster integrifolius")) %>%
  mutate(Species=str_replace(Species,"Cyanus montanus", "Centaurea montana")) %>%
  #  mutate(Species=str_replace(Species,"Empetrum hermaphroditum", "Empetrum nigrum")) %>%
  mutate(Species=str_replace(Species,"Erysimum virgatum", "Erysimum strictum")) %>%
  #  mutate(Species=str_replace(Species,"Festuca prolifera", "Festuca rubra")) %>%
  mutate(Species=str_replace(Species,"Festuca trachyphylla", "Festuca brevipila")) %>%
  mutate(Species=str_replace(Species,"Galium album", "Galium mollugo")) %>%
  #  mutate(Species=str_replace(Species,"Galium elongatum", "Galium palustre")) %>%
  mutate(Species=str_replace(Species,"Helictotrichon pratense", "Avenula pratensis")) %>%
  mutate(Species=str_replace(Species,"Helictotrichon pubescens", "Avenula pubescens")) %>%
  #  mutate(Species=str_replace(Species,"Hieracium alpina", "Hieracium Alpina")) %>%
  #  mutate(Species=str_replace(Species,"Hieracium alpinum", "Hieracium Alpina")) %>%
  #  mutate(Species=str_replace(Species,"Hieracium hieracium", "Hieracium Hieracium")) %>%
  #  mutate(Species=str_replace(Species,"Hieracium hieracioides", "Hieracium umbellatum")) %>%
  mutate(Species=str_replace(Species,"Hieracium murorum", "Hieracium Hieracium")) %>%
  mutate(Species=str_replace(Species,"Hieracium vulgatum", "Hieracium Vulgata")) %>%
  #  mutate(Species=str_replace(Species,"Hieracium oreadea", "Hieracium Oreadea")) %>%
  #  mutate(Species=str_replace(Species,"Hieracium prenanthoidea", "Hieracium Prenanthoidea")) %>%
  #  mutate(Species=str_replace(Species,"Hieracium vulgata", "Hieracium Vulgata")) %>%
  #  mutate(Species=str_replace(Species,"Hieracium pilosella", "Pilosella officinarum")) %>%
  #  mutate(Species=str_replace(Species,"Hieracium vulgatum", "Hieracium umbellatum")) %>%
  #  mutate(Species=str_replace(Species,"Hierochloã« alpina", "Hierochloë alpina")) %>%
  #  mutate(Species=str_replace(Species,"Hierochloã« hirta", "Hierochloë hirta")) %>%
  #  mutate(Species=str_replace(Species,"Hierochloã« odorata", "Hierochloë odorata")) %>%
  mutate(Species=str_replace(Species,"Hylotelephium maximum", "Sedum telephium")) %>%
  #  mutate(Species=str_replace(Species,"Listera cordata", "Neottia cordata")) %>%
  #  mutate(Species=str_replace(Species,"Leontodon autumnalis", "Scorzoneroides autumnalis")) %>%
  mutate(Species=str_replace(Species,"Lepidotheca suaveolens", "Matricaria discoidea")) %>%
  #  mutate(Species=str_replace(Species,"Loiseleuria procumbens", "Kalmia procumbens")) %>%
  mutate(Species=str_replace(Species,"Malus ×domestica", "Malus domestica")) %>%
  #  mutate(Species=str_replace(Species,"Mycelis muralis", "Lactuca muralis")) %>%
  #  mutate(Species=str_replace(Species,"Omalotheca supina", "Gnaphalium supinum")) %>%
  #  mutate(Species=str_replace(Species,"Omalotheca norvegica", "Gnaphalium norvegicum")) %>%
  #  mutate(Species=str_replace(Species,"Omalotheca sylvatica", "Gnaphalium sylvaticum")) %>%
  #  mutate(Species=str_replace(Species,"Oreopteris limbosperma", "Thelypteris limbosperma")) %>%
  #  mutate(Species=str_replace(Species,"Oxycoccus microcarpus", "Vaccinium microcarpum")) %>%
  #  mutate(Species=str_replace(Species,"Oxycoccus palustris", "Vaccinium oxycoccos")) %>%
  #  mutate(Species=str_replace(Species,"Phalaris minor", "Phalaris arundinacea")) %>%
  mutate(Species=str_replace(Species,"Pilosella peletariana", "Pilosella officinarum")) %>%
  #  mutate(Species=str_replace(Species,"Pinus unicinata", "Pinus mugo")) %>%
  #  mutate(Species=str_replace(Species,"Poa alpigena", "Poa pratensis")) %>%
  mutate(Species=str_replace(Species,"Poa angustifolia", "Poa pratensis")) %>%
  mutate(Species=str_replace(Species,"Poa humilis", "Poa pratensis")) %>%
  #  mutate(Species=str_replace(Species,"Pyrola grandiflora", "Pyrola rotundifolia")) %>%
  mutate(Species=str_replace(Species,"Rosa dumalis", "Rosa vosagiaca")) %>%
  #  mutate(Species=str_replace(Species,"Rumex alpestris", "Rumex acetosa")) %>%
  mutate(Species=str_replace(Species,"Sorbus hybrida", "Hedlundia hybrida")) %>%
  mutate(Species=str_replace(Species,"Spergularia salina", "Spergularia marina")) %>%
  #  mutate(Species=str_replace(Species,"Syringa emodi", "Syringa vulgaris")) %>%
  #  mutate(Species=str_replace(Species,"Taraxacum crocea", "Taraxacum officinale")) %>%
  #  mutate(Species=str_replace(Species,"Taraxacum croceum", "Taraxacum officinale")) %>%
  #  mutate(Species=str_replace(Species,"Trientalis europaea", "Lysimachia europaea")) %>%
  mutate(Species=str_replace(Species,"Trifolium pallidum", "Trifolium pratense"))

# merge species data with indicators
GRUK.species.ind <- merge(x=GRUK.species[,c("Species", "art_dekning", "ParentGlobalID","PolygonID","RuteID")], 
                          y= ind.dat[,c("species","CC", "SS", "RR","Light", "Nitrogen", "Soil_disturbance")],
                          by.x="Species", by.y="species", all.x=T)
summary(GRUK.species.ind)
# checking which species didn't find a match
unique(GRUK.species.ind[is.na(GRUK.species.ind$Light & 
                                is.na(GRUK.species.ind$RR)),'Species'])

## GRUK ruter data handling
names(GRUK.ruter)
# make coordinates numeric
GRUK.ruter <- GRUK.ruter %>% 
  mutate( UTM33_E_ne=as.numeric(UTM33_E_ne) ) %>%
  mutate( UTM33_N_ne=as.numeric(UTM33_N_ne) ) %>%
  mutate( UTM33_E_sw=as.numeric(UTM33_E_sw) ) %>%
  mutate( UTM33_N_sw=as.numeric(UTM33_N_sw) )

# calculate central coordinates for each plot
GRUK.ruter <- GRUK.ruter %>% 
  mutate(UTM33_N = (UTM33_N_ne + UTM33_N_sw)/2) %>%
  mutate(UTM33_E = (UTM33_E_ne + UTM33_E_sw)/2)

# some of the calculations throw NA's because there's only one set of coordinates, coalesce that set into the calculation column 
GRUK.ruter <- GRUK.ruter %>% 
  mutate (UTM33_N = coalesce(UTM33_N,UTM33_N_ne) ) %>%
  mutate (UTM33_E = coalesce(UTM33_E,UTM33_E_ne) )
GRUK.ruter <- GRUK.ruter %>% 
  mutate (UTM33_N = coalesce(UTM33_N,UTM33_N_sw) ) %>%
  mutate (UTM33_E = coalesce(UTM33_E,UTM33_E_sw) )


## GRUK sirkler data handling

## GRUK polygoner data handling

## merge information on mapping units and condition variables from GRUK.sirkler into GRUK.ruter
names(GRUK.ruter)
names(GRUK.sirkler)
GRUK.variables <- merge(x=GRUK.ruter[,c(2,4:11,19:23,25,52:53)], 
                        y=GRUK.sirkler[,c("GlobalID",
                                          "Kartleggingsenhet 1:5000",
                                          "Spor etter slitasje og slitasjebetinget erosjon (%)",
                                          "Dekning % av nakent berg",
                                          "Total dekning % av vedplanter i feltsjikt",
                                          "Dekning % av busker i busksjikt",
                                          "Dekning % av tresjikt",
                                          "Dekning % av problemarter",
                                          "Total dekning % av fremmede arter")], 
                        by.x="GlobalID", by.y="GlobalID", all.x=T)
summary(GRUK.variables)

## merge information on condition and quality from GRUK.polygoner into GRUK.variables
# transform GRUK.variables into spatial object
GRUK.variables <- st_as_sf(GRUK.variables, coords = c("UTM33_E","UTM33_N"),remove=F, crs = 25833)

# transform GRUK.polygoner into spatial object
GRUK.polygoner <- st_as_sf(GRUK.polygoner, wkt = "WKT" ,remove=F, crs = 25833)

tm_shape(GRUK.polygoner) +
  tm_graticules() +
  tm_polygons("PolygonID") +
  tm_shape(GRUK.variables) +
  tm_dots("RuteID")

# run a spatial join to get columns from GRUK.polygoner into GRUK.variables
#GRUK.variables <- st_join(GRUK.variables,GRUK.polygoner[,c(3:4,9,15,18,20,22,24,60)])
#names(GRUK.variables)[1:33]<-c("GlobalID","PolygonID.x","RuteID","RuteID_loknr",
#                               "Dekning_karplanter_feltsjikt","Dekning_moser","Dekning_lav","Dekning_strø",
#                               "Dekning_bar_substrat","Precision","UTM33_E_ne","UTM33_N_ne",
#                               "UTM33_E_sw","UTM33_N_sw","areal(m2)","UTM33_N","UTM33_E","Kartleggingsenhet_1til5000",
#                               "erosjon_prosent","Dekning_nakentberg",
#                               "Totaldekning_vedplanter_feltsjikt","Dekning_busker_busksjikt","Dekning_tresjikt",
#                               "Dekning_problemarter","Totaldekning_fremmedearter","LokalitetID","PolygonID.y",
#                               "Kartleggingsdato","Lokalitetskvalitet","Kommune","Tilstand","Naturmangfold","NiNKartleggingsenheter")

# check how good the spatial join worked
#cbind(GRUK.variables$PolygonID.x,GRUK.variables$PolygonID.y)
#GRUK.variables[7,]
#GRUK.polygoner[GRUK.polygoner$PolygonID=="46-2",]
# some points could not be matched to polygons -> merge by PolygonID instead, drop geometry of GRUK.variables first
GRUK.variables <- st_drop_geometry(GRUK.variables)
names(GRUK.variables)
names(GRUK.polygoner)
GRUK.variables <- merge(x=GRUK.variables, 
                        y=GRUK.polygoner[,c(3:4,9,15,18,20,22,24)], 
                        by.x="PolygonID", by.y="PolygonID", all.x=T)
summary(GRUK.variables) 
summary(as.factor(GRUK.variables$Tilstand)) # no unexpected NA's

# edit the column names
names(GRUK.variables)[1:32]<-c("PolygonID","GlobalID","RuteID","RuteID_loknr",
                               "Dekning_karplanter_feltsjikt","Dekning_moser",
                               "Dekning_lav","Dekning_strø",
                               "Dekning_bar_substrat","Precision",
                               "UTM33_E_ne","UTM33_N_ne","UTM33_E_sw","UTM33_N_sw",
                               "areal_m2","UTM33_N","UTM33_E","Kartleggingsenhet_1til5000",
                               "erosjon_prosent","Dekning_nakentberg",
                               "Totaldekning_vedplanter_feltsjikt","Dekning_busker_busksjikt",
                               "Dekning_tresjikt","Dekning_problemarter","Totaldekning_fremmedearter",
                               "LokalitetID","Aar","Lokalitetskvalitet","Kommune","Tilstand",
                               "Naturmangfold","NiNKartleggingsenheter")



## adding information on ecosystem and condition variables to species+indicator data
names(GRUK.species.ind)
names(GRUK.variables)
GRUK.species.ind <- merge(x=GRUK.species.ind, 
                          y=GRUK.variables[,-c(1,3)], 
                          by.x="ParentGlobalID", by.y="GlobalID", all.x=T)
summary(GRUK.species.ind)
# fixing variable types
GRUK.species.ind <- GRUK.species.ind %>% 
  mutate(Species = as.factor(Species)) %>%
  mutate(areal_m2 = as.numeric(areal_m2)) %>%
  mutate(Kartleggingsenhet_1til5000 = as.factor(Kartleggingsenhet_1til5000)) %>%
  mutate(Dekning_nakentberg = as.numeric(Dekning_nakentberg)) %>%
  mutate(Dekning_problemarter = as.numeric(Dekning_problemarter)) %>%
  mutate(Lokalitetskvalitet = as.factor(Lokalitetskvalitet)) %>%
  mutate(Kommune = as.factor(Kommune)) %>%
  mutate(Tilstand = as.factor(Tilstand)) %>%
  mutate(Naturmangfold = as.factor(Naturmangfold)) %>%
  mutate(NiNKartleggingsenheter = as.factor(NiNKartleggingsenheter))
summary(GRUK.species.ind)

# trimming away the points without information on NiN, species or cover  
GRUK.species.ind <- GRUK.species.ind[!is.na(GRUK.species.ind$Species),]
GRUK.species.ind <- GRUK.species.ind[!is.na(GRUK.species.ind$art_dekning),]
# no NA's for kartleggingsenhet

#rm(GRUK.species)
#rm(GRUK.ruter)


summary(GRUK.species.ind)
head(GRUK.species.ind)




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### double check the Coa and Aosa species in the ANO fuzzy matching

#### ANO monitoring data



# load in ano data
ano_sp <- st_read("P:/41201785_okologisk_tilstand_2022_2023/data/ANO/naturovervaking_eksport.gdb",
                  layer="ANO_Art", quiet = T)
ano_geo <- st_read("P:/41201785_okologisk_tilstand_2022_2023/data/ANO/naturovervaking_eksport.gdb",
                   layer="ANO_SurveyPoint", quiet = T)

# extract species data
ano_sp_processed <- ano_sp |> 
  tibble() |> 
  mutate(scientific_name_original = str_replace_all(art_navn, "_", " ")) |> 
  separate_wider_delim(scientific_name_original, delim = " ", names = c("genus", "species", "subspecies"), too_few = "align_start", cols_remove = FALSE) |> 
  mutate(scientific_name_for_matching = if_else(
    # select subspecies
    !is.na(subspecies) & !grepl("agg.", scientific_name_original) & !grepl("sp.$", scientific_name_original), paste(genus, species, "subsp.", subspecies), 
    # select aggregates 
    # CURRENTLY EXCLUDING AGGREGATE
    #if_else(!is.na(subspecies) & grepl("agg.", scientific_name_original), paste(genus, species, subspecies),
            # and species identified to species level
            paste(genus, species)#)
  ),
  scientific_name_for_matching = if_else(grepl("^sp.$", species), paste0(genus), scientific_name_for_matching)) |> 
  # capitalise first letter of genus
  mutate(scientific_name_for_matching = str_to_sentence(scientific_name_for_matching))

ano_species_distinct <- distinct(ano_sp_processed, art_navn, .keep_all = TRUE)

ano_species_distinct |> 
  group_by(scientific_name_for_matching) |> 
    filter(n() > 1)

# check for incongruencies
ano_sp_prepared <- WFO.prepare(ano_species_distinct$scientific_name_for_matching)


# standardise names to the WFO backbone
ano_sp_matched <- WFO.match(spec.data = ano_species_distinct$scientific_name_for_matching, WFO.data = wfo_backbone, Fuzzy = 2)

# filter for species with matches in WFO
# 1. new.accepted == TRUE
# 2. old.name is populated
# 3. otherwise scientificName is accepted.

ano_sp_matched_accepted <- ano_sp_matched |> 
  mutate(accepted_name = if_else(
    New.accepted == TRUE & !Old.name == "", Old.name, scientificName
  )) |> 
  distinct(accepted_name, spec.name.ORIG)


# run through backbone clean
ano_xtras <- ano_species_distinct |> 
  semi_join(ano_sp_matched_accepted |> 
                                    group_by(accepted_name) |> 
                                    filter(n() > 1), by = join_by(scientific_name_for_matching == spec.name.ORIG)) |> 
  mutate(scientific_name_for_matching_edit = str_remove(scientific_name_for_matching, "subsp\\..*")) 

ano_xtras_matched <- WFO.match.fuzzyjoin(spec.data = ano_xtras$scientific_name_for_matching_edit, WFO.data = wfo_backbone)

result <- ano_xtras_matched |>
  group_by(spec.name.ORIG) |>
  mutate(
    has_accepted_old = any(
      New.accepted == TRUE &
        !is.na(Old.name) &
        Old.name != ""
    ),
    
    scientific_candidates = list(unique(
      scientificName[
        New.accepted == FALSE &
          (is.na(Old.name) | Old.name == "")
      ]
    )),
    
    n_scientific_candidates = length(scientific_candidates[[1]]),
    
    multiple_scientific_matches =
      !has_accepted_old &
      n_scientific_candidates > 1,
    
    accepted_name = if_else(
      has_accepted_old,
      
      Old.name[
        New.accepted == TRUE &
          !is.na(Old.name) &
          Old.name != ""
      ][1],
      
      scientificName[
        New.accepted == FALSE &
          (is.na(Old.name) | Old.name == "")
      ][1]
    )
  ) |>
  reframe(
    accepted_name = first(accepted_name),
    scientific_candidates = first(scientific_candidates),
    n_scientific_candidates = first(n_scientific_candidates),
    multiple_scientific_matches = first(multiple_scientific_matches)
  )


ano_xtras

# join misfits back onto ano accepted list
ano_sp_matched_accepted <- ano_sp_matched_accepted |> 
  filter(!is.na(accepted_name)) |> 
  bind_rows(ano_xtras_matched_accepted) |> 
  tibble()


# fit back onto the ano sp data
ano_species_distinct |> 
  full_join(ano_sp_matched_accepted, by = join_by(scientific_name_for_matching == spec.name.ORIG))


















#head(ANO.sp)
#head(ANO.geo)

## fix NiN information
ANO.geo$hovedtype_rute <- substr(ANO.geo$kartleggingsenhet_1m2,1,3) # take the 3 first characters
ANO.geo$hovedtype_rute <- gsub("-", "", ANO.geo$hovedtype_rute) # remove hyphen
unique(as.factor(ANO.geo$hovedtype_rute))

#ANO.geo$hovedoekosystem_rute <- ANO.geo$hovedtype_rute
ANO.geo <- ANO.geo %>% mutate(hovedoekosystem_rute=recode(hovedtype_rute, 
                                                          "T4"="Forest", "T30"="Forest",
                                                          "T3"="Mountain", "T7"="Mountain", "T14"="Mountain", "T22"="Mountain",
                                                          "V1"="Wetland", "V2"="Wetland", "V3"="Wetland", "V4"="Wetland", "V5"="Wetland", "V6"="Wetland", "V7"="Wetland", "V8"="Wetland", 
                                                          "T31"="Seminat", "T32"="Seminat", "T33"="Seminat", "T34"="Seminat","V9"="Seminat", "V10"="Seminat",
                                                          "T2"="Natopen", "T8"="Natopen", "T11"="Natopen", "T12"="Natopen","T13"="Natopen", "T15"="Natopen","T16"="Natopen", "T18"="Natopen","T21"="Natopen", "T24"="Natopen", "T29"="Natopen"
))
unique(as.factor(ANO.geo$hovedoekosystem_rute))

## fix NiN-variables
colnames(ANO.geo)
colnames(ANO.geo)[42:47] <- c("groeftingsintensitet",
                              "bruksintensitet",
                              "beitetrykk",
                              "slatteintensitet",
                              "tungekjoretoy",
                              "slitasje")
#head(ANO.geo)

# remove variable code in the data
ANO.geo$groeftingsintensitet <- gsub("7GR-GI_", "", ANO.geo$groeftingsintensitet) 
unique(ANO.geo$groeftingsintensitet)
ANO.geo$groeftingsintensitet <- gsub("X", "NA", ANO.geo$groeftingsintensitet)
unique(ANO.geo$groeftingsintensitet)
ANO.geo$groeftingsintensitet <- as.numeric(ANO.geo$groeftingsintensitet)
unique(ANO.geo$groeftingsintensitet)

ANO.geo$bruksintensitet <- gsub("7JB-BA_", "", ANO.geo$bruksintensitet) 
unique(ANO.geo$bruksintensitet)
ANO.geo$bruksintensitet <- gsub("X", "NA", ANO.geo$bruksintensitet)
unique(ANO.geo$bruksintensitet)
ANO.geo$bruksintensitet <- as.numeric(ANO.geo$bruksintensitet)
unique(ANO.geo$bruksintensitet)

ANO.geo$beitetrykk <- gsub("7JB-BT_", "", ANO.geo$beitetrykk) 
unique(ANO.geo$beitetrykk)
ANO.geo$beitetrykk <- gsub("X", "NA", ANO.geo$beitetrykk)
unique(ANO.geo$beitetrykk)
ANO.geo$beitetrykk <- as.numeric(ANO.geo$beitetrykk)
unique(ANO.geo$beitetrykk)

ANO.geo$slatteintensitet <- gsub("7JB-SI_", "", ANO.geo$slatteintensitet) 
unique(ANO.geo$slatteintensitet)
ANO.geo$slatteintensitet <- gsub("X", "NA", ANO.geo$slatteintensitet)
unique(ANO.geo$slatteintensitet)
ANO.geo$slatteintensitet <- as.numeric(ANO.geo$slatteintensitet)
unique(ANO.geo$slatteintensitet)

ANO.geo$tungekjoretoy <- gsub("7TK_", "", ANO.geo$tungekjoretoy) 
unique(ANO.geo$tungekjoretoy)
ANO.geo$tungekjoretoy <- gsub("X", "NA", ANO.geo$tungekjoretoy)
unique(ANO.geo$tungekjoretoy)
ANO.geo$tungekjoretoy <- as.numeric(ANO.geo$tungekjoretoy)
unique(ANO.geo$tungekjoretoy)

ANO.geo$slitasje <- gsub("7SE_", "", ANO.geo$slitasje) 
unique(ANO.geo$slitasje)
ANO.geo$slitasje <- gsub("X", "NA", ANO.geo$slitasje)
unique(ANO.geo$slitasje)
ANO.geo$slitasje <- as.numeric(ANO.geo$slitasje)
unique(ANO.geo$slitasje)

## check that every point is present only once
#length(levels(as.factor(ANO.geo$ano_flate_id)))
#length(levels(as.factor(ANO.geo$ano_punkt_id)))
summary(as.factor(ANO.geo$ano_punkt_id))
# there's many double presences, probably some wrong registrations of point numbers

# we filter out everything that is not forest
ANO.forest <- ANO.geo %>% dplyr::filter(hovedoekosystem_rute == "Forest")
## add region information
nor <- st_read(here::here("data/outlineOfNorway_EPSG25833.shp"),
               quiet = T)%>%
  st_as_sf() %>%
  st_transform(crs = st_crs(ANO.forest))

reg <- st_read(here::here("data/regions.shp"),
               quiet = T) %>%
  st_as_sf() %>%
  st_transform(crs = st_crs(ANO.forest))

# change region names to something R-friendly
# reg$region
reg$region <- c("Northern.Norway","Central.Norway","Eastern.Norway","Western.Norway","Southern.Norway")

regnor <- st_intersection(reg,nor)

ANO.forest = st_join(ANO.forest, regnor, left = TRUE, join = st_nearest_feature)

### fix species names
ANO.sp$Species <- ANO.sp$art_navn
unique(as.factor(ANO.sp$Species))
ANO.sp[,'Species'] <- word(ANO.sp[,'Species'], 1,2) # lose subspecies
ANO.sp$Species <- str_to_title(ANO.sp$Species) # make first letter capital
ANO.sp$Species <- gsub("( .*)","\\L\\1",ANO.sp$Species,perl=TRUE) # make capital letters after hyphen to lowercase
ANO.sp$Species <- gsub("( .*)","\\L\\1",ANO.sp$Species,perl=TRUE) # make capital letters after space to lowercase

## merge species data with indicators
ANO.sp.ind <- merge(x=ANO.sp[,c("Species", "art_dekning", "parentglobalid")], 
                    y= ind.dat[,c("species", "Moisture", "Moisture")],
                    by.x="Species", by.y="species", all.x=T)
summary(ANO.sp.ind)


## checking which species didn't find a match
unique(ANO.sp.ind[is.na(ANO.sp.ind$Moisture),'Species'])

# fix species name issues
ind.dat <- ind.dat %>% 
  mutate(species=str_replace(species,"Aconitum lycoctonum", "Aconitum septentrionale")) %>% 
  mutate(species=str_replace(species,"Carex simpliciuscula", "Kobresia simpliciuscula")) %>%
  mutate(species=str_replace(species,"Carex myosuroides", "Kobresia myosuroides")) %>%
  mutate(species=str_replace(species,"Clinopodium acinos", "Acinos arvensis")) %>%
  mutate(species=str_replace(species,"Artemisia rupestris", "Artemisia norvegica")) %>%
  mutate(species=str_replace(species,"Cherleria biflora", "Minuartia biflora")) %>%
  mutate(species=str_replace(species,"Rosa vosagica", "Rosa vosagiaca"))

ANO.sp <- ANO.sp %>% 
  mutate(Species=str_replace(Species,"Agrostis hyemalis", "Agrostis scabra")) %>%
  mutate(Species=str_replace(Species,"Antennaria lapponica", "Antennaria alpina")) %>%
  mutate(Species=str_replace(Species,"Antennaria porsildii", "Antennaria alpina")) %>%
  mutate(Species=str_replace(Species,"Arctous alpinus", "Arctous alpina")) %>%
  mutate(Species=str_replace(Species,"Betula tortuosa", "Betula pubescens")) %>%
  mutate(Species=str_replace(Species,"Blysmopsis rufa", "Blysmus rufus")) %>%
  mutate(Species=str_replace(Species,"Cardamine nymanii", "Cardamine pratensis")) %>%
  mutate(Species=str_replace(Species,"Carex adelostoma", "Carex buxbaumii")) %>%
  mutate(Species=str_replace(Species,"Carex concolor", "Carex aquatilis")) %>%
  mutate(Species=str_replace(Species,"Carex leersii", "Carex echinata")) %>%
  mutate(Species=str_replace(Species,"Carex myosuroides", "Kobresia myosuroides")) %>%
  mutate(Species=str_replace(Species,"Carex paupercula", "Carex magellanica")) %>%
  mutate(Species=str_replace(Species,"Carex simpliciuscula", "Kobresia simpliciuscula")) %>%
  mutate(Species=str_replace(Species,"Carex viridula", "Carex flava")) %>%
  mutate(Species=str_replace(Species,"Chamaepericlymenum suecicum", "Cornus suecia")) %>%
  mutate(Species=str_replace(Species,"Cicerbita alpina", "Lactuca alpina")) %>%
  mutate(Species=str_replace(Species,"Cornus suecia", "Cornus suecica")) %>%
  mutate(Species=str_replace(Species,"Cotoneaster scandinavicus", "Cotoneaster integerrimus")) %>%
  mutate(Species=str_replace(Species,"Dactylorhiza viridis", "Coeloglossum viride")) %>%
  mutate(Species=str_replace(Species,"Diphasiastrum alpinum", "Lycopodium alpinum")) %>%
  mutate(Species=str_replace(Species,"Diphasiastrum complanatum", "Lycopodium complanatum")) %>%
  mutate(Species=str_replace(Species,"Dryopteris affinis", "Dryopteris filix-mas")) %>%
  mutate(Species=str_replace(Species,"Empetrum hermaphroditum", "Empetrum nigrum")) %>%
  mutate(Species=str_replace(Species,"Elymus alaskanus", "Elymus kronokensis")) %>%
  mutate(Species=str_replace(Species,"Festuca prolifera", "Festuca rubra")) %>%
  mutate(Species=str_replace(Species,"Galium album", "Galium mollugo")) %>%
  mutate(Species=str_replace(Species,"Galium elongatum", "Galium palustre")) %>%
  mutate(Species=str_replace(Species,"Helictotrichon pratense", "Avenula pratensis")) %>%
  mutate(Species=str_replace(Species,"Helictotrichon pubescens", "Avenula pubescens")) %>%
  mutate(Species=str_replace(Species,"Hieracium alpina", "Hieracium Alpina")) %>%
  mutate(Species=str_replace(Species,"Hieracium alpinum", "Hieracium Alpina")) %>%
  mutate(Species=str_replace(Species,"Hieracium hieracium", "Hieracium Hieracium")) %>%
  mutate(Species=str_replace(Species,"Hieracium hieracioides", "Hieracium umbellatum")) %>%
  mutate(Species=str_replace(Species,"Hieracium murorum", "Hieracium Vulgata")) %>%
  mutate(Species=str_replace(Species,"Hieracium oreadea", "Hieracium Oreadea")) %>%
  mutate(Species=str_replace(Species,"Hieracium prenanthoidea", "Hieracium Prenanthoidea")) %>%
  mutate(Species=str_replace(Species,"Hieracium vulgata", "Hieracium Vulgata")) %>%
  mutate(Species=str_replace(Species,"Hieracium pilosella", "Pilosella officinarum")) %>%
  mutate(Species=str_replace(Species,"Hieracium vulgatum", "Hieracium umbellatum")) %>%
  mutate(Species=str_replace(Species,"Hierochloã« alpina", "Hierochloë alpina")) %>%
  mutate(Species=str_replace(Species,"Hierochloã« hirta", "Hierochloë hirta")) %>%
  mutate(Species=str_replace(Species,"Hierochloã« odorata", "Hierochloë odorata")) %>%
  mutate(Species=str_replace(Species,"Huperzia appressa", "Huperzia selago")) %>%
  mutate(Species=str_replace(Species,"Huperzia arctica", "Huperzia selago")) %>%
  mutate(Species=str_replace(Species,"Hylotelephium maximum", "Sedum telephium")) %>%
  mutate(Species=str_replace(Species,"Listera cordata", "Neottia cordata")) %>%
  mutate(Species=str_replace(Species,"Leontodon autumnalis", "Scorzoneroides autumnalis")) %>%
  mutate(Species=str_replace(Species,"Loiseleuria procumbens", "Kalmia procumbens")) %>%
  mutate(Species=str_replace(Species,"Minuartia rubella", "Sabulina rubella")) %>%
  mutate(Species=str_replace(Species,"Minuartia stricta", "Sabulina stricta")) %>%
  mutate(Species=str_replace(Species,"Mycelis muralis", "Lactuca muralis")) %>%
  mutate(Species=str_replace(Species,"Omalotheca supina", "Gnaphalium supinum")) %>%
  mutate(Species=str_replace(Species,"Omalotheca norvegica", "Gnaphalium norvegicum")) %>%
  mutate(Species=str_replace(Species,"Omalotheca sylvatica", "Gnaphalium sylvaticum")) %>%
  mutate(Species=str_replace(Species,"Oreopteris limbosperma", "Thelypteris limbosperma")) %>%
  mutate(Species=str_replace(Species,"Oxycoccus microcarpus", "Vaccinium microcarpum")) %>%
  mutate(Species=str_replace(Species,"Oxycoccus palustris", "Vaccinium oxycoccos")) %>%
  mutate(Species=str_replace(Species,"Phalaris minor", "Phalaris arundinacea")) %>%
  mutate(Species=str_replace(Species,"Pinus unicinata", "Pinus mugo")) %>%
  mutate(Species=str_replace(Species,"Poa alpigena", "Poa pratensis")) %>%
  mutate(Species=str_replace(Species,"Poa angustifolia", "Poa pratensis")) %>%
  mutate(Species=str_replace(Species,"Poa ×jemtlandica", "Poa alpina")) %>%
  mutate(Species=str_replace(Species,"Potentilla anserina", "Argentina anserina")) %>%
  mutate(Species=str_replace(Species,"Potentilla arenosa", "Potentilla nivea")) %>%
  mutate(Species=str_replace(Species,"Pyrola grandiflora", "Pyrola rotundifolia")) %>%
  mutate(Species=str_replace(Species,"Rubus fruticosus", "Rubus plicatus")) %>%
  mutate(Species=str_replace(Species,"Rumex alpestris", "Rumex acetosa")) %>%
  mutate(Species=str_replace(Species,"Stellaria uliginosa", "Stellaria alsine")) %>%
  mutate(Species=str_replace(Species,"Syringa emodi", "Syringa vulgaris")) %>%
  mutate(Species=str_replace(Species,"Taraxacum crocea", "Taraxacum officinale")) %>%
  mutate(Species=str_replace(Species,"Taraxacum croceum", "Taraxacum officinale")) %>%
  mutate(Species=str_replace(Species,"Trientalis europaea", "Lysimachia europaea")) %>%
  mutate(Species=str_replace(Species,"Trifolium pallidum", "Trifolium pratense")) %>%
  mutate(Species=str_replace(Species,"Veratrum lobelianum", "Veratrum album"))

## merge species data with indicators
ANO.sp.ind <- merge(x=ANO.sp[,c("Species", "art_dekning", "parentglobalid")], 
                    y= ind.dat[,c("species", "Moisture", "Moisture")],
                    by.x="Species", by.y="species", all.x=T)
summary(ANO.sp.ind)
# checking which species didn't find a match
unique(ANO.sp.ind[is.na(ANO.sp.ind$Moisture),'Species'])
# don't find synonyms for these in the ind lists

## trimming away the points without information on NiN, species or cover
ANO.sp.ind <- ANO.sp.ind[!is.na(ANO.sp.ind$Species),]
ANO.sp.ind <- ANO.sp.ind[!is.na(ANO.sp.ind$art_dekning),]


summary(ANO.sp.ind)
#head(ANO.sp.ind)
rm(ANO.sp)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#### reference data - data handling

### generalized species lists for forest, mountain, wetland, and semi-natural ecosystems
str(Eco_State)

# species
Eco_State$Concept_Data$Species$Species_List$species
# environments
t(Eco_State$Concept_Data$Env$Env_Data)
# abundances
t(Eco_State$Concept_Data$Species$Species_Data)

## transposing abundance data
NiN.sp <- t(Eco_State$Concept_Data$Species$Species_Data)
NiN.sp <- as.data.frame(NiN.sp)
NiN.sp$sp <- as.factor(as.vector(Eco_State$Concept_Data$Species$Species_List$species))
# only genus and species name
NiN.sp$sp <- word(NiN.sp$sp, 1,2)
NiN.sp$spgr <- as.factor(as.vector(Eco_State$Concept_Data$Species$Species_List$art.code))
# if relevant, trimming to desired species groups (for forests e.g. removing trees)
#NiN.sp <- NiN.sp[NiN.sp$spgr!="a1a",]

## environment data
NiN.env <- Eco_State$Concept_Data$Env$Env_Data

## merging with indicator values
NiN.sp.ind <- merge(NiN.sp,ind.dat, by.x="sp", by.y="species", all.x=T)
summary(NiN.sp.ind)

NiN.sp.ind[NiN.sp.ind==999] <- NA

## checking which species didn't find a match
unique(NiN.sp.ind[is.na(NiN.sp.ind$Moisture) & NiN.sp.ind$spgr %in% list("a1a","a1b","a1c"),'sp'])

## fix species name issues
#ind.dat <- ind.dat %>% 
#  mutate(species=str_replace(species,"Aconitum lycoctonum", "Aconitum septentrionale")) %>% 
#  mutate(species=str_replace(species,"Carex simpliciuscula", "Kobresia simpliciuscula")) %>%
#  mutate(species=str_replace(species,"Carex myosuroides", "Kobresia myosuroides")) %>%
#  mutate(species=str_replace(species,"Clinopodium acinos", "Acinos arvensis")) %>%
#  mutate(species=str_replace(species,"Artemisia rupestris", "Artemisia norvegica")) %>%
#  mutate(species=str_replace(species,"Cherleria biflora", "Minuartia biflora"))

NiN.sp <- NiN.sp %>% 
  mutate(sp=str_replace(sp,"Aconitum lycoctonum", "Aconitum septentrionale")) %>% 
  mutate(sp=str_replace(sp,"Anagallis arvensis", "Lysimachia arvensis")) %>% 
  mutate(sp=str_replace(sp,"Anagallis minima", "Lysimachia minima")) %>% 
  mutate(sp=str_replace(sp,"Arctous alpinus", "Arctous alpina")) %>%
  mutate(sp=str_replace(sp,"Betula tortuosa", "Betula pubescens")) %>%
  mutate(sp=str_replace(sp,"Blysmopsis rufa", "Blysmus rufus")) %>%
  mutate(sp=str_replace(sp,"Chamerion angustifolium", "Chamaenerion angustifolium")) %>%
  mutate(sp=str_replace(sp,"Cardamine nymanii", "Cardamine pratensis")) %>%
  mutate(sp=str_replace(sp,"Carex adelostoma", "Carex buxbaumii")) %>%
  mutate(sp=str_replace(sp,"Carex leersii", "Carex echinata")) %>%
  mutate(sp=str_replace(sp,"Carex paupercula", "Carex magellanica")) %>%
  mutate(sp=str_replace(sp,"Carex simpliciuscula", "Kobresia simpliciuscula")) %>%
  mutate(sp=str_replace(sp,"Carex _vacillans", "Carex vacillans")) %>%
  mutate(sp=str_replace(sp,"Carex viridula", "Carex flava")) %>%
  mutate(sp=str_replace(sp,"Chamaepericlymenum suecicum", "Cornus suecia")) %>%
  mutate(sp=str_replace(sp,"Cornus suecia", "Cornus suecica")) %>%
  mutate(sp=str_replace(sp,"Cicerbita alpina", "Lactuca alpina")) %>%
  mutate(sp=str_replace(sp,"Dactylorhiza fuchsii", "Dactylorhiza maculata")) %>%
  mutate(sp=str_replace(sp,"Dactylorhiza sphagnicola", "Dactylorhiza majalis")) %>%
  mutate(sp=str_replace(sp,"Diphasiastrum alpinum", "Lycopodium alpinum")) %>%
  mutate(sp=str_replace(sp,"Diphasiastrum complanatum", "Lycopodium complanatum")) %>%
  mutate(sp=str_replace(sp,"Elymus alaskanus", "Elymus kronokensis")) %>%
  mutate(sp=str_replace(sp,"Empetrum hermaphroditum", "Empetrum nigrum")) %>%
  mutate(sp=str_replace(sp,"Erigeron acer", "Erigeron acris")) %>%
  mutate(sp=str_replace(sp,"Erigeron eriocephalus", "Erigeron uniflorus")) %>%
  mutate(sp=str_replace(sp,"Festuca altissima", "Drymochloa sylvatica")) %>%
  mutate(sp=str_replace(sp,"Festuca prolifera", "Festuca rubra")) %>%
  mutate(sp=str_replace(sp,"Galium album", "Galium mollugo")) %>%
  mutate(sp=str_replace(sp,"Galium elongatum", "Galium palustre")) %>%
  mutate(sp=str_replace(sp,"Glaux maritima", "Lysimachia maritima")) %>%
  mutate(sp=str_replace(sp,"Helictotrichon pratense", "Avenula pratensis")) %>%
  mutate(sp=str_replace(sp,"Helictotrichon pubescens", "Avenula pubescens")) %>%
  mutate(sp=str_replace(sp,"Hieracium alpina", "Hieracium Alpina")) %>%
  mutate(sp=str_replace(sp,"Hieracium alpinum", "Hieracium Alpina")) %>%
  mutate(sp=str_replace(sp,"Hieracium aurantiacum", "Pilosella aurantiaca")) %>%
  mutate(sp=str_replace(sp,"Hieracium dovrense", "Hieracium Alpestria")) %>%
  mutate(sp=str_replace(sp,"Hieracium hieracium", "Hieracium Hieracium")) %>%
  mutate(sp=str_replace(sp,"Hieracium hieracioides", "Hieracium umbellatum")) %>%
  mutate(sp=str_replace(sp,"Hieracium lactucella", "Pilosella lactucella")) %>%
  mutate(sp=str_replace(sp,"Hieracium murorum", "Hieracium Vulgata")) %>%
  mutate(sp=str_replace(sp,"Hieracium oreadea", "Hieracium Oreadea")) %>%
  mutate(sp=str_replace(sp,"Hieracium prenanthoidea", "Hieracium Prenanthoidea")) %>%
  mutate(sp=str_replace(sp,"Hieracium vulgata", "Hieracium Vulgata")) %>%
  mutate(sp=str_replace(sp,"Hieracium pilosella", "Pilosella officinarum")) %>%
  mutate(sp=str_replace(sp,"Hieracium vulgatum", "Hieracium umbellatum")) %>%
  mutate(sp=str_replace(sp,"Hierochloã« alpina", "Hierochloë alpina")) %>%
  mutate(sp=str_replace(sp,"Hierochloã« hirta", "Hierochloë hirta")) %>%
  mutate(sp=str_replace(sp,"Hierochloã« odorata", "Hierochloë odorata")) %>%
  mutate(sp=str_replace(sp,"Huperzia appressa", "Huperzia selago")) %>%
  mutate(sp=str_replace(sp,"Hylotelephium maximum", "Hylotelephium telephium")) %>%
  mutate(sp=str_replace(sp,"Lappula myosotis", "Lappula squarrosa")) %>%
  mutate(sp=str_replace(sp,"Lepidotheca suaveolens", "Matricaria discoidea")) %>%
  mutate(sp=str_replace(sp,"Listera cordata", "Neottia cordata")) %>%
  mutate(sp=str_replace(sp,"Listera ovata", "Neottia ovata")) %>%
  mutate(sp=str_replace(sp,"Leontodon autumnalis", "Scorzoneroides autumnalis")) %>%
  mutate(sp=str_replace(sp,"Loiseleuria procumbens", "Kalmia procumbens")) %>%
  mutate(sp=str_replace(sp,"Logfia arvensis", "Filago arvensis")) %>%
  mutate(sp=str_replace(sp,"Mentha _verticillata", "Mentha verticillata")) %>%
  mutate(sp=str_replace(sp,"Minuartia rubella", "Sabulina rubella")) %>%
  mutate(sp=str_replace(sp,"Minuartia stricta", "Sabulina stricta")) %>%
  mutate(sp=str_replace(sp,"Mycelis muralis", "Lactuca muralis")) %>%
  mutate(sp=str_replace(sp,"Omalotheca supina", "Gnaphalium supinum")) %>%
  mutate(sp=str_replace(sp,"Omalotheca norvegica", "Gnaphalium norvegicum")) %>%
  mutate(sp=str_replace(sp,"Omalotheca sylvatica", "Gnaphalium sylvaticum")) %>%
  mutate(sp=str_replace(sp,"Ononis arvensis", "Ononis spinosa")) %>%
  mutate(sp=str_replace(sp,"Oreopteris limbosperma", "Thelypteris limbosperma")) %>%
  mutate(sp=str_replace(sp,"Oxycoccus microcarpus", "Vaccinium microcarpum")) %>%
  mutate(sp=str_replace(sp,"Oxycoccus palustris", "Vaccinium oxycoccos")) %>%
  mutate(sp=str_replace(sp,"Phalaris minor", "Phalaris arundinacea")) %>%
  mutate(sp=str_replace(sp,"Phalaroides arundinacea", "Phalaris arundinacea")) %>%
  mutate(sp=str_replace(sp,"Pinus unicinata", "Pinus mugo")) %>%
  mutate(sp=str_replace(sp,"Platanthera montana", "Platanthera chlorantha")) %>%
  mutate(sp=str_replace(sp,"Poa alpigena", "Poa pratensis")) %>%
  mutate(sp=str_replace(sp,"Poa angustifolia", "Poa pratensis")) %>%
  mutate(sp=str_replace(sp,"Poa laxa", "Poa flexuosa")) %>%
  mutate(sp=str_replace(sp,"Poa _herjedalica", "Poa herjedalica")) %>%
  mutate(sp=str_replace(sp,"Poa _jemtlandica", "Poa jemtlandica")) %>%
  mutate(sp=str_replace(sp,"Poa jemtlandica", "Poa alpina")) %>%
  mutate(sp=str_replace(sp,"Poa lindebergii", "Poa arctica")) %>%
  mutate(sp=str_replace(sp,"Potentilla anserina", "Argentina anserina")) %>%
  mutate(sp=str_replace(sp,"Pyrola grandiflora", "Pyrola rotundifolia")) %>%
  mutate(sp=str_replace(sp,"Rhamnus catharticus", "Rhamnus cathartica")) %>%
  mutate(sp=str_replace(sp,"Rumex alpestris", "Rumex acetosa")) %>%
  mutate(sp=str_replace(sp,"Salix _fragilis", "Salix fragilis")) %>%
  mutate(sp=str_replace(sp,"Saxifraga _opdalensis", "Saxifraga opdalensis")) %>%
  mutate(sp=str_replace(sp,"Sorbus hybrida", "Hedlundia hybrida")) %>%
  mutate(sp=str_replace(sp,"Spergularia salina", "Spergularia marina")) %>%
  mutate(sp=str_replace(sp,"Syringa emodi", "Syringa vulgaris")) %>%
  mutate(sp=str_replace(sp,"Taraxacum crocea", "Taraxacum officinale")) %>%
  mutate(sp=str_replace(sp,"Taraxacum croceum", "Taraxacum officinale")) %>%
  mutate(sp=str_replace(sp,"Taraxacum erythrospermum", "Taraxacum officinale")) %>%
  mutate(sp=str_replace(sp,"Taraxacum hamatum", "Taraxacum officinale")) %>%
  mutate(sp=str_replace(sp,"Trientalis europaea", "Lysimachia europaea")) %>%
  mutate(sp=str_replace(sp,"Trifolium pallidum", "Trifolium pratense")) %>%
  mutate(sp=str_replace(sp,"Vicia orobus", "Vicia cassubica"))


## merge species data with indicators
NiN.sp.ind <- merge(NiN.sp,ind.dat, by.x="sp", by.y="species", all.x=T)
summary(NiN.sp.ind)

NiN.sp.ind[NiN.sp.ind==999] <- NA

# checking which species didn't find a match
unique(NiN.sp.ind[is.na(NiN.sp.ind$Moisture) & NiN.sp.ind$spgr %in% list("a1a","a1b","a1c"),'sp'])
# ok now

## matching with NiN ecosystem types - for wetlands, mountains, and forests
# NB! beware of rogue spaces in the 'Nature_type' & 'Sub_Type' variables, e.g. "Spring_Forest " & "Mountain "
NiN.forest <- NiN.sp.ind[,c("sp",paste(NiN.env[NiN.env$Nature_Type %in% c("Forest"),"ID"]),colnames(ind.dat)[c(16,18)])]   # Moisture, Moisture
NiN.forest[1,]
names(NiN.forest)

cbind(colnames(NiN.forest),
      c("",
        
        "T4-C1","T4-C5","T4-C9","T4-C13",
        rep("",4),
        "T4-C2","T4-C6","T4-C10","T4-C14",
        "T4-C3","T4-C7","T4-C11","T4-C15",
        "T4-C4","T4-C8","T4-C12","T4-C16",
        "T4-C17","T4-C18a","T4-C19a","T4-C18b",
        "T4-C19b","T4-C20",
        rep("",16),
        "T30-C1","T30-C2","T30-C3",
        rep("",3),
        "T30-C4",
        rep("",8),
        
        rep("",2) # indicators
      )
)

NiN.forest <- NiN.forest[,c(1,2:5,10:27,44:46,50,    # forest types
                            59:60       # indicators
)]

colnames(NiN.forest)[2:27] <- c(
  
  'T4-C1','T4-C5','T4-C9','T4-C13',
  'T4-C2','T4-C6','T4-C10','T4-C14',
  'T4-C3','T4-C7','T4-C11','T4-C15',
  'T4-C4','T4-C8','T4-C12','T4-C16',
  'T4-C17','T4-C18a','T4-C19a','T4-C18b',
  'T4-C19b','T4-C20',
  'T30-C1','T30-C2','T30-C3',
  'T30-C4'
  
)
#head(NiN.forest)


# translating the abundance classes into %-cover
coverscale <- data.frame(orig=0:6,
                         cov=c(0,1/32,1/8,3/8,0.6,4/5,1)
)

NiN.forest.cov <- NiN.forest
colnames(NiN.forest.cov)
for (i in 2:27) {
  NiN.forest.cov[,i] <- coverscale[,2][ match(NiN.forest[,i], 0:6 ) ]
}

NiN.forest.cov$sp <- as.factor(NiN.forest.cov$sp)

```

This leaves us with the monitoring data including plant indicators (ANO.sp.ind) and the reference data including plant indicators (NiN.forest.cov):
  
  <!-- Print head of data set with horizontal scrolling -->
  ```{r ano-sp-tab, echo=F}
head(ANO.sp.ind) %>%
  kable("html",caption = "ANO species occurrence data set with attached plant trait data.") %>% 
  kable_styling("striped") %>% scroll_box(width = "100%")
```

```{r wetland_mountain_forest_seminatural-ref-data-tab, echo=FALSE}
head(NiN.forest.cov) %>% 
  kable("html", caption = "Reference data set for wetland, mountain, forest, and semi-natural ecosystems") %>% 
  kable_styling("striped") %>% scroll_box(width = "100%")
```

<br />
  
  For each ecosystem type with a NiN species list, we can calculate a community weighted mean (CWM) for the relevant functional plant indicators.
In order to get distributions of these metrics rather than one single value (for comparison with the empirical testing data) the NiN lists can be bootstrapped.

<br />
  
  ##### Bootstrap functions for frequency abundance and cumulative density
  - function to calculate community weighted means (indBoot.freq() ) or cumulative density (indBoot.HeatOverhang() ) of selected indicator values (ind)
- for species lists (sp) with given abundances in percent (or on a scale from 0 to 1) in one or more 'sites' (abun)
- with a given number of iterations (iter),
- with species given a certain minimum abundance occurring in all bootstraps (obl), and
- with a given re-sampling ratio of the original species list (rat)
- in every bootstrap iteration the abundance of the sampled species can be randomly changed by a limited amount if wished by introducing a re-sampling of abundance values from adjacent abundance steps with a certain probability (var.abun)

```{r bootstrapping}
indBoot.freq <- function(sp,abun,ind,iter,obl,rat=2/3,var.abun=F) {
  
  ind.b <- matrix(nrow=iter,ncol=length(colnames(abun)))
  colnames(ind.b) <- colnames(abun)
  ind.b <- as.data.frame(ind.b)  
  
  ind <- as.data.frame(ind)
  ind.list <- as.list(1:length(colnames(ind)))
  names(ind.list) <- colnames(ind)
  
  for (k in 1:length(colnames(ind)) ) {
    ind.list[[k]] <- ind.b }
  
  for (j in 1:length(colnames(abun)) ) {
    
    dat <- cbind(sp,abun[,j],ind)
    dat <- dat[dat[,2]>0,]            # only species that are present in the ecosystem
    dat <- dat[!is.na(dat[,3]),]      # only species that have indicator values
    
    for (i in 1:iter) {
      
      speciesSample <- sample(dat$sp[dat[,2] < obl], size=round( (length(dat$sp)-length(dat$sp[dat[,2]>=obl])) *rat,0), replace=F)  
      dat.b <- rbind(dat[dat[,2] >= obl,],
                     dat[match(speciesSample,dat$sp),]
      )
      
      if (var.abun==T) {
        for (m in 1:nrow(coverscale[-1,]) ) {
          xxx <- dat.b[dat.b[,2]==coverscale[-1,][m,2],2]
          if ( m==1 ) { dat.b[dat.b[,2]==coverscale[-1,][m,2],2] <- sample( c(0.01,coverscale[2:7,2]), prob = c(0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0) ,size=length(xxx),replace=T) }
          if ( m==2 ) { dat.b[dat.b[,2]==coverscale[-1,][m,2],2] <- sample( c(0.01,coverscale[2:7,2]), prob = c(0.2, 0.3, 0.5, 0.0, 0.0, 0.0, 0.0) ,size=length(xxx),replace=T) }
          if ( m==3 ) { dat.b[dat.b[,2]==coverscale[-1,][m,2],2] <- sample( c(0.01,coverscale[2:7,2]), prob = c(0.0, 0.2, 0.3, 0.5, 0.0, 0.0, 0.0) ,size=length(xxx),replace=T) }
          if ( m==4 ) { dat.b[dat.b[,2]==coverscale[-1,][m,2],2] <- sample( c(0.01,coverscale[2:7,2]), prob = c(0.0, 0.0, 0.2, 0.3, 0.5, 0.0, 0.0) ,size=length(xxx),replace=T) }
          if ( m==5 ) { dat.b[dat.b[,2]==coverscale[-1,][m,2],2] <- sample( c(0.01,coverscale[2:7,2]), prob = c(0.0, 0.0, 0.0, 0.2, 0.3, 0.5, 0.0) ,size=length(xxx),replace=T) }
          if ( m==6 ) { dat.b[dat.b[,2]==coverscale[-1,][m,2],2] <- sample( c(0.01,coverscale[2:7,2]), prob = c(0.0, 0.0, 0.0, 0.0, 0.2, 0.3, 0.5) ,size=length(xxx),replace=T) }
        }
        dat.b[!is.na(dat.b[,2]) & dat.b[,2]<=(0),2] <- 0.01
        dat.b[!is.na(dat.b[,2]) & dat.b[,2]>1,2] <- 1
      }
      
      for (k in 1:length(colnames(ind))) {
        
        if ( nrow(dat.b)>2 ) {
          
          ind.b <- sum(dat.b[!is.na(dat.b[,2+k]),2] * dat.b[!is.na(dat.b[,2+k]),2+k] , na.rm=T) / sum(dat.b[!is.na(dat.b[,2+k]),2],na.rm=T)
          ind.list[[k]][i,j] <- ind.b
          
        } else {ind.list[[k]][i,j] <- NA}
        
      }
      
      #      print(paste(i,"",j)) 
    }
    
  }
  return(ind.list)
}
```

```{r}
#colnames(NiN.forest.cov)
# 1st column is the species
# 2nd-27th column is the abundances of sp in different ecosystem types
# 28th-29th column is the indicator values of the respective species
# we choose 1000 iterations
# species with abundance 0.8 or more (dominant species) must be included in each sample
# each sample re-samples 1/3 of the number of species
# the abundance of the re-sampled species may vary (see bootstrap function for details)
```

Running the bootstraps:
  ```{r boot2, eval = F}
forest.ref.cov <- indBoot.freq(sp=NiN.forest.cov[,1],abun=NiN.forest.cov[,2:27],ind=NiN.forest.cov[,28:29],
                               iter=1000,obl=0.8,rat=1/2,var.abun=T)

# fixing NaNs
for (i in 1:length(forest.ref.cov) ) {
  for (j in 1:ncol(forest.ref.cov[[i]]) ) {
    v <- forest.ref.cov[[i]][,j]
    v[is.nan(v)] <- NA
    forest.ref.cov[[i]][,j] <- v
  }
}

#saveRDS(forest.ref.cov, paste0(here::here(),"/data/cache/forest.ref.cov.RDS"))

```

```{r, include = F}
# Data from cache
forest.ref.cov<-readRDS(paste0(here::here(), "/data/cache/forest.ref.cov.RDS"))
```

```{r}
head(forest.ref.cov[[1]]) %>%
  kable("html", caption = "Table showing the first 6 rows of the bootstrapped data set.") %>% kable_styling("striped") %>% scroll_box(width = "100%")
```

This results in an R-list for forest ecosystem, with a data frame for every selected plant indicator and each data frame with as many columns as there are NiN species lists and as many rows as there were iterations in the bootstrap.

Next, we need to derive scaling values from these bootstrap-lists (the columns) for every mapping unit in NiN. Here, we define things in the following way:
  
  - Median = reference values
- 0.025 and 0.975 quantiles = lower and upper limit values
- min and max of the respective indicator's scale = min/max values



```{r scalingValues, attr.output='style="max-height: 300px;"', results='hide'}

# every NiN-type is represented by one 'generalisert artsliste'
# some NiN-types are represented by two such species lists
# in some cases two NiN-types are represented by the same species list
#head(forest.ref.cov[[1]])
forest.ref.cov[[1]][0,]

# NiN-types where each type is represented by one species list (including when one species list represents two NiN-types), i.e. excluding types with a- and b-suffix
names(forest.ref.cov[["Moisture"]])
x <- c(1:17,22:26)

# checking the actual NiN-types in the forest lists
forest.NiNtypes <- colnames(forest.ref.cov[["Moisture"]])
forest.NiNtypes[-x] <- substr(forest.NiNtypes[-x], 1, nchar(forest.NiNtypes[-x])-1)
forest.NiNtypes

# 2 indicator-value indicators: Moisture, "Moisture"
indEll.n=2
# creating a table to hold:
# for every indicator: the 0.5 quantile (median), 0.025 quantile and  0.975 quantile for each NiN-type (3 metrics -> 3 columns per indicator), for Heat_requirement in alpine/arctic systems we do 0.05 quantile and  0.95 quantile and ignore the former as it is a one-sided indicator
# for every nature type (nrows)
tab <- matrix(ncol=3*indEll.n, nrow=24 ) # 24 basic ecosystem types
# coercing the values into the table

myQuantiles <- c(0.025, 0.5, 0.975)

for (i in 1:length(x) ) {
  tab[i,1:3] <- quantile(as.matrix(forest.ref.cov[["Moisture"]][,x[i]]),probs=myQuantiles,na.rm=T)
  tab[i,4:6] <- quantile(as.matrix(forest.ref.cov[["Moisture"]][,x[i]]),probs=myQuantiles,na.rm=T)
}

tab <- as.data.frame(tab)
tab$NiN <- NA
tab$NiN[1:length(x)] <- names(forest.ref.cov[[1]])[x]
tab


# NiN-types represented by several species lists
forest.NiNtypes2 <- forest.NiNtypes[-x]
unique(forest.NiNtypes2)
grep(pattern=unique(forest.NiNtypes2)[1], x=forest.NiNtypes) # finds columns in e.g. colnames(forest.ref.cov[["Light"]]) that match the first NiN-type


for (i in 1:length(unique(forest.NiNtypes2)) ) {
  
  tab[length(x)+i,1:3] <- quantile(as.matrix(forest.ref.cov[["Moisture"]][,grep(pattern=unique(forest.NiNtypes2)[i], x=forest.NiNtypes)]),probs=myQuantiles,na.rm=T)
  tab[length(x)+i,4:6] <- quantile(as.matrix(forest.ref.cov[["Moisture"]][,grep(pattern=unique(forest.NiNtypes2)[i], x=forest.NiNtypes)]),probs=myQuantiles,na.rm=T)

  tab$NiN[length(x)+i] <- unique(forest.NiNtypes2)[i]
  
}

tab

# when species lists represent several NiN-types
# does not apply here
#tab$NiN
#tab <- rbind(tab,tab[c(27:30,42,44,47,55:57,60,62,71,83,85),])
#tab$NiN[c(27:30,42,44,47,55:57,60,62,71,83,85,88:102)] <- c("T3-C3","T3-C9","T3-C5","T3-C8",
#                                                            "T7-C13","T7-C7","T22-C2","T32-C1",
#                                                            "T32-C3","T32-C7","T45-C1","V10-C1","T22-C1",
#                                                            "T32-C5","T32-C21",
#                                                            
#                                                            "T3-C6","T3-C12","T3-C5","T3-C11",
#                                                            "T7-C14","T7-C9","T22-C4","T32-C2",
#                                                            "T32-C4","T32-C8","T45-C2","V10-C2","T22-C3",
#                                                            "T32-C20","T32-C6")
#tab

# making it a proper data frame
dim(tab)
round(tab[,1:6],digits=2)

colnames(tab) <- c("Moist_q2.5","Moist_q50","Moist_q97.5",
                   "Moisture_q2.5","Moisture_q50","Moisture_q97.5",

                   "NiN")
summary(tab)
tab$NiN <- gsub("C", "C-", tab$NiN) # add extra hyphen after C for NiN-types
tab


# restructuring into separate indicators for lower (q2.5) and higher (q97.5) than reference value (=median, q50)
y.Moist <- numeric(length=nrow(tab)*2)
y.Moist[((1:dim(tab)[1])*2)-1] <- tab$Moist_q2.5 
y.Moist[((1:dim(tab)[1])*2)] <- tab$Moist_q97.5 

y.Moisture <- numeric(length=nrow(tab)*2)
y.Moisture[((1:dim(tab)[1])*2)-1] <- tab$Moisture_q2.5 
y.Moisture[((1:dim(tab)[1])*2)] <- tab$Moisture_q97.5 

# creating final objects holding the reference and limit values for all indicators

# ref object for indicators
forest.ref.cov.val <- data.frame( grunn=c(rep(rep(tab$NiN,each=2),indEll.n)),
                                                    county=rep('all',(nrow(tab)*2*indEll.n)),
                                                    region=rep('all',(nrow(tab)*2*indEll.n)),
                                                    Ind=c(rep(c('Moist1','Moist2'),nrow(tab)),
                                                          rep(c('Moisture1','Moisture2'),nrow(tab))
                                                          
                                                    ),
                                                    Rv=c(rep(tab$Moist_q50,each=2),
                                                         rep(tab$Moisture_q50,each=2)
                                                    ),
                                                    Gv=c(y.Moist,y.Moisture),
                                                    maxmin=c(rep(c(1,12),nrow(tab)), # 12 levels of moisture
                                                             rep(c(1,9),nrow(tab))  # 9 levels of Moisture
                                                    )
)

forest.ref.cov.val
forest.ref.cov.val$grunn <- as.factor(forest.ref.cov.val$grunn)
forest.ref.cov.val$Ind <- as.factor(forest.ref.cov.val$Ind)
summary(forest.ref.cov.val)


```

```{r}
forest.ref.cov.val %>%
  kable("html", caption = "Reference values (Rv and threshold values (Gv) for each indicator and nature type combination") %>% kable_styling("striped") %>% scroll_box(width = "100%", height = "300px")
```

Once test data and the scaling values from the reference data are in place, we can calculate CWMs of the selected indicator for the ANO community data and scale them against the scaling values from the reference distribution. Note that we scale each ANO plot's CWM against either the lower threshold value and the min value OR the upper threshold value and the max value based on whether the CWM is smaller or higher than the reference value. Since the scaled values for both sides range between 0 and 1, we generate separate lower and upper indicators for each plant functional indicator type. An ANO plot can only have a scaled value in either the lower or the upper indicator (the other one will be 'NA'), except for the unlikely event that the CWM exactly matches the reference value, in which case both lower and upper indicator will receive a scaled indicator value of 1.

For scaling ANO-data against the reference we use the normalisation function from package ecTools.

```{r}
#| warning: false

#### calculating scaled and non-truncated values for the indicators based on the dataset ####
# drop geometry
ANO.forest <- ANO.forest %>%
  mutate(
    X = st_coordinates(.)[, 1],
    Y = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

ANO.forest <- ANO.forest |>
  mutate(
    kartleggingsenhet_1m2 = stringr::str_remove(kartleggingsenhet_1m2, "\\s.*$")
  )

for (i in 1:nrow(ANO.forest) ) {  
  tryCatch({
    #print(i)
    #print(paste(ANO.forest$ano_flate_id[i]))
    #print(paste(ANO.forest$ano_punkt_id[i]))
    #    ANO.forest$Hovedoekosystem_sirkel[i]
    #    ANO.forest$Hovedoekosystem_rute[i]
    
    
    
    # if the ANO.hovedtype exists in the forest reference
    if ( ANO.forest$hovedtype_rute[i] %in% unique(sub("\\-.*", "", forest.ref.cov.val$grunn)) ) {
      
      # if there is any species present in current ANO point  
      if ( length(ANO.sp.ind[ANO.sp.ind$parentglobalid==as.character(ANO.forest$globalid[i]),'Species']) > 0 ) {
        
        
        # Moisture
        dat <- ANO.sp.ind[ANO.sp.ind$parentglobalid==as.character(ANO.forest$globalid[i]),c('art_dekning','Moisture')]
        ANO.forest[i,'richness'] <- nrow(dat)
        dat <- dat[!is.na(dat$Moisture),]
        
        if ( nrow(dat)>0 ) {
          
          val <- sum(dat[,'art_dekning'] * dat[,'Moisture'],na.rm=T) / sum(dat[,'art_dekning'],na.rm=T)
          # lower part of distribution
          ref <- forest.ref.cov.val[forest.ref.cov.val$Ind=='Moisture1' & forest.ref.cov.val$grunn==as.character(ANO.forest[i,"kartleggingsenhet_1m2"]),'Rv']
          lim <- forest.ref.cov.val[forest.ref.cov.val$Ind=='Moisture1' & forest.ref.cov.val$grunn==as.character(ANO.forest[i,"kartleggingsenhet_1m2"]),'Gv']
          maxmin <- forest.ref.cov.val[forest.ref.cov.val$Ind=='Moisture1' & forest.ref.cov.val$grunn==as.character(ANO.forest[i,"kartleggingsenhet_1m2"]),'maxmin']
          # normalisation
          ANO.forest[i,'Moisture1'] <- ec_normalise(variable = val, x0 = maxmin, x60 = lim, x100 = ref, fun = "linear")
          
          # upper part of distribution
          ref <- forest.ref.cov.val[forest.ref.cov.val$Ind=='Moisture2' & forest.ref.cov.val$grunn==as.character(ANO.forest[i,"kartleggingsenhet_1m2"]),'Rv']
          lim <- forest.ref.cov.val[forest.ref.cov.val$Ind=='Moisture2' & forest.ref.cov.val$grunn==as.character(ANO.forest[i,"kartleggingsenhet_1m2"]),'Gv']
          maxmin <- forest.ref.cov.val[forest.ref.cov.val$Ind=='Moisture2' & forest.ref.cov.val$grunn==as.character(ANO.forest[i,"kartleggingsenhet_1m2"]),'maxmin']
          # normalisation
          ANO.forest[i,'Moisture2'] <- ec_normalise(variable = val, x0 = maxmin, x60 = lim, x100 = ref, fun = "linear")
          
        }
        
        
      }
      
    }
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# generate overall Moisture indicator
ANO.forest$Moisture <- ANO.forest$Moisture1*ANO.forest$Moisture2

ANO.forest$Moisture1[ANO.forest$Moisture1==1] <- NA
ANO.forest$Moisture2[ANO.forest$Moisture2==1] <- NA

# add geometry again
ANO.forest <- ANO.forest %>%
  st_as_sf(coords = c("X", "Y"), crs = st_crs(ANO.geo))

```

```{r}
head(ANO.forest) %>%
  kable("html", caption = "ANO plots for ecosystems with scaled indicator values in the far-right columns.") %>% 
  kable_styling("striped") %>% scroll_box(width = "100%")
```


We can visualise the normalisation with a generic figure (@fig-norm).

```{r fig-norm}
#| fig-cap: 'A figure illustrating the normalisation process.The placement of the x60 values are arbitrary in this figure and does not give a representative picture of the relationship between the variable and indicator values.'

var <- seq(0, 100, by = 0.5)
var2 <- seq(101, 200, by = 0.5)

ind <- ec_normalise(var, x0 = 0, x100 = 100, x60 = 50)
ind2 <- ec_normalise(var2, x100 = 101, x0 = 200, x60 = 150)

ec_norm_plot(
  variable = c(var, var2),
  indicator = c(ind, ind2)
) +
  scale_x_continuous(
    breaks = c(0, 50, 100, 150, 200),
    labels = c("lowest possible\nvalue", "2.5%", "50%", "97.5%", "highest possible\nvalue"),
    expand = expansion(add = 30)
  ) +
  labs(x = "Variable values as percentiles of \nbootstrapped reference values")
```

```{r}
#| eval: false
#| include: false

ec_norm_plot(
  variable = c(var, var2),
  indicator = c(ind, ind2)
) +
  scale_x_continuous(
    breaks = c(0, 50, 100, 150, 200),
    labels = c("lowest possible\nvalue", "2.5%", "50%", "97.5%", "highest possible\nvalue"),
    expand = expansion(add = 30)
  ) +
  labs(x = "Variabelverdier som persentiler ift. \nbootstrappede referanse verdier",
       y = "Indikatorverdi")

ggsave(here::here("img/norm_plot.png"))
```

#### Calculating national and regional indices

We can calculate national and region-wise indicator values using a mixed-effects, zero-one inflated beta-regression analysis. Since the indicator values are on a scale from 0-1 the data per definition don't follow a Gaussian distribution, hence the zero-one inflated beta-regression. Since the underlying monitoring data are spatially aggregated (several points per site), not all points are independent, hence the mixed-effects model.

The uncertainty around the national and region indicator values is derived from the same statistical model.

The national and regional estimates are computed separately.

```{r tmbs}
#| warning: false
## Fit ordered beta regression null-model with a random effect
# Norway

Moisture.no <- glmmTMB(
  Moisture ~ 1 + (1|ano_flate_id),
  data = ANO.forest,
  family = ordbeta
)

Moisture1.no <- glmmTMB(
  Moisture1 ~ 1 + (1|ano_flate_id),
  data = ANO.forest[!is.na(ANO.forest$Moisture1),],
  family = ordbeta
)

Moisture2.no <- glmmTMB(
  Moisture2 ~ 1 + (1|ano_flate_id),
  data = ANO.forest[!is.na(ANO.forest$Moisture2),],
  family = ordbeta
)

summary(Moisture.no)
summary(Moisture1.no)
summary(Moisture2.no)

# regions
Moisture.re <- glmmTMB(
  Moisture ~ 0 + region + (1|ano_flate_id),
  data = ANO.forest,
  family = ordbeta
)

Moisture1.re <- glmmTMB(
  Moisture1 ~ 0 + region + (1|ano_flate_id),
  data = ANO.forest[!is.na(ANO.forest$Moisture1),],
  family = ordbeta
)

Moisture2.re <- glmmTMB(
  Moisture2 ~ 0 + region + (1|ano_flate_id),
  data = ANO.forest[!is.na(ANO.forest$Moisture2),],
  family = ordbeta
)

summary(Moisture.re)
summary(Moisture1.re)
summary(Moisture2.re)
```


## 10. Results

The model estimates and standard errors are on a logit scale, but we can calculate indicator values and confidence intervals (CI) by transforming back from logit to values between 0-1
```{r}

NO_FUMO_004_table <- tibble(
  Region = c(
    "Norway",
    "Central Norway",
    "Eastern Norway",
    "Northern Norway",
    "Southern Norway",
    "Western Norway"
  ),
  indicator_value = c(
    plogis(summary(Moisture.no)$coef$cond[, 1]),
    plogis(summary(Moisture.re)$coef$cond[, 1])
  ),
  CI_upper = c(
    plogis(
      summary(Moisture.no)$coef$cond[, 1] +
        1.96 * summary(Moisture.no)$coef$cond[, 2]
    ),
    plogis(
      summary(Moisture.re)$coef$cond[, 1] +
        1.96 * summary(Moisture.re)$coef$cond[, 2]
    )
  ),
  CI_lower = as.numeric(
    plogis(
      summary(Moisture.no)$coef$cond[, 1] -
        1.96 * summary(Moisture.no)$coef$cond[, 2]
    ),
    plogis(
      summary(Moisture.re)$coef$cond[, 1] -
        1.96 * summary(Moisture.re)$coef$cond[, 2]
    )
  ),
  upper_indicator_value = c(
    plogis(summary(Moisture2.no)$coef$cond[, 1]),
    plogis(summary(Moisture2.re)$coef$cond[, 1])
  ),
  n_upper_indicator_value = c(
    summary(Moisture2.no)$nobs,
    as.vector(table(Moisture2.re$frame$region))
  ),
  lower_indicator_value = c(
    plogis(summary(Moisture1.no)$coef$cond[, 1]),
    plogis(summary(Moisture1.re)$coef$cond[, 1])
  ),
  n_lower_indicator_value = c(
    summary(Moisture1.no)$nobs,
    as.vector(table(Moisture1.re$frame$region))
  )
)
NO_FUMO_004_table
```

To express the results as a distribution, we need to simulate mean estimates based on the models' intercept and associated standard error.
```{r}

## Norway
# Fixed-effect intercept and SE
b0 <- fixef(Moisture.no)$cond["(Intercept)"]
se_b0 <- sqrt(vcov(Moisture.no)$cond["(Intercept)", "(Intercept)"])
# Dispersion parameter (For glmmTMB this is usually on the response/positive scale via sigma() )
phi <- sigma(Moisture.no)


# draw estimates and backtransform from logit-scale
eta_draws <- rnorm(1000, mean = b0, sd = se_b0)
mu_draws  <- plogis(eta_draws)


## regions
b_re <- fixef(Moisture.re)$cond
vc_re <- vcov(Moisture.re)$cond
phi <- sigma(Moisture.re)

eta_draws_reg <- MASS::mvrnorm(
  n = 1000,
  mu = b_re,
  Sigma = vc_re
)

mu_draws_reg <- plogis(eta_draws_reg)


# combining Moisture_no_sim and Moisture_re_sim
mu <- mu_draws_reg |>
  as_tibble() |>
  add_column(Norway = mu_draws) |>
  pivot_longer(everything(), names_to = "part") |>
  mutate(part = case_when(
    part == "regionCentral.Norway" ~ "C",
    part == "regionEastern.Norway" ~ "E",
    part == "regionNorthern.Norway" ~ "N",
    part == "regionWestern.Norway" ~ "W",
    part == "regionSouthern.Norway" ~ "S",
    .default = part
  )) |>
  add_column(year = "2024")
```


```{r tbl-results}
#| tbl-cap: 'Indicator values and 95% confidence intervals fro NO_FUMO_004, calculated for Norwegian alpine ecosystems.'
mu |>
  group_by(part) |>
  summarise(
    indicator = mean(value),
    low = round(quantile(value, probs = 0.025), 2),
    high = round(quantile(value, probs = .975), 2)
  ) |>
  unite(col = "95% CI", c(low, high), sep = " - ") |>
  gt::gt()
```

## 11. Export file

<!--# 
  
  Optional: Display the code (don't execute it) or the workflow for exporting the indicator values to file. Ideally the indicator values are exported as a georeferenced shape or raster file with indicators values, reference values and errors. You can also chose to export the raw (un-normalised or unscaled variable) as a separate product. You should not save large sptaial output data on GitHub. You can use eval=FALSE to avoid code from being executed (example below - delete if not relevant) 

-->

```{r export}
#| eval: false
arrow::write_parquet(mu,"data/results_NO_FUMO_004.parquet")
```

Demonstrating how to read back the data:

```{r readingBack}
arrow::read_parquet(here::here("data/results_NO_FUMO_004.parquet"))
```

## 12. Author contributions



**CRediT statement**

* Conceptualization: [JT]
* Methodology: [JT]
* Data curation: [JT, ALK]
* Formal analysis: [JT]
* Investigation: [JT]




**MeRIT attribution**

This part supplements the CRediT statement with more granularity for method section.

**Study design**
The study was designed by [JT].

**Data processing**
Raw data were curated and cleaned by [JT and AK].

**Statistical analysis**
All analyses were conducted by [JT] , with verification and minor edits by [AK].

**Visualization**
Figures were produced by [AK].

**Reproducibility**
All code was reviewed and validated by [AK].



::: {.callout-tip collapse="true"}
## Session Info
<!-- You can leave this last part as it is, unless you have a good reason to alter it. -->

This workflow uses renv. Run `renv::restore()` to get started.

:::