
# Install and load packages --------------------------------------------------------

# CRAN packages used
.cran.pkgs <- c("readxl", "tidyverse", "RColorBrewer")

# install CRAN packages if not done yet
.to_install <- !.cran.pkgs %in% installed.packages()
sapply(.cran.pkgs[.to_install], install.packages)


# Bioconductor package used (install if not done yet)
if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install("KEGGREST")


# github package used (install if not done yet)
if (!require("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("omixer/omixer-rpmR")


# load all packages
.all.pkgs <- c("KEGGREST", "omixerRpm", .cran.pkgs)
.loaded.pkgs <- sapply(.all.pkgs, require, character.only = T)


# check if all packages were loaded
if(all(.loaded.pkgs)){
    cat("All packages are installed and loaded")
  } else {
  cat(paste(names(.loaded.pkgs[.loaded.pkgs == FALSE]), "was not loaded"))
}




# read input tables in -----------------------------------------------------

target_carbohydrates <- read_csv(file = "inputs/carbohydrates_and_degrading_enzymes.csv") %>% 
  filter(!is.na(carbohydrate_kegg_id))


ec_to_ko <- read_delim("internal_files/ec_to_ko_map.txt", delim = " ") %>% 
  rename_with(.cols = everything(), .fn = ~paste0(., "_code")) %>% 
  mutate(across(.cols = everything(), 
                .fns = ~gsub(x = ., pattern = "ec:|ko:", replacement = "")))


microorganism_table <- read_tsv(file = "internal_files/kegg_organisms.txt", col_names = c("t_code", "abbreviation", "organism", "taxonomy")) %>% 
  filter(!grepl(x = taxonomy, pattern = "Eukaryotes"))



# functions ---------------------------------------------------------------

rx_from_compound <- function(compound){
  
  # Rx where compound participates
  reaction <- as.character(keggGet(compound)[[1]]$REACTION)
  
  # Tidy the output
  reaction_list <- ifelse(is.null(reaction),
                     # if no reaction
                     yes = NA,
                     # change space for ";" when more than 1 reaction per compound
                     no = gsub(pattern = " ", replacement = ";", x = reaction))

  return(reaction_list)
  
}

definition_from_rx <- function(rx){
  
  rx_definition_list <- ifelse(is.na(rx),
                               yes = NA,
                               no = keggGet(rx)[[1]]$DEFINITION)
  
  return(rx_definition_list)
}

ec_from_rx <- function(rx){
  
  enzyme_list <- ifelse(is.na(rx),
                   yes = NA,
                   no = keggGet(rx)[[1]]$ENZYME %>% paste(collapse = ";"))
  
  return(enzyme_list)
  
}

ko_from_ec <- function(ec, ko_map){
  
  kos <- ifelse(is.na(ec),
                yes = NA,
                no = filter(.data = ko_map, ec_code %in% ec) %>% 
                            pull(ko_code) %>% 
                            paste(collapse = ";")) 
  
  clean_kos <- ifelse(kos == "", yes = NA, no = kos)
  
  return(clean_kos)
  
}

abbreviation_from_ko <- function(ko){

  abbreviation_list <- ifelse(is.na(ko),
                              yes = NA,
                              no = keggGet(ko)[[1]]$GENES %>% # Which organisms have the KO of interest
                                   sapply(str_split, pattern = ":") %>% # From the list of organisms, we will keep just the abbreviation (3 letters code)
                                   sapply('[', 1) %>%
                                   paste(collapse = ";") %>% 
                                   as.character()
                              )
  #abbreviation_list <- as.character(abbreviation_list)
  return(abbreviation_list)
}


genome_wise_omixerRpm <- function(x, 
                                  module.db, 
                                  minimum.coverage, ...){
  
  
  # for this use case, annotation will always be 1, as we are analysed modules in single genomes
  genomes_neuroactive_potentials <- data.frame()
  
  # run rpm for each column (genome) independently
  for(col in 2:ncol(x)){
    
    genome_GBM <- rpm(x = x[,c(1, col)], 
                       module.db = module.db,
                       minimum.coverage = minimum.coverage)
    
    
    # get the abundances of the modules
    genome_neuroactive_potential <- genome_GBM@abundance
    
    
    genome_modules_names <- genome_GBM@db@module.names %>%
      as.data.frame() %>%
      rownames_to_column("code") %>%
      mutate(full_module = paste0(code, "|", V2)) %>%
      filter(code %in% genome_GBM@annotation$Module) %>%
      pull(full_module)
    
    
    
    # rename the rownames with human-redable pathways
    rownames(genome_neuroactive_potential) <- genome_modules_names
      
    
     
    # merge the neuroactive potential of the current genome with the dataframe
    # storing the neuroactive potentials of the previously analysed genomes
    genomes_neuroactive_potentials <- merge(x = genome_neuroactive_potential,
                                            y = genomes_neuroactive_potentials,
                                            by = 0, all = TRUE) %>% 
                                      transform(row.names = Row.names,
                                                Row.names = NULL)
    
    
    # change NAs to 0
    genomes_neuroactive_potentials[is.na(genomes_neuroactive_potentials)] <- 0
    
  }
  
  return(genomes_neuroactive_potentials)
}



# find probiotic candidates based on target prebiotics --------------------

prebiotic_degrading_enzymes <- target_carbohydrates %>%

  # get reactions (number and definitions) where the prebiotics of interest participate
  mutate(rx = purrr::map_chr(.x = carbohydrate_kegg_id,
                             .f = rx_from_compound)) %>% 
  separate_longer_delim(cols = rx, delim = ";") %>% 
  mutate(definition = purrr::map_chr(.x = rx,
                                     .f = definition_from_rx)) %>%  
  
  # get the EC of the enzymes involved in those reactions
  mutate(ec = purrr::map_chr(.x = rx,
                             .f = ec_from_rx)) %>% 
  separate_longer_delim(cols = ec, delim = ";") %>%
  
  # merge ECs got from the enzymes and those manually annotated
  mutate(ec = coalesce(enzyme, ec)) %>% 
  
  # get the KO of the enzymes involved in those reactions
  mutate(ko =  purrr::map_chr(.x = ec,
                              .f = ~ko_from_ec(ec = .x, ko_map = ec_to_ko))) %>% 
  separate_longer_delim(cols = ko, delim = ";")



prebiotic_degrading_bacteria <- prebiotic_degrading_enzymes %>%

  # get the taxonomic info of the bacteria having any of the enzymes identified before
  mutate(abbreviation = purrr::map(.x = ko,
                                   .f = abbreviation_from_ko)) %>% 
  mutate(abbreviation = unlist(abbreviation)) %>% 
  separate_longer_delim(abbreviation, delim = ";") %>% 
  mutate(abbreviation = tolower(abbreviation)) %>% 
  inner_join(x = .,
             y = microorganism_table,
             by = "abbreviation") 


# fix tanomy annotation to have the names at all levels
prebiotic_degrading_bacteria <- prebiotic_degrading_bacteria %>% 
  mutate(taxonomy = paste(taxonomy, word(organism, start = 2, end = -1, sep = " "), sep = ";"))





# output taxas capable of degrading target prebiotics ---------------------

if(dir.exists("outputs")){
  write_tsv(x = prebiotic_degrading_bacteria,
            file = "outputs/prebiotic_degradation_bacteria.tsv")
} else {
  
  dir.create("outputs")
  write_tsv(x = prebiotic_degrading_bacteria,
            file = "outputs/prebiotic_degradation_bacteria.tsv")
}




# find probiotic candidates based on GBMs pontential ----------------------

# codes used in KEGG to id the bacteria of interest
bacteria_codes <- prebiotic_degrading_bacteria %>% 
  select(t_code) %>% 
  unique() %>% 
  pull()


bacterial_genomes <- bacteria_codes %>% 
  purrr::map(keggLink, target = "ko") %>% 
  purrr::map_df(function(x) data.frame(ko = x,
                                       abbreviation = names(x),
                                       presence = 1,
                                       row.names = NULL)) %>% 
  
  mutate(ko = gsub(x = ko, pattern = ".*:", replacement = ""),
         abbreviation = gsub(x = abbreviation, pattern = ":.*", replacement = ""))


omixer_input <- bacterial_genomes %>% 
  as_tibble() %>% 
  inner_join(x = .,
             y = unique(select(.data = prebiotic_degrading_bacteria, abbreviation, taxonomy)),
             by = "abbreviation") %>% 
  select(!abbreviation) %>% 
  unique() %>% 
  pivot_wider(id_cols = ko,
              names_from = "taxonomy",
              values_from = "presence",
              values_fill = 0) %>% 
  rename(entry = ko)

omixer_input <- janitor::clean_names(omixer_input)



# Determine neuroactive potential of probiotics ---------------------------

library(omixerRpm)

# load the database with the bacterial neuroactive pathways
db <- omixerRpm::loadDB(omixerRpm::listDB()[grepl(x = omixerRpm::listDB(), pattern = "GBMs")])


bacterial_neuroactive_potential <- omixer_input %>% 
  genome_wise_omixerRpm(x = ., module.db = db, minimum.coverage = 1)




# # infer GBMs from genomic content
# GBM <- rpm(x = omixer_input,
#            module.db = db, 
#            annotation = 1,
#            minimum.coverage = 1)
# 
# # get the abundances of the modules
# neuroactive_potential_bacteria <- GBM@abundance
# 
# modules_names <- GBM@db@module.names %>% 
#   as.data.frame() %>% 
#   rownames_to_column("code") %>% 
#   mutate(full_module = paste0(code, "|", V2)) %>% 
#   filter(code %in% GBM@annotation$Module) %>% 
#   pull(full_module)
# 
# # rename the rownames with human-redable pathways
# rownames(neuroactive_potential_bacteria) <- modules_names
# 
# # round the GBM values to the floor
# neuroactive_potential_bacteria <- neuroactive_potential_bacteria %>% floor()




# neuroactive potential of prebiotic-degrading taxa -----------------------

if(dir.exists("outputs")){
  write_tsv(x = bacterial_neuroactive_potential,
            file = "outputs/bacterial_neuroactive_potential.tsv")
} else {
  
  dir.create("outputs")
  write_tsv(x = bacterial_neuroactive_potential,
            file = "outputs/bacterial_neuroactive_potential.tsv")
}
