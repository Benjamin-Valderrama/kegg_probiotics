
<!-- README.md is generated from README.Rmd. Please edit that file -->

# KEGG probiotics

<!-- badges: start -->
<!-- badges: end -->

The goal of kegg probiotics is to leverage the data available in the
[KEGG database](https://www.genome.jp/kegg/) to make finding bacterial
genomes with potential to degrade carbohydrates of interest while
metabolising neuroactive compounds an easier task. To do so, here we
share the script `find_probiotic_candidates.R`, which was successfully
used to produce the *in-silico* analysis conducted in the published
article: [Paper name](link%20to%20paper).

## Overview

The scripts starts from a list of KEGG compound IDs of carbohydrates (or
internal carbohydrate subunits of a polymer), and performs the following
steps using the information stored in KEGG:

1.  Connects carbohydrate’s compounds IDs to IDs of reactions where it
    participats
2.  Connects the reaction IDs to Enzymes (EC codes). If EC codes are
    provided as an additional optional input for some carbohydrates,
    those EC numbers are integrated into the next step.
3.  Connects EC numbers to Kegg Orthologs (KOs)
4.  Connects KOs to all bacterial genomes where that KO was found.

## Input example

The scripts takes a tabular input data frame in the following form:

| carbohydrate_name                                   | carbohydrate_kegg_id | enzyme   |
|:----------------------------------------------------|:---------------------|:---------|
| 1-Kestose                                           | C03661               | NA       |
| Arabinoxylan‐oligosaccharides (AXOS) : Arabinoxylan | C01889               | 3.2.1.55 |
| Beta-glucan                                         | C00551               | NA       |
| Chitooligosaccharides (COS) : Chitositan            | C00734               | NA       |
| Fructooligosaccharide : Sucrose                     | C00089               | NA       |
| Fructooligosaccharide : Inuline                     | C08659               | 3.2.1.7  |

where:

- `carbohydrate_name` Identifies the target carbohydrate in an
  human-friendly way (not used in the script) \[optional colum\].
- `carbohydrate_kegg_id` ID of the carbohydrate (or a subunit within a
  polymer) used in the KEGG db for its identification.
- `carbohydrate_name` ECs of enzymes known to interact with the
  carbohydrate of interest \[optional column\].

## Output example

There are two main outputs produced by the script.

1)  A tabular output with the bacterial genomes showing enzymes capable
    of degrading the carbohydrates of interest

| carbohydrate_name                                   | carbohydrate_kegg_id | rx     | definition                                 | ec       | ko     | abbreviation | t_code | organism                                       | taxonomy                                                                                   |
|:----------------------------------------------------|:---------------------|:-------|:-------------------------------------------|:---------|:-------|:-------------|:-------|:-----------------------------------------------|:-------------------------------------------------------------------------------------------|
| Arabinoxylan‐oligosaccharides (AXOS) : Arabinoxylan | C01889               | NA     | NA                                         | 3.2.1.55 | K01209 | esa          | T00586 | Cronobacter sakazakii ATCC BAA-894             | prokaryotes_bacteria_gammaproteobacteria_enterobacteria_cronobacter_sakazakii_atcc_baa-894 |
| Arabinoxylan‐oligosaccharides (AXOS) : Arabinoxylan | C01889               | NA     | NA                                         | 3.2.1.55 | K01209 | csk          | T02119 | Cronobacter sakazakii ES15                     | prokaryotes_bacteria_gammaproteobacteria_enterobacteria_cronobacter_sakazakii_es15         |
| Chitooligosaccharides (COS) : Chitositan            | C00734               | R02333 | Chitin + n H2O \<=\> Chitosan + n Acetate  | 3.5.1.41 | K01452 | bcr          | T00799 | Bacillus cereus AH187                          | prokaryotes_bacteria_firmicutes_bacilli_bacillus_cereus_ah187                              |
| Chitooligosaccharides (COS) : Chitositan            | C00734               | R02333 | Chitin + n H2O \<=\> Chitosan + n Acetate  | 3.5.1.41 | K01452 | btk          | T00179 | Bacillus thuringiensis serovar konkukian 97-27 | prokaryotes_bacteria_firmicutes_bacilli_bacillus_thuringiensis_serovar_konkukian_97-27     |
| Fructooligosaccharide : Sucrose                     | C00089               | R00801 | Sucrose + H2O \<=\> D-Fructose + D-Glucose | 3.2.1.10 | K01182 | kpe          | T00772 | Klebsiella variicola 342                       | prokaryotes_bacteria_gammaproteobacteria_enterobacteria_klebsiella_variicola_342           |
| Fructooligosaccharide : Sucrose                     | C00089               | R00801 | Sucrose + H2O \<=\> D-Fructose + D-Glucose | 3.2.1.10 | K01182 | kpk          | T03368 | Klebsiella variicola KP5-1                     | prokaryotes_bacteria_gammaproteobacteria_enterobacteria_klebsiella_variicola_kp5-1         |

2)  A tabular output with the bacterial genomes showing enzymes involved
    in the metabolism of neuroactive compounds

| gbm                             | prokaryotes_bacteria_bacteria_incertae_sedis_bacteria_candidate_phyla_division_tm6_bacterium_gw2011_gwf2_28_16 | prokaryotes_bacteria_bacteria_incertae_sedis_chazhemtobacterium_chazhemtobacterium_aquaticus | prokaryotes_bacteria_firmicutes_clostridia_alkaliphilus_oremlandii | prokaryotes_bacteria_firmicutes_others_veillonella_rodentium | prokaryotes_bacteria_firmicutes_others_sporanaerobacter_sp_jn_28 | prokaryotes_bacteria_firmicutes_clostridia_aminipila_terrae |
|:--------------------------------|---------------------------------------------------------------------------------------------------------------:|---------------------------------------------------------------------------------------------:|-------------------------------------------------------------------:|-------------------------------------------------------------:|-----------------------------------------------------------------:|------------------------------------------------------------:|
| MGB004\|Kynurenine synthesis    |                                                                                                              0 |                                                                                            0 |                                                                  0 |                                                            0 |                                                                0 |                                                           0 |
| MGB005\|Tryptophan synthesis    |                                                                                                              0 |                                                                                            0 |                                                                  0 |                                                            1 |                                                                0 |                                                           0 |
| MGB006\|Glutamate synthesis I   |                                                                                                              0 |                                                                                            0 |                                                                  0 |                                                            1 |                                                                1 |                                                           1 |
| MGB007\|Glutamate synthesis II  |                                                                                                              0 |                                                                                            0 |                                                                  1 |                                                            1 |                                                                0 |                                                           1 |
| MGB008\|Acetylcholine synthesis |                                                                                                              0 |                                                                                            0 |                                                                  0 |                                                            0 |                                                                0 |                                                           0 |
| MGB009\|Histamine synthesis     |                                                                                                              0 |                                                                                            0 |                                                                  0 |                                                            0 |                                                                0 |                                                           0 |

## Comments on the usage of the software

While the script doesn’t need the user to provide a list of already
known enzymes that can interact with carbohydrates of interest, adding
them can increase the performance of the script and make the search of
suitable candidates more thorough. Although the curation of such a list
is project specific, there are some resources the authors of this
repository find useful:

- Supplementary table 1 from the article [Metatranscriptomic holobiont
  analysis of carbohydrate-active enzymes in the millipede *Telodeinopus
  aoutii* (Diplopoda,
  Spirostreptida)](https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2022.931986/full)
- Supplementary table 1 & 2 from the article [The microbiome’s fiber
  degradation profile and its relationship with the host
  diet](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-022-01461-6)
- Fiber_Enzyme_matrix.csv from the
  [Github](https://github.com/borenstein-lab/IFDP) repository associated
  with the above-mentioned article.

Finally, the script here provided `find_probiotic_candidates.R` uses the
R package [`omixerRpm`](https://github.com/omixer/omixer-rpmR) which is
an R interface to the software
[Omixer-RPM](https://github.com/raeslab/omixer-rpm), which is
distributed under an academic Non-commercial Software License Agreement
for Omixer-RPM. Details on the conditions of that licence can be found
in the software’s [GitHub
account](https://github.com/raeslab/omixer-rpm?tab=License-1-ov-file).
