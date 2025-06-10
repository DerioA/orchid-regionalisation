

### send these code lines to the consol. This will upload (and install if necessary) needed packages

if(!any(names(installed.packages()[,1])=="maptools")) install.packages("maptools")
library(maptools)

##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Below are several examples on how to explore the database and extract specific component
##### ##### ##### ##### ##### ##### ##### ##### ##### 

### Show all family in RAINBIO database
sort(as.character(unique(RAINBIO[,"family"])))

### Show all genera of a given family
show_genera_of_family("Ebenaceae")  ## Provide a list of genera ordered by total number of records
show_genera_of_family("Ebenaceae", order="alphabetical") ## Provide a list of genera by alphabetical order

### Show all species of a given genus, with corresponding number of records
show_species_of_genus('Asystasia')   ## Provide a list of species ordered by total number of records
show_species_of_genus('Asystasia', order="alphabetical") ## Provide a list of species by alphabetical order

### Show all infra-specific levels of given species
show_infra_sp_of_species('Greenwayodendron','suaveolens')

### List all countries represented in the database
sort(unique(as.character(RAINBIO[,"country"])))

### Extract data for a given species
EXTRACT <- Extract_Data(RAINBIO=RAINBIO, GENUS="Diospyros", SP="iturensis")
nrow(EXTRACT)
head(EXTRACT)

### Extract data for multiple species within a genus
EXTRACT <- Extract_Data(RAINBIO=RAINBIO, GENUS="Lophira", SP=c("alata","lanceolata"))
nrow(EXTRACT)
head(EXTRACT)

### Extract data for one country
EXTRACT <- Extract_Data(RAINBIO=RAINBIO, COUNTRY=c("Togo"))
nrow(EXTRACT)
head(EXTRACT)

### Extract data for multiple countries
EXTRACT <- Extract_Data(RAINBIO=RAINBIO, COUNTRY=c("Equatorial Guinea","Benin"))
nrow(EXTRACT)
head(EXTRACT)

### Extract data for multiple species in one country
EXTRACT <- Extract_Data(RAINBIO=RAINBIO, COUNTRY="Gabon", GENUS="Diospyros", SP=c("fragrans", "iturensis", "gabunensis", "iturensis"))
nrow(EXTRACT)
head(EXTRACT)

### Extract data for multiple species in multiple countries
EXTRACT <- Extract_Data(RAINBIO=RAINBIO, GENUS="Diospyros", SP=c("fragrans", "iturensis"), 
                        COUNTRY=c("Gabon","Congo - Kinshasa","Central African Republic"))
nrow(EXTRACT)
head(EXTRACT)


### Explore and extract RAINBIO database for each child database individually
sort(unique(RAINBIO[,"institutionCode"])) ### List of all child dataset

EXTRACT <- Extract_childatabase(RAINBIO, child="BRLU")
head(EXTRACT) ## Insight of the extract. The last column contain the original ID of the child dataset

### Extract only records with missing/problematic georeferences
EXTRACT <- Extract_childatabase(RAINBIO, child="BRLU",only.georef.problems=T)
head(EXTRACT) ## Insight of the extract. The last column contain the original ID of the child dataset


### Export in csv format the Extract result
write.csv(EXTRACT, "extract.csv")



