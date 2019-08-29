##----------------------
## Load the data for regional trend analysis  
## CM: 29/08/2019
##
##----------------------
## From Ray's RAMCore.rdata and updated RAM

## load Ray's version
load("../../data/RAMCore.rdata") 

variables <- c("BvB", "UvU", "MCatch")

## replace some of the region names
Region[Region == "Atlantic Ocean"] <- "Atlantic Ocean Tunas"
Region[Region == "European Union"] <- "European Union non Mediterranean"
Region[Region == "Europe non EU"] <- "Norway, Iceland, Faroe Islands"
Region[Region == "Indian Ocean"] <- "Indian Ocean Tunas"
Region[Region == "Pacific Ocean"] <- "Pacific Ocean Tunas"
Region[Region == "Russia Japan"] <- "Northwest Pacific"

region.df <- data.frame(stockid = colnames(BvT), region = Region, stringsAsFactors = FALSE)

library(plyr)

## create dataframe objects from the matrices
for(i in 1:length(variables)){
    print(i)
    var <- variables[i]
    df <- adply(get(var), c(1, 2))
    names(df) <- c("year", "stockid", var)
    df$year <- as.numeric(as.character(df$year))
    df <- merge(df, region.df, all.x = TRUE)
    df <- df[with(df, order(region, stockid, year)), ]
    assign(paste0(var, ".df"), df)
    rm(list = c("var", "df"))
}
