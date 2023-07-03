########################################################################
############################ PHYLOSEQ ##################################

#instalar phyloseq 
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

#activar paquetes a utilizar

library(phyloseq)
library(Biostrings)
library(ggplot2)
library(tidyverse)

#como el nombre de las secuencias no tienen 
#informacion sobre el tipo de muestra, 
#por separado se hara la lista de metadatos 
#para poder juntarlos con las secuencias a 


### metadatos para phyloseq

muestras<-c("DDZQ79_2","DDZQ79_1","DDZQ78_2","DDZQ78_1","DDZQ77_2",
            "DDZQ77_1","DDZQ76_2","DDZQ76_1","DDZQ75_2","DDZQ75_1",
            "DDZQ005_2","DDZQ74_2","DDZQ74_1","DDZQ73_2","DDZQ73_1",
            "DDZQ72_2","DDZQ72_1","DDZQ71_2","DDZQ71_1","DDZQ70_2",
            "DDZQ70_1","DDZQ005_1","DDZQ69_2","DDZQ69_1","DDZQ68_2",
            "DDZQ68_1","DDZQ67_2","DDZQ67_1","DDZQ64_1","DDZQ25_2",
            "N53","N52","N51","N37","DDZQ010_2","N36","N27","N26","N25",
            "N24","N23","N22","N21","N20","N19","DDZQ010_1","N18","N17",
            "N16","N14","N13","N12","N11","N10","N09","N08","DDZQ009_2",
            "N07","N06","DDZQ93_1","DDZQ91_1","DDZQ82_2","DDZQ82_1",
            "DDZQ80_2","DDZQ041_1","DDZQ038_2","DDZQ038_1","DDZQ037_2",
            "DDZQ037_1","DDZQ033_2","DDZQ033_1","DDZQ031_2","DDZQ031_1",
            "DDZQ030_2","DDZQ030_1","DDZQ026_2","DDZQ026_1","DDZQ022_2",
            "DDZQ022_1","DDZQ019_2","DDZQ019_1","DDZQ018_2","DDZQ018_1",
            "DDZQ016_2","DDZQ016_1","DDZQ012_2","DDZQ012_1","DDZQ011_2",
            "DDZQ011_1","N57","N56","N55","N54","DDZQ25_1","DDZQ064_2",
            "DDZQ004_2","DDZQ062_2","DDZQ062_1","DDZQ061_2","DDZQ061_1",
            "DDZQ059_2","DDZQ059_1","DDZQ056_2","DDZQ056_1","DDZQ055_2",
            "DDZQ004_1","DDZQ055_1","DDZQ054_2","DDZQ054_1","DDZQ053_2",
            "DDZQ053_1","DDZQ051_2","DDZQ051_1","DDZQ049_2","DDZQ049_1",
            "DDZQ048_2","DDZQ048_1","DDZQ046_2","DDZQ046_1","DDZQ041_2",
            "DDZQ80_1","DDZQ009_1")


grupos<-c("MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD",
          "MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD",
          "MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","CS","CS","CS","CS",
          "MDD","CS","CS","CS","CS","CS","CS","CS","CS","CS","CS","MDD","CS",
          "CS","CS","CS","CS","CS","CS","CS","CS","CS","MDD","CS","CS","MDD",
          "MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD",
          "MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD",
          "MDD","MDD","MDD","MDD","MDD","MDD","MDD","CS","CS","CS","CS","MDD",
          "MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD",
          "MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD","MDD",
          "MDD","MDD","MDD","MDD","MDD","MDD")

tiempo<-c("08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks",
          "0_Weeks","08_Weeks","0_Weeks","08_Weeks","08_Weeks","0_Weeks","08_Weeks",
          "0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks",
          "0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks",
          "0_Weeks","08_Weeks","0_Weeks","0_Weeks","0_Weeks","0_Weeks","08_Weeks",
          "0_Weeks","0_Weeks","0_Weeks","0_Weeks","0_Weeks","0_Weeks","0_Weeks",
          "0_Weeks","0_Weeks","0_Weeks","0_Weeks","0_Weeks","0_Weeks","0_Weeks",
          "0_Weeks","0_Weeks","0_Weeks","0_Weeks","0_Weeks","0_Weeks","0_Weeks",
          "08_Weeks","0_Weeks","0_Weeks","0_Weeks","0_Weeks","08_Weeks","0_Weeks",
          "08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks",
          "0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks",
          "08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks",
          "0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks","0_Weeks","0_Weeks",
          "0_Weeks","0_Weeks","0_Weeks","08_Weeks","08_Weeks","08_Weeks","0_Weeks",
          "08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks",
          "0_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks",
          "0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks","08_Weeks","0_Weeks",
          "08_Weeks","0_Weeks","0_Weeks")


## empezamos a trabajar con el objeto seqtab.nochim, 
#el cual fue el ultimo objeto creado en DADA2, 
# que es la tabla ASV sin las secuencias quimericas


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa)) #objeto de physoleq creado con el objeto de salida de DADA2

df_muestras<-sample_data(data.frame(samples=muestras,weeks=tiempo, groups=grupos,
                                    row.names=sample_names(ps),
                                    stringsAsFactors=FALSE)) #creando un objeto para los metadatos

ps_1 = merge_phyloseq(ps, df_muestras) #union de metadatos y objeto phyloseq


dna <- Biostrings::DNAStringSet(taxa_names(ps_1))
names(dna) <- taxa_names(ps_1)
ps <- merge_phyloseq(ps_1, dna) 
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps #de esta manera guardamos las secuencias de ADN de los ASV en una ranura refseq de este objeto y
#tambien se cambio el nombre los taxones a una cadena corta para que aparezcan asi en las tablas y figuras


########### Diversidad alfa ########### 

plot_richness(ps, x="weeks", measures=c("Shannon", "Simpson"), color="groups")


########### Ordenada (diversidad beta) ########### 

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu)) 
#transforma los datos en proporciones apropiadas para las distancias 
#Bray-Curtis (metodo para cuantificar la disimilitud composicional entre dos sitios diferentes)

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray") 
#NMDS es un enfoque que produce una ordenacion basada en una matriz de distancia

plot_ordination(ps.prop, ord.nmds.bray, color="groups")

# Observacion en grafico de barras

#Para genero
top222 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:222]
ps.top222 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top222 <- prune_taxa(top222, ps.top222)
plot_bar(ps.top222, x="weeks", fill="Genus") + facet_wrap(~groups, scales="free_x")

#para familia
plot_bar(ps.top222, x="weeks", fill="Family") + facet_wrap(~groups, scales="free_x")


                                     
########### TABLA BIOM ##########
library(biomformat)
otu <-t(as(otu_table(ps),"matrix"))
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"~/Desktop/otu_biom.biom")
