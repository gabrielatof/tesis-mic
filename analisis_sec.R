########################################################################################
################################## TUTORIAL DADA2 ##################################


##################################### CARGAR DATOS #####################################

library(dada2) # activar paquete a utilizar

path <- "~/tesis/01_Data/secuencias/sra" # direccionar en un objeto todas las secuencias 

head (list.files(path)) # verificar que se hayan cargado correctamente 

############################### ASIGNAR SECUENCIAS F y R ###############################

fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE)) #secuencia foward

fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE)) #secuencia reverse

###### Extraccion de los sample names #####

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) 

head (sample.names)

################################## ANALISIS DE CALIDAD ##################################

#Observar la calidad de las lecturas; Esto se hace para saber si necesitamos recortar o filtrar las secuencias (F y R)

plotQualityProfile(fnFs[1:6]) # empieza a caer en la posicion 280. Hay que quitar 20 nucleotidos


plotQualityProfile(fnRs[1:6]) #empieza a caer en la posicion 220, hay que quitar 80 nucleotidos

#Antes de cortar, primero hay que asignarles el nombre a cada muestra, el nombre filtrado
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names

names(filtRs) <- sample.names

#Se establecer los parametros para recortar la secuencia y el número máximo de errores permitidos en una lectura
#cambiamos con los valores de "truncLen" que son las posiciones de corte, lo demas se deja en los valores estandar 

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

 #### tasas de error ####

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

#Hay que visualizae estas tasas de error
plotErrors(errF, nominalQ=TRUE)

plotErrors(errR, nominalQ=TRUE)

####Inferencia de muestra#### 

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Guardamos el objeto, ya que tardo en cargar

save(dadaFs, file = "dadasFs1.RData")
save(dadaRs, file = "dadasRs1.RData")


################################## TABLA ASV ##################################

###### Combinar las lecturas #####

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

save(mergers, file = "merger.RData")

#### Construir la taba de secuencias ASV ####

seqtab <- makeSequenceTable(mergers)

dim (seqtab)

save (seqtab, file ="seqtabs.RData")

table (nchar(getSequences(seqtab))) #nspeccionar la distribucion del largo de las secuencias

############################ DADA TO FASTA #########################

#Generamos ya la tabla ASV pero esta no esta en el formato fasta para que lo podamos leer en PICRUSt2
#por esto, tenemos que buscar una funcion que lo haga, la siguiente funcion pasa de los formatos de 
#DADA2 a formato fasta, es una funcion que se encuentra en el paquete: "metagMisc" creada por vmikk.
#este paquete lo podemos encontrar en GitHub (https://github.com/vmikk/metagMisc)


dada_to_fasta <- function(seqtab, out = "DADA2.fasta", hash = "sha1", ...){
  
  seq_uniq <- dada2::getUniques(seqtab) 
  
  if(hash == "sha1"){ hh <- openssl::sha1(names(seq_uniq)) }
  if(hash == "sha256"){ hh <- openssl::sha256(names(seq_uniq)) }
  if(hash == "md5"){ hh <- openssl::md5(names(seq_uniq)) }
  
  seq_names <- paste(as.character(hh),
                     ";size=",
                     seq_uniq,
                     ";",
                     sep="")
  
  dada2::uniquesToFasta(seq_uniq, fout = out, ids = seq_names, ...)
  
  invisible(seq_names)
}

# La funcion ya quedo, y ya podemos utilizarla. El argumento es nuestra seqtab que creamos anteriormente 
# el hash es el predeterminado y para el formato de salida colocamos el ".fna"
dada_to_fasta(seqtab, out = "DADA2.fna", hash = "sha1")

##### eliminar quimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread=TRUE, verbose = TRUE)

save (seqtab.nochim, file="seqtabnochim.RData")

dim(seqtab.nochim) 

sum(seqtab.nochim)/sum(seqtab) 

##### verificar la calidad: Seguimiento de lecturas a través de la canalización ######

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

########################## ASIGNAR TAXONOMIA #############################

taxa <- assignTaxonomy(seqtab.nochim, "~/tesis/01_Data/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

#dada2 implemento un paso extra para poder asignar a nivel de especie, descargando una base extra

taxa <- addSpecies(taxa, "~/tesis/01_Data/silva_species_assignment_v138.1.fa.gz") #no ocrrio por Error: vector memory exhausted (limit reached?)

taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print) #llegamos a nivel de genero 

########################################################################
############################ TABLA BIOM ##################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomformat")

library (biomformat)

biomtable <- make_biom(taxa, sample_metadata = NULL, observation_metadata = NULL,
          id = NULL, matrix_element_type = "int")

save (biomtable, file ="table.biom")

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
