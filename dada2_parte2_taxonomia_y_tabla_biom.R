########################################################################################
################################## ASIGNACION DE TAXONOMIA ##################################

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
