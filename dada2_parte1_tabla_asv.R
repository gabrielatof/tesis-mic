########################################################################################
################################## OBTENCION DE ASV's ##################################


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

#### aqui termina tutorial de DADA2

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
