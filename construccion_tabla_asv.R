##### CONSTRUCCION DE TABLA ASV #####

library(dada2) 

##### 1. Cargar los datos #####

path <- "~/tesis/01_Data/secuencias/sra" #cambiamos al directorio que contiene todas las secuencias para redireccionarlas a un objeto

#Renombrar los datos y asignar las secuencias F y R
#Las secuencias fastq foward y reverse tienen el formato _1.fastq y _2.fastq   

fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE)) #secuencia foward
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE)) #secuencia reverse

#Extraccion de los sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

##### 2. Analisis de calidad ####

#Observar la calidad de las lecturas; Esto se hace para saber si necesitamos recortar o filtrar las secuencias (F y R)

plotQualityProfile(fnFs[1:2]) 

plotQualityProfile(fnRs[1:12])

#### 3. Filtrar y recortar #### 

#Primero hay que asignarles el nombre a cada muestra, el nombre filtrado
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#y lo siguiente a realizar es establecer los parametros para recortar la secuencia
# lo unico que ccambiamos con los valore de 'truncLen' que son las posiciones de corte, lo demas 
#lo dejamos en los valores predeterminados  

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

#La funcion filterandtrim se encarga de remover los  primers 
#y recortar / filtrar las secuencias para la calidad


#### 4. Tasas de error #### 

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

#Esta funcion se encarga de generar un modelo de error de nuestros datos

#Hay que visualizae estas tasas de error
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Se muestran las tasas de error para cada transición posible (A→C, A→G, …). 
#Los puntos son las tasas de error observadas para cada puntaje de calidad de consenso. 
#La línea negra muestra las tasas de error estimadas después de la convergencia del algoritmo. 
#La línea roja muestra las tasas de error esperadas según la definición nominal de la puntuación Q. 

#### 5. Inferencia de muestra #### 

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#esta función es para inferir las variantes de secuencia de amplicón  en lecturas
#directas (F) e inversas (R) de forma independiente


#despues hay que ver que infirio el aloritmo para ambas secuencias
dadaFs[[1]]
dadaRs[[1]]

##### 6. Combinar las lecturas #####
#hay que combinar ambas secuencias, eliminado asi el ruido y tener las secuencias completos 

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

##### 7. Construir la taba de secuencias ASV ####
seqtab <- makeSequenceTable(mergers)

dim (seqtab)

#hay que inspeccionar la distribucion del largo de las secuencias
table (nchar(getSequences(seqtab)))

##### 8. Eliminar quimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread=TRUE, verbose = TRUE)

dim(seqtab.nochim) 

sum(seqtab.nochim)/sum(seqtab)

#### 9. Seguimiento de lecturas a traves de la canalizacion #####

# Para verificar la calidad y ver si se conservaron la mayoria de las lecturas al pasar por todos los pasos

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)