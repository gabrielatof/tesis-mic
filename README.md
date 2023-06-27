# ANÁLISIS DE PERFILES FUNCIONALES DE LA MICROBIOTA INTESTINAL DE PACIENTES CON TRASTORNO DEPRESIVO MAYOR

En este trabajo se seleccionaron se incluyeron las secuencias del gen 16S del ARNr de los 63 pacientes hospitalizados y diagnosticados con Trastorno Depresivo Mayor (MDD), y los 30 controles saludables (CS) del trabajo publicado de [Dong et al. (2022)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9354493/). Los archivos SRA sin procesar se extrajeron de la NCBI [Short Read Archive](https://www.ncbi.nlm.nih.gov/sra), bajo el ID BioProject [PRJNA778934](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA778934). Se seleccionaron 122 archivos, exluyendo aquellos que no contenían información sobre el tipo de paciente y/o réplica de la muestra después del periodo de tratamiento.

El objetivo de este trabajo es predecir los perfiles funcionales de la microbiota intestinal asociado a pacientes con MDD a partir de secuencias del gen 16S del ARNr mencionadas anteriormente, así como identificar la diversidad bacteriana presente en pacientes con MDD y en los CS.

## CÓDIGOS

**1. Descargar base de datos**

Los archivos se descargaron con ayuda de la herramienta  `SRA Toolkit`. Para este paso consultar el código [sratoolkit.Rmd](https://github.com/gabrielatof/tesis/blob/main/sratoolkit.Rmd)

2. **Procesamiento de datos con DADA2 y construcción de tabla ASV**

Los archivos FASTQ se procesaron en R studio con el paquete `DADA2` para obtener la tabla de ASV. Este paso se puede consultar en [dada2_parte1_tabla_asv.R](https://github.com/gabrielatof/tesis/blob/main/dada2_parte1_tabla_asv.R)

- **Tabla BIOM**

Después se proceso esta tabla y se le asigno la taxonomía y se generó la tabla BIOM. Consultar en [dada2_parte2_taxonomia_y_tabla_biom.R](https://github.com/gabrielatof/tesis/blob/main/dada2_parte2_taxonomia_y_tabla_biom.R)

3. **Análisis del microbioma de pacientes con MDD y CS**

Con el paquete `phyloseq` se construyó un objeto, el cual va a ser utilizado para todos los analisis del microbioma posteriores. Consultar en [phyloseq.R](https://github.com/gabrielatof/tesis/blob/main/phyloseq.R), para los análisis de diversidad alfa y beta así como análisis de abundancia consultar en [analisis_microbioma.R](url)


4. **Análisis estadístico**

5. **Análisis de perfiles funcionales con PICRUSt2**

Este análisis se realizo con el sofware `PICRUSt2` que se ejecuta en Python. Se necesitan tres archivos de entrada que se generaron en los pasos anteriores: tabla asv, tabla biom y tabla de metadatos. El código se puede consular en [picrust2.py](https://github.com/gabrielatof/tesis/blob/main/picrust2.py)


## CONTACTO 

Este repositorio forma parte de la realización de una tesis de licenciatura por parte del **Laboratorio de Biología Cuantitativa y Sistemas Complejos** de la Licenciatura en Microbiología de la Universidad Autónoma de Querétaro. 

**Tesista:** [Gabriela Montserrat Torres Fernández](https://github.com/gabrielatof) 

**Correo de contacto:** gabrielatofdz@gmail.com 
