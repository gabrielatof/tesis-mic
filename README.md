# ANÁLISIS DE PERFILES FUNCIONALES DE LA MICROBIOTA INTESTINAL DE PACIENTES CON TRASTORNO DEPRESIVO MAYOR

En este trabajo se incluyeron las secuencias del gen 16S del ARNr de los 63 pacientes hospitalizados y diagnosticados con Trastorno Depresivo Mayor (MDD), y los 30 controles saludables (CS) del trabajo publicado de [Dong et al. (2022)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9354493/). Los archivos SRA sin procesar se extrajeron de la NCBI [Short Read Archive](https://www.ncbi.nlm.nih.gov/sra), bajo el ID BioProject [PRJNA778934](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA778934). Se seleccionaron 122 archivos, excluyendo aquellos que no contenían información sobre el tipo de paciente y/o réplica de la muestra después del periodo de tratamiento.

El objetivo de este trabajo es predecir los perfiles funcionales de la microbiota intestinal asociado a pacientes con MDD a partir de secuencias del gen 16S del ARNr mencionadas anteriormente, así como identificar la diversidad bacteriana presente en pacientes con MDD y en los CS.

## CÓDIGOS

1. [sratoolkit.Rmd](https://github.com/gabrielatof/tesis/blob/main/sratoolkit.Rmd)

2. [tabla de metadatos](https://docs.google.com/spreadsheets/d/1EDimsD27WBPn68wyx7j0EZHouHNaX7XzLdx-4fQm128/edit?usp=sharing)

3. [dada2.R](https://github.com/gabrielatof/tesis-mic/blob/main/dada2.R)

4. [phyloseq.R](https://github.com/gabrielatof/tesis/blob/main/phyloseq.R)

5. [tabla BIOM](https://github.com/gabrielatof/tesis/blob/main/phyloseq.R)

6. [picrust2.py](https://github.com/gabrielatof/tesis/blob/main/picrust2.py)

Consultar la sección de [wiki](https://github.com/gabrielatof/tesis-mic/wiki/PICRUSt2) para el tutorial de instalación de PICRUSt2

7. [analisis_microbioma.R](https://github.com/gabrielatof/tesis/blob/main/analisis_microbioma.Rmd)

  7.1. [analisis_microbioma_resumido.R](url) 



## METODOLOGÍA

![Diagrama](https://github.com/gabrielatof/tesis-mic/blob/main/diagrama%20tesis.png)

## CONTACTO 

Este repositorio forma parte de la realización de una tesis de licenciatura por parte del **Laboratorio de Biología Cuantitativa y Sistemas Complejos** de la Licenciatura en Microbiología de la Universidad Autónoma de Querétaro. 

**Tesista:** [Gabriela Montserrat Torres Fernández](https://github.com/gabrielatof) 

**Correo de contacto:** gabrielatofdz@gmail.com 

