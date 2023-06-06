#PICRUSt2 
#Para analisis de perfiles funcionales de microbiota de 
#personas con Trastorno Depresivo Mayor (MDD) y Controles Saludables (CS)

#necesitamos 3 archivos: la DADA2.fna (los ASV , metadata.tsv (la tabla de los metadatos) y table.biom

#paso 1: crear una nueva carpeta para poner los objetos de salida de PICRUSt2
mkdir picrust2_out_pipeline 
cd picrust2_out_pipeline

#paso 2: colocar las lecturas del ASV en el arbol de referencia 
place_seqs.py -s ../DADA2.fna -o out.tre -p 1 \
              --intermediate intermediate/place_seqs

#paso 3: predicción del estado oculto de familia de genes 
# Este paso es para ver que tan similar es cada ASV a una secuencia de referencia existente
hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 1 -n # secuencia 16S de referencia más cercana
hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 1 # predice el numero de EC por asv 
# Este paso puede tardar dependiendo la capacidad de la computadora (30 horas en una macOS  Monterey 12.3.1 - 16 RAM)

#para visualizar las tablas: 
zless -S marker_predicted_and_nsti.tsv.gz 
#la primera columna es el nombre de los ASV, seguido de el numero previsto de copias de 16S por ASV
#seguido de el valor NSTI (valor del índice de taxón de secuencia más cercana) por ASV.

zless -S EC_predicted.tsv.gz 
#Esta tabla muestra el numero de copias previsto de todos los numeros de clasificacion de enzimas 
#para casa ASV.

#paso 4: predicción de los metageomas 
