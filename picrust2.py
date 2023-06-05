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
