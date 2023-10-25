# Analisis de funcionalidad

```{r}
#paquetes a utilizar
library(microeco)
library(file2meco)
library(agricolae)
library(phyloseq)
library(microbiome)
```


## Semana 0
```{r}
# convertir de un objeto phyloseq a uno de microtable 
meco_dataset <- phyloseq2meco(p_0week)
meco_dataset
# Functional enrichment
t1 <- trans_func$new(meco_dataset)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
# use list to prepare data
tmp <- list()
# transpose res_spe_func_perc to be a data.frame like taxonomic abundance
tmp$func <- as.data.frame(t(t1$res_spe_func_perc), check.names = FALSE)
# assign the list as taxa_abund in your microtable object
meco_dataset$taxa_abund <- tmp
 # use trans_diff class to perform differential test
# Acá debes usar una comparación adecuada, por ejemplo: estado de salud, es decir una columna
# de tu sample_data(GlobalPatterns), en este caso usé SampleType
t2 <- trans_diff$new(dataset = meco_dataset, method = "anova", group = "groups", taxa_level = "all")
t2$plot_diff_abund(add_sig = T, color_values =c("#EE7600", "#008B8B")) + ggplot2::ylab("Relative abundance (%)")
ggsave("~/tesis/03_Output/funcion.png", dpi=600,width = 21,
  height = 13)
```

## pacientes con MDD semana 0 y 8
```{r}
# convertir de un objeto phyloseq a uno de microtable 
meco_dataset <- phyloseq2meco(p_mdd)
meco_dataset
# Functional enrichment
t1 <- trans_func$new(meco_dataset)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
# use list to prepare data
tmp <- list()
# transpose res_spe_func_perc to be a data.frame like taxonomic abundance
tmp$func <- as.data.frame(t(t1$res_spe_func_perc), check.names = FALSE)
# assign the list as taxa_abund in your microtable object
meco_dataset$taxa_abund <- tmp
 # use trans_diff class to perform differential test
t2 <- trans_diff$new(dataset = meco_dataset, method = "anova", group = "weeks", taxa_level = "all")
t2$plot_diff_abund(add_sig = T, color_values =c("#00868B", "#8B3A62")) + ggplot2::ylab("Relative abundance (%)")
ggsave("~/tesis/03_Output/funcion_dd.png", dpi=600,width = 21,
  height = 13)
```
