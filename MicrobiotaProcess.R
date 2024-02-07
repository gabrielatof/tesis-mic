#### Abundancias relativas con MicrobiotaProcess

transformar el objeto de phyloseq a MPSE
```{r}
p_0week #CS y MDD en semana 0
psd4_mdd #MDD semana 0 y semana 8

p_0_mpse <- p_0week %>% as.MPSE()

p_mdd_mpse <- psd4_mdd %>% as.MPSE()
```

Graficar top 20
### CS y MDD Week 0
```{r}
## a nivel familia 
p1 <- p_0_mpse %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=groups, 
    taxa.class = Family, 
    topn = 20,
    relative = TRUE
  )

## a nivel genero
p2 <- p_0_mpse %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=groups, 
    taxa.class = Genus, 
    topn = 20,
    relative = TRUE
  )

```

### MMD Week 0 y 08 
```{r}
## a nivel familia 
p3 <- p_mdd_mpse %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=groups, 
    taxa.class = Family, 
    topn = 20,
    relative = TRUE
  )

## a nivel genero
p4 <- p_mdd_mpse %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=groups, 
    taxa.class = Genus, 
    topn = 20,
    relative = TRUE
  )
```
