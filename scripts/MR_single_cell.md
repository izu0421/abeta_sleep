---
title: "MR analysis of KCNAB1 and KCNAB2 on AD risk and sleep"
author: "Yizhou Yu"
date: 'updated: <i>2023-07-27</i></h4>'
output:
  html_document:
    df_print: paged
    keep_md: yes
  word_document: default
  pdf_document: default
---

Load library 

```r
library(tidyverse)
```

```
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.4.0      ✔ purrr   0.3.5 
## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
## ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
## ✔ readr   2.1.3      ✔ forcats 0.5.2 
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library(TwoSampleMR) # Allows access to GWAS catalogues for analyses 
```

```
## TwoSampleMR version 0.5.6 
## [>] New: Option to use non-European LD reference panels for clumping etc
## [>] Some studies temporarily quarantined to verify effect allele
## [>] See news(package='TwoSampleMR') and https://gwas.mrcieu.ac.uk for further details
## 
## 
## Warning:
## You are running an old version of the TwoSampleMR package.
## This version:   0.5.6
## Latest version: 0.5.7
## Please consider updating using devtools::install_github('MRCIEU/TwoSampleMR')
```

## Curate exposure data from SC eQTL


```r
genes_interest = c("KCNAB1","KCNAB2")
```


```r
sc_path = "dt/sc_brain_data_robust_005_fixed/"
files = list.files(sc_path)

byrois_sceqtl = data.frame()

for (file in files){
  byrois_sceqtl_tmp = read.csv(paste0(sc_path,file))
  byrois_sceqtl_tmp$tissue_name = str_split_fixed(file,"\\.",2)[,1]
  
  byrois_sceqtl_tmp$gene = str_split_fixed(byrois_sceqtl_tmp$geneid.exposure,"_",2)[,1]
  
  byrois_sceqtl_keep = byrois_sceqtl_tmp[byrois_sceqtl_tmp$gene %in% genes_interest,]
  
  byrois_sceqtl = rbind(byrois_sceqtl,byrois_sceqtl_keep)
}

byrois_sceqtl$gene = str_split_fixed(byrois_sceqtl$geneid.exposure,"_",2)[,1]

write.csv(byrois_sceqtl,"dt_out/sc_exposure_dt.csv",row.names = F)
```


```r
# Format Data
exp_format_nonclump <- byrois_sceqtl %>% 
  rename (effect_allele.exposure=effect_allele, other_allele.exposure=other_allele) %>%
  mutate(exposure = tissue_name, mr_keep.exposure = "TRUE", pval_origin.exposure = "reported", id.exposure = paste0(tissue_name,"_",geneid.exposure), eaf.exposure = NA) %>%
  select(SNP, beta.exposure, exposure, mr_keep.exposure, pval.exposure, pval_origin.exposure, id.exposure, eaf.exposure, se.exposure, effect_allele.exposure, other_allele.exposure)
```



```r
#Then Clumping (Not done here)
#cortex_exp_format_clump <- clump_data(cortex_exp_format, clump_r2 = 0.001) # Default is 0.001 but adjusting requires justification and leads to bias/double counting

#Curate Corresponding Outcome Summary Data 
## Read in Formatted ADGWAS data 
outcome_filename = "dt/format_AD_GWAS.tsv.gz"
ADGWAS_out_format <- read.table(outcome_filename, sep ="\t", header= TRUE, stringsAsFactors=FALSE)

#Then Select SNPs!
#ADGWAS_out_format_select_nonclump <- ADGWAS_out_format[ADGWAS_out_format$SNP %in% exp_format_nonclump$SNP,] #Selecting on non-clumped exposures
```

Harmonise the exposure and outcome datasets, with both Exposure & Outcome dataframes using the same reference alleles (rsids).


```r
exp_format = exp_format_nonclump
exp_format = subset(exp_format, pval.exposure < 0.01)
exp_format <- clump_data(exp_format, clump_r2 = 0.95,clump_kb = 100,pop = "EUR")
```

```
## API: public: http://gwas-api.mrcieu.ac.uk/
```

```
## Please look at vignettes for options on running this locally if you need to run many instances of this command.
```

```
## Clumping Excitatory_KCNAB2_ENSG00000069424, 34 variants, using EUR population reference
```

```
## Removing 24 of 34 variants due to LD with other variants or absence from LD reference panel
```

```
## Clumping Excitatory_KCNAB1_ENSG00000169282, 26 variants, using EUR population reference
```

```
## Removing 14 of 26 variants due to LD with other variants or absence from LD reference panel
```

```
## Clumping Inhibitory_KCNAB2_ENSG00000069424, 91 variants, using EUR population reference
```

```
## Removing 63 of 91 variants due to LD with other variants or absence from LD reference panel
```

```
## Clumping Inhibitory_KCNAB1_ENSG00000169282, 9 variants, using EUR population reference
```

```
## Removing 1 of 9 variants due to LD with other variants or absence from LD reference panel
```

```
## Clumping Oligodendrocytes_KCNAB1_ENSG00000169282, 312 variants, using EUR population reference
```

```
## Removing 234 of 312 variants due to LD with other variants or absence from LD reference panel
```

```
## Clumping OPCs_KCNAB1_ENSG00000169282, 15 variants, using EUR population reference
```

```
## Removing 10 of 15 variants due to LD with other variants or absence from LD reference panel
```

```r
harmonised_data_nonclump <- harmonise_data(
  exposure_dat = exp_format, 
  outcome_dat = ADGWAS_out_format,
  action = 1)
```

```
## Harmonising Oligodendrocytes (Oligodendrocytes_KCNAB1_ENSG00000169282) and AD (ADGWAS)
```

```
## Harmonising Inhibitory (Inhibitory_KCNAB2_ENSG00000069424) and AD (ADGWAS)
```

```
## Harmonising Excitatory (Excitatory_KCNAB2_ENSG00000069424) and AD (ADGWAS)
```

```
## Harmonising Excitatory (Excitatory_KCNAB1_ENSG00000169282) and AD (ADGWAS)
```

```
## Harmonising OPCs (OPCs_KCNAB1_ENSG00000169282) and AD (ADGWAS)
```

```
## Harmonising Inhibitory (Inhibitory_KCNAB1_ENSG00000169282) and AD (ADGWAS)
```

```r
# Power Prune?
###harmonised_data_nonclump<-power_prune(harmonised_data_nonclump, method=1)

#Save Harmonised data
write.table(harmonised_data_nonclump, file = paste0("dt_out/",paste(genes_interest,collapse = "_"),"_AD_harmoniseddata.tsv"), sep="\t", row.names=FALSE)
```

## Perform MR analyses, sensitivity analyses, and compile reports. 


```r
mr_results<-mr(harmonised_data_nonclump, method_list=c("mr_egger_regression", "mr_weighted_median", "mr_ivw", "mr_wald_ratio")) #TIME CONSUMING
```

```
## Analysing 'Excitatory_KCNAB1_ENSG00000169282' on 'ADGWAS'
```

```
## Analysing 'Excitatory_KCNAB2_ENSG00000069424' on 'ADGWAS'
```

```
## Analysing 'Inhibitory_KCNAB1_ENSG00000169282' on 'ADGWAS'
```

```
## Analysing 'Inhibitory_KCNAB2_ENSG00000069424' on 'ADGWAS'
```

```
## Analysing 'Oligodendrocytes_KCNAB1_ENSG00000169282' on 'ADGWAS'
```

```
## Analysing 'OPCs_KCNAB1_ENSG00000169282' on 'ADGWAS'
```

```r
mr_results<-generate_odds_ratios(mr_results)

write.table(mr_results, file = paste0("dt_out/",paste(genes_interest,collapse = "_"),"_1mrresults.tsv"), sep="\t", row.names=FALSE)

#mr_pleiotropy<-mr_pleiotropy_test(harmonised_data_nonclump) #SAVE CSV File
#write.table(mr_pleiotropy, file = paste0(outputs,tissue_name,"_2mrpleiotropy.tsv"), sep="\t", row.names=FALSE)
#mr_heterogeneity<-mr_heterogeneity(harmonised_data_nonclump, method_list=c("mr_egger_regression", "mr_ivw")) #SAVE CSV File
#write.table(mr_heterogeneity, file = paste0(outputs,tissue_name,"_3mrheterogeneity.tsv"), sep="\t", row.names=FALSE)
#res_single<-mr_singlesnp(harmonised_data_nonclump)
#res_loo<-mr_leaveoneout(harmonised_data_nonclump)
```


```r
subset(mr_results,exposure == "Inhibitory")
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["id.exposure"],"name":[1],"type":["chr"],"align":["left"]},{"label":["id.outcome"],"name":[2],"type":["chr"],"align":["left"]},{"label":["outcome"],"name":[3],"type":["chr"],"align":["left"]},{"label":["exposure"],"name":[4],"type":["chr"],"align":["left"]},{"label":["method"],"name":[5],"type":["chr"],"align":["left"]},{"label":["nsnp"],"name":[6],"type":["int"],"align":["right"]},{"label":["b"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["se"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["pval"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["lo_ci"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["up_ci"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["or"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["or_lci95"],"name":[13],"type":["dbl"],"align":["right"]},{"label":["or_uci95"],"name":[14],"type":["dbl"],"align":["right"]}],"data":[{"1":"Inhibitory_KCNAB1_ENSG00000169282","2":"ADGWAS","3":"AD","4":"Inhibitory","5":"MR Egger","6":"8","7":"-0.12185198","8":"0.11040358","9":"0.312011774","10":"-0.33824300","11":"0.094539039","12":"0.8852794","13":"0.7130220","14":"1.0991521","_rn_":"7"},{"1":"Inhibitory_KCNAB1_ENSG00000169282","2":"ADGWAS","3":"AD","4":"Inhibitory","5":"Weighted median","6":"8","7":"0.02112287","8":"0.04101203","9":"0.606524446","10":"-0.05926071","11":"0.101506455","12":"1.0213475","13":"0.9424610","14":"1.1068371","_rn_":"8"},{"1":"Inhibitory_KCNAB1_ENSG00000169282","2":"ADGWAS","3":"AD","4":"Inhibitory","5":"Inverse variance weighted","6":"8","7":"0.01561378","8":"0.03119652","9":"0.616724963","10":"-0.04553141","11":"0.076758961","12":"1.0157363","13":"0.9554896","14":"1.0797818","_rn_":"9"},{"1":"Inhibitory_KCNAB2_ENSG00000069424","2":"ADGWAS","3":"AD","4":"Inhibitory","5":"MR Egger","6":"27","7":"-0.22493022","8":"0.08400825","9":"0.012913187","10":"-0.38958640","11":"-0.060274044","12":"0.7985719","13":"0.6773370","14":"0.9415065","_rn_":"10"},{"1":"Inhibitory_KCNAB2_ENSG00000069424","2":"ADGWAS","3":"AD","4":"Inhibitory","5":"Weighted median","6":"27","7":"-0.03399659","8":"0.01870166","9":"0.069088798","10":"-0.07065184","11":"0.002658653","12":"0.9665748","13":"0.9317862","14":"1.0026622","_rn_":"11"},{"1":"Inhibitory_KCNAB2_ENSG00000069424","2":"ADGWAS","3":"AD","4":"Inhibitory","5":"Inverse variance weighted","6":"27","7":"-0.03563546","8":"0.01337538","9":"0.007715843","10":"-0.06185121","11":"-0.009419715","12":"0.9649920","13":"0.9400227","14":"0.9906245","_rn_":"12"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>



```r
plt = mr_scatter_plot(mr_results, harmonised_data_nonclump)
```


```r
plt$Inhibitory_KCNAB2_ENSG00000069424.ADGWAS + theme_classic()
```

![](MR_single_cell_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
ggsave("fig/mr_scatter_plot_Inhibitory_KCNAB2_ENSG00000069424.ADGWAS.pdf", width = 9, height = 4)
```



