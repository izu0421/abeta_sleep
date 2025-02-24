---
title: "MR analysis of KCNAB1 and KCNAB2 on AD risk using pQTL"
author: "Yizhou Yu"
date: 'updated: <i>2023-07-28</i></h4>'
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

"KCNAB1","KCNAB2"


```r
genes_interest = c("Q14722","Q13303")
```


```r
exp_full = read.csv("dt/Index_ROSMAP_pQTLs_robust005_exp_format_nonclump_rdy4MR_fixedDirection.csv")


exp = exp_full[exp_full$id.exposure %in% genes_interest,]

exp_format = exp
exp_format = subset(exp_format, pval.exposure < 0.05)
exp_format <- clump_data(exp_format, clump_r2 = 0.95,clump_kb = 100,pop = "EUR")
```

```
## API: public: http://gwas-api.mrcieu.ac.uk/
```

```
## Please look at vignettes for options on running this locally if you need to run many instances of this command.
```

```
## Clumping Q13303, 37 variants, using EUR population reference
```

```
## Removing 21 of 37 variants due to LD with other variants or absence from LD reference panel
```

```
## Clumping Q14722, 10 variants, using EUR population reference
```

```
## Removing 4 of 10 variants due to LD with other variants or absence from LD reference panel
```


```r
#Curate Corresponding Outcome Summary Data 
## Read in Formatted ADGWAS data 
outcome_filename = "dt/format_AD_GWAS.tsv.gz"
ADGWAS_out_format <- read.table(outcome_filename, sep ="\t", header= TRUE, stringsAsFactors=FALSE)
```

Harmonise the exposure and outcome datasets, with both Exposure & Outcome dataframes using the same reference alleles (rsids).


```r
harmonised_data_nonclump <- harmonise_data(
  exposure_dat = exp_format, 
  outcome_dat = ADGWAS_out_format,
  action = 1)
```

```
## Harmonising ROSMAP_pQTLs_DLPFC (Q14722) and AD (ADGWAS)
```

```
## Harmonising ROSMAP_pQTLs_DLPFC (Q13303) and AD (ADGWAS)
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
## Analysing 'Q13303' on 'ADGWAS'
```

```
## Analysing 'Q14722' on 'ADGWAS'
```

```r
mr_results<-generate_odds_ratios(mr_results)

write.table(mr_results, file = paste0("dt_out/",paste(genes_interest,collapse = "_"),"_1mrresults.tsv"), sep="\t", row.names=FALSE)

mr_results
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["id.exposure"],"name":[1],"type":["chr"],"align":["left"]},{"label":["id.outcome"],"name":[2],"type":["chr"],"align":["left"]},{"label":["outcome"],"name":[3],"type":["chr"],"align":["left"]},{"label":["exposure"],"name":[4],"type":["chr"],"align":["left"]},{"label":["method"],"name":[5],"type":["chr"],"align":["left"]},{"label":["nsnp"],"name":[6],"type":["int"],"align":["right"]},{"label":["b"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["se"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["pval"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["lo_ci"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["up_ci"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["or"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["or_lci95"],"name":[13],"type":["dbl"],"align":["right"]},{"label":["or_uci95"],"name":[14],"type":["dbl"],"align":["right"]}],"data":[{"1":"Q13303","2":"ADGWAS","3":"AD","4":"ROSMAP_pQTLs_DLPFC","5":"MR Egger","6":"16","7":"-2.6130423","8":"1.4327153","9":"8.959165e-02","10":"-5.4211643","11":"0.1950796","12":"0.07331117","13":"0.004421995","14":"1.2154078"},{"1":"Q13303","2":"ADGWAS","3":"AD","4":"ROSMAP_pQTLs_DLPFC","5":"Weighted median","6":"16","7":"-1.2550556","8":"0.3276931","9":"1.281574e-04","10":"-1.8973341","11":"-0.6127772","12":"0.28505999","13":"0.149967886","14":"0.5418440"},{"1":"Q13303","2":"ADGWAS","3":"AD","4":"ROSMAP_pQTLs_DLPFC","5":"Inverse variance weighted","6":"16","7":"-1.3467706","8":"0.2091576","9":"1.202457e-10","10":"-1.7567195","11":"-0.9368216","12":"0.26007881","13":"0.172610181","14":"0.3918714"},{"1":"Q14722","2":"ADGWAS","3":"AD","4":"ROSMAP_pQTLs_DLPFC","5":"MR Egger","6":"6","7":"0.2737981","8":"1.0803050","9":"8.124177e-01","10":"-1.8435997","11":"2.3911959","12":"1.31494930","13":"0.158246755","14":"10.9265537"},{"1":"Q14722","2":"ADGWAS","3":"AD","4":"ROSMAP_pQTLs_DLPFC","5":"Weighted median","6":"6","7":"-0.3216262","8":"0.3115937","9":"3.019796e-01","10":"-0.9323498","11":"0.2890973","12":"0.72496912","13":"0.393627687","14":"1.3352217"},{"1":"Q14722","2":"ADGWAS","3":"AD","4":"ROSMAP_pQTLs_DLPFC","5":"Inverse variance weighted","6":"6","7":"-0.2232269","8":"0.2298312","9":"3.314166e-01","10":"-0.6736960","11":"0.2272422","12":"0.79993335","13":"0.509820806","14":"1.2551339"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>



```r
plt = mr_scatter_plot(mr_results, harmonised_data_nonclump)
```


```r
plt$Q13303.ADGWAS + theme_classic()
```

![](MR_pQTL_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
ggsave("fig/mr_scatter_plot_Q13303.ADGWAS.pdf", width = 9, height = 4)
```


```r
plt$Q14722.ADGWAS + theme_classic()
```

![](MR_pQTL_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
ggsave("fig/mr_scatter_plot_Q14722.ADGWAS.pdf", width = 9, height = 4)
```

