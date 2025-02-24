---
title: "Metabolomic comparison of Abeta42 vs Abeta42Arc"
author: "Yizhou Yu"
date: "updated: <i>Dec-10-2022</i></h4>"
output:
  html_document:
    df_print: paged
    keep_md: yes
---

# Data processing and curation


```r
library(readxl)
library(ggplot2)
library(DEqMS)
```

```
## Loading required package: matrixStats
```

```
## Loading required package: limma
```

```r
library(ggrepel)
```

## Load data


```r
dt=read.csv("dt/3319961623_TMT_Proteins.txt", sep = "\t")
```


```r
gene_info_dt = dt[,c("Accession","Biological.Process",
                     "Molecular.Function","Gene.Symbol","WikiPathways")]
```



```r
colnames(dt)
```

```
##   [1] "Proteins.Unique.Sequence.ID"                    
##   [2] "Checked"                                        
##   [3] "Tags"                                           
##   [4] "Protein.FDR.Confidence.Combined"                
##   [5] "Master"                                         
##   [6] "Protein.Group.IDs"                              
##   [7] "Accession"                                      
##   [8] "Description"                                    
##   [9] "Sequence"                                       
##  [10] "FASTA.Title.Lines"                              
##  [11] "Exp.q.value.Combined"                           
##  [12] "Contaminant"                                    
##  [13] "Sum.PEP.Score"                                  
##  [14] "Number.of.Decoy.Protein.Combined"               
##  [15] "Coverage.in.Percent"                            
##  [16] "Number.of.Peptides"                             
##  [17] "Number.of.PSMs"                                 
##  [18] "Number.of.Protein.Unique.Peptides"              
##  [19] "Number.of.Unique.Peptides"                      
##  [20] "Number.of.AAs"                                  
##  [21] "MW.in.kDa"                                      
##  [22] "calc.pI"                                        
##  [23] "Score.Sequest.HT.Sequest.HT"                    
##  [24] "Coverage.in.Percent.by.Search.Engine.Sequest.HT"
##  [25] "Number.of.PSMs.by.Search.Engine.Sequest.HT"     
##  [26] "Number.of.Peptides.by.Search.Engine.Sequest.HT" 
##  [27] "Biological.Process"                             
##  [28] "Cellular.Component"                             
##  [29] "Molecular.Function"                             
##  [30] "Pfam.IDs"                                       
##  [31] "GO.Accessions"                                  
##  [32] "Entrez.Gene.ID"                                 
##  [33] "Gene.Symbol"                                    
##  [34] "Gene.ID"                                        
##  [35] "Reactome.Pathway.Accessions"                    
##  [36] "Reactome.Pathway.Accessions.All"                
##  [37] "Reactome.Pathways"                              
##  [38] "WikiPathway.Accessions"                         
##  [39] "WikiPathways"                                   
##  [40] "Ensembl.Gene.ID"                                
##  [41] "Number.of.Protein.Annotation.Groups"            
##  [42] "Number.of.Protein.Pathway.Groups"               
##  [43] "Abundance.Ratio.ARC..42"                        
##  [44] "Abundance.Ratio.PLUS..42"                       
##  [45] "Abundance.Ratio.PLUS..ARC"                      
##  [46] "Abundance.Ratio.log2.ARC..42"                   
##  [47] "Abundance.Ratio.log2.PLUS..42"                  
##  [48] "Abundance.Ratio.log2.PLUS..ARC"                 
##  [49] "Abundance.Ratio.P.Value.ARC..42"                
##  [50] "Abundance.Ratio.P.Value.PLUS..42"               
##  [51] "Abundance.Ratio.P.Value.PLUS..ARC"              
##  [52] "Abundance.Ratio.Adj.P.Value.ARC..42"            
##  [53] "Abundance.Ratio.Adj.P.Value.PLUS..42"           
##  [54] "Abundance.Ratio.Adj.P.Value.PLUS..ARC"          
##  [55] "Abundances.Grouped.42"                          
##  [56] "Abundances.Grouped.ARC"                         
##  [57] "Abundances.Grouped.PLUS"                        
##  [58] "Abundances.Grouped.CV.in.Percent.42"            
##  [59] "Abundances.Grouped.CV.in.Percent.ARC"           
##  [60] "Abundances.Grouped.CV.in.Percent.PLUS"          
##  [61] "Abundances.Grouped.Count.42"                    
##  [62] "Abundances.Grouped.Count.ARC"                   
##  [63] "Abundances.Grouped.Count.PLUS"                  
##  [64] "Abundance.F1.127N.Sample.42"                    
##  [65] "Abundance.F1.128C.Sample.42"                    
##  [66] "Abundance.F1.129N.Sample.42"                    
##  [67] "Abundance.F1.130C.Sample.42"                    
##  [68] "Abundance.F1.133C.Sample.42"                    
##  [69] "Abundance.F1.126.Sample.ARC"                    
##  [70] "Abundance.F1.128N.Sample.ARC"                   
##  [71] "Abundance.F1.129C.Sample.ARC"                   
##  [72] "Abundance.F1.132C.Sample.ARC"                   
##  [73] "Abundance.F1.133N.Sample.ARC"                   
##  [74] "Abundance.F1.127C.Sample.PLUS"                  
##  [75] "Abundance.F1.130N.Sample.PLUS"                  
##  [76] "Abundance.F1.131N.Sample.PLUS"                  
##  [77] "Abundance.F1.131C.Sample.PLUS"                  
##  [78] "Abundance.F1.132N.Sample.PLUS"                  
##  [79] "Abundances.Normalized.F1.127N.Sample.42"        
##  [80] "Abundances.Normalized.F1.128C.Sample.42"        
##  [81] "Abundances.Normalized.F1.129N.Sample.42"        
##  [82] "Abundances.Normalized.F1.130C.Sample.42"        
##  [83] "Abundances.Normalized.F1.133C.Sample.42"        
##  [84] "Abundances.Normalized.F1.126.Sample.ARC"        
##  [85] "Abundances.Normalized.F1.128N.Sample.ARC"       
##  [86] "Abundances.Normalized.F1.129C.Sample.ARC"       
##  [87] "Abundances.Normalized.F1.132C.Sample.ARC"       
##  [88] "Abundances.Normalized.F1.133N.Sample.ARC"       
##  [89] "Abundances.Normalized.F1.127C.Sample.PLUS"      
##  [90] "Abundances.Normalized.F1.130N.Sample.PLUS"      
##  [91] "Abundances.Normalized.F1.131N.Sample.PLUS"      
##  [92] "Abundances.Normalized.F1.131C.Sample.PLUS"      
##  [93] "Abundances.Normalized.F1.132N.Sample.PLUS"      
##  [94] "Abundances.Count.F1.127N.Sample.42"             
##  [95] "Abundances.Count.F1.128C.Sample.42"             
##  [96] "Abundances.Count.F1.129N.Sample.42"             
##  [97] "Abundances.Count.F1.130C.Sample.42"             
##  [98] "Abundances.Count.F1.133C.Sample.42"             
##  [99] "Abundances.Count.F1.126.Sample.ARC"             
## [100] "Abundances.Count.F1.128N.Sample.ARC"            
## [101] "Abundances.Count.F1.129C.Sample.ARC"            
## [102] "Abundances.Count.F1.132C.Sample.ARC"            
## [103] "Abundances.Count.F1.133N.Sample.ARC"            
## [104] "Abundances.Count.F1.127C.Sample.PLUS"           
## [105] "Abundances.Count.F1.130N.Sample.PLUS"           
## [106] "Abundances.Count.F1.131N.Sample.PLUS"           
## [107] "Abundances.Count.F1.131C.Sample.PLUS"           
## [108] "Abundances.Count.F1.132N.Sample.PLUS"           
## [109] "Found.in.Fraction.F11"                          
## [110] "Found.in.Fraction.F12"                          
## [111] "Found.in.Fraction.F13"                          
## [112] "Found.in.Fraction.F14"                          
## [113] "Found.in.Fraction.F15"                          
## [114] "Found.in.Fraction.F16"                          
## [115] "Found.in.Fraction.F17"                          
## [116] "Found.in.Fraction.F18"                          
## [117] "Found.in.Fraction.F19"                          
## [118] "Found.in.Fraction.F110"                         
## [119] "Found.in.Fraction.F111"                         
## [120] "Found.in.Fraction.F112"                         
## [121] "Found.in.Fraction.F113"                         
## [122] "Found.in.Sample.in.S2.F1.127N.Sample.42"        
## [123] "Found.in.Sample.in.S5.F1.128C.Sample.42"        
## [124] "Found.in.Sample.in.S6.F1.129N.Sample.42"        
## [125] "Found.in.Sample.in.S9.F1.130C.Sample.42"        
## [126] "Found.in.Sample.in.S15.F1.133C.Sample.42"       
## [127] "Found.in.Sample.in.S1.F1.126.Sample.ARC"        
## [128] "Found.in.Sample.in.S4.F1.128N.Sample.ARC"       
## [129] "Found.in.Sample.in.S7.F1.129C.Sample.ARC"       
## [130] "Found.in.Sample.in.S13.F1.132C.Sample.ARC"      
## [131] "Found.in.Sample.in.S14.F1.133N.Sample.ARC"      
## [132] "Found.in.Sample.in.S3.F1.127C.Sample.PLUS"      
## [133] "Found.in.Sample.in.S8.F1.130N.Sample.PLUS"      
## [134] "Found.in.Sample.in.S10.F1.131N.Sample.PLUS"     
## [135] "Found.in.Sample.in.S11.F1.131C.Sample.PLUS"     
## [136] "Found.in.Sample.in.S12.F1.132N.Sample.PLUS"     
## [137] "Found.in.Sample.Group.42"                       
## [138] "Found.in.Sample.Group.ARC"                      
## [139] "Found.in.Sample.Group.PLUS"                     
## [140] "Number.of.Protein.Groups"                       
## [141] "Modifications"
```


```r
dt_cleaning = dt[dt$Master == "IsMasterProtein" & 
                   dt$Protein.FDR.Confidence.Combined != "Low" ,]
dt_cleaning = dt_cleaning[,c(7,17,79:93)]
dt_cleaning = na.omit(dt_cleaning)
dt_cleaning = subset(dt_cleaning, Accession != "sp")
```


```r
print("Total proteins detected: ")
```

```
## [1] "Total proteins detected: "
```

```r
print(nrow(dt_cleaning))
```

```
## [1] 5717
```

## Differential expression analysis 


```r
dt_de = dt_cleaning[,c(3:ncol(dt_cleaning))]
row.names(dt_de)<-dt_cleaning$Accession
```

### Log norm

```r
dt_de.log = log(dt_de)

boxplot(dt_de.log)
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 1 is not drawn
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 2 is not drawn
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 3 is not drawn
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 4 is not drawn
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 5 is not drawn
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 6 is not drawn
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 9 is not drawn
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 10 is not drawn
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 11 is not drawn
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 12 is not drawn
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 13 is not drawn
```

```
## Warning in bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group
## == : Outlier (-Inf) in boxplot 14 is not drawn
```

![](Abeta42vsArc_proteomics_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

### Design table 


```r
cond = as.factor(c("Ab42","Ab42","Ab42","Ab42","Ab42",
                   "AbArc","AbArc","AbArc","AbArc","AbArc",
                   "control","control","control","control","control"))

# The function model.matrix is used to generate the design matrix
design = model.matrix(~0+cond) # 0 means no intercept for the linear model
colnames(design) = gsub("cond","",colnames(design))
```


```r
x <- c("Ab42-control","AbArc-control")
contrast =  makeContrasts(contrasts=x,levels=design)
```

### DEqMS analysis


```r
fit1 <- lmFit(dt_de.log, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)
fit3$count = dt_cleaning$Number.of.PSMs
fit4 = spectraCounteBayes(fit3)
```


```r
DEqMS.results_Ab42 = outputResult(fit4,coef_col = 1)
DEqMS.results_AbArc = outputResult(fit4,coef_col = 2)
```

## Abeta42 vs Arc 


```r
comparison_Ab42 = DEqMS.results_Ab42[,c("gene","logFC","sca.P.Value","sca.adj.pval")]
comparison_AbArc = DEqMS.results_AbArc[,c("gene","logFC","sca.P.Value","sca.adj.pval")]
colnames(comparison_Ab42)<-c("gene","Ab42logFC","Ab42.P.Value","Ab42.adj.pval")
colnames(comparison_AbArc)<-c("gene","AbArclogFC","AbArc.P.Value","AbArc.adj.pval")
```


```r
merged_comparison = merge(comparison_Ab42,comparison_AbArc,
                          by = "gene")
write.csv(merged_comparison,
          "dt_out/merged_comparison.csv", row.names = F)
```

### Plot 


```r
ggplot(merged_comparison, aes(x=Ab42logFC, y=AbArclogFC)) + 
  geom_hline(yintercept = c(-0.3,0.3),linetype="dashed") + 
  geom_vline(xintercept = c(-0.3,0.3),linetype="dashed") + 
  geom_point(aes(colour = Ab42logFC >0.3 | 
                   Ab42logFC < -0.3 |
                 AbArclogFC > 0.3 |
                 AbArclogFC < -0.3),
             size=2.5) + 
  scale_colour_manual(name = '', values = setNames(c('#ED2024','#939598'),c(T, F)))+
  geom_smooth(method = "lm")+
  theme_classic() +theme(axis.text.x=element_text(colour="black"),
          axis.text.y=element_text(colour="black"))
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

![](Abeta42vsArc_proteomics_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
ggsave("fig/proteomics_comparison.pdf", width = 10, height = 4)
```

```
## `geom_smooth()` using formula = 'y ~ x'
```



### Analysis of uncoordinated expression


```r
merged_comparison_addLab = merge(merged_comparison,
                                 gene_info_dt, by.x = "gene",
                                 by.y = "Accession")
```

#### Coordinated changes


```r
ggplot(merged_comparison, aes(x=Ab42logFC, y=AbArclogFC)) + 
  geom_hline(yintercept = c(-0.3,0.3),linetype="dashed") + 
  geom_vline(xintercept = c(-0.3,0.3),linetype="dashed") + 
  geom_point(aes(colour = Ab42logFC >0.3 | 
                   Ab42logFC < -0.3 |
                 AbArclogFC > 0.3 |
                 AbArclogFC < -0.3),
             size=2.5) + 
  scale_colour_manual(name = '', values = setNames(c('#ED2024','#939598'),c(T, F)))+
  geom_smooth(method = "lm")+
  geom_label_repel(data=merged_comparison_addLab[merged_comparison$Ab42logFC > 2 | merged_comparison$Ab42logFC < -2 |merged_comparison$AbArclogFC > 2 | merged_comparison$AbArclogFC < -2,],aes(label = Gene.Symbol))+
  theme_classic() +theme(axis.text.x=element_text(colour="black"),
          axis.text.y=element_text(colour="black"))
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

![](Abeta42vsArc_proteomics_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
ggsave("fig/proteomics_comparison_withLab.pdf", width = 10, height = 4)
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

#### Uncoordinated changes


```r
merged_comparison_analysis = subset(merged_comparison_addLab,
                                    Ab42.P.Value < 0.05 & 
                                      AbArc.P.Value<0.05)
merged_comparison_analysis$ArcUp_42down = ifelse(
  merged_comparison_analysis$Ab42logFC < 0 & 
    merged_comparison_analysis$AbArclogFC>0,
  1,0
)
merged_comparison_analysis$ArcDown_42up = ifelse(
  merged_comparison_analysis$AbArclogFC < 0 & 
    merged_comparison_analysis$Ab42logFC>0,
  1,0
)
```


```r
summary(merged_comparison_analysis[,c("ArcDown_42up",
                                      "ArcUp_42down")])
```

```
##   ArcDown_42up       ArcUp_42down     
##  Min.   :0.000000   Min.   :0.000000  
##  1st Qu.:0.000000   1st Qu.:0.000000  
##  Median :0.000000   Median :0.000000  
##  Mean   :0.009667   Mean   :0.006445  
##  3rd Qu.:0.000000   3rd Qu.:0.000000  
##  Max.   :1.000000   Max.   :1.000000
```


```r
merged_comparison_analysis$both = merged_comparison_analysis$ArcDown_42up + merged_comparison_analysis$ArcUp_42down
summary(as.factor(merged_comparison_analysis$both))
```

```
##   0   1 
## 916  15
```


```r
merged_comparison_diff = subset(merged_comparison_analysis,
                                        both ==1)
write.csv(merged_comparison_diff,
          "dt_out/merged_comparison_bothChanged.csv", row.names = F)
```


```r
ggplot(merged_comparison, aes(x=Ab42logFC, y=AbArclogFC)) + 
  geom_hline(yintercept = c(-0.3,0.3),linetype="dashed") + 
  geom_vline(xintercept = c(-0.3,0.3),linetype="dashed") + 
  geom_point(aes(colour = Ab42logFC >0.3 | 
                   Ab42logFC < -0.3 |
                 AbArclogFC > 0.3 |
                 AbArclogFC < -0.3),
             size=2.5) + 
  scale_colour_manual(name = '', values = setNames(c('#ED2024','#939598'),c(T, F)))+
  geom_smooth(method = "lm")+
  geom_label_repel(data=merged_comparison_diff[!duplicated(merged_comparison_diff$ Gene.Symbol),],aes(label = Gene.Symbol))+
  theme_classic() +theme(axis.text.x=element_text(colour="black"),
          axis.text.y=element_text(colour="black"))
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

```
## Warning: ggrepel: 15 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

![](Abeta42vsArc_proteomics_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
ggsave("fig/proteomics_comparison_withLab_subset.pdf", width = 10, height = 4)
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

```
## Warning: ggrepel: 15 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

**Up in Arc, down in 42**

*These are not FDR-significant*
CG15093, enable 3-hydroxyisobutyrate dehydrogenase activity. Predicted to be involved in valine catabolic process. Located in mitochondrion. Is expressed in adult head; adult heart; and organism. Orthologous to human HIBADH

3-hydroxy-2-methylpropanoate + NAD(+) <=> 2-methyl-3-oxopropanoate + H(+) + NADH

Mal-A6 : maltose alpha-glucosidase activity --> Hydrolysis of terminal, non-reducing (1->4)-linked alpha-D-glucose residues with release of alpha-D-glucose --> + glucose metabolism?


Hk - Hyperkinetic: !!! beta-subunit of Sh K[+] channels, NADPH binds to Hk NADP+-boud form of Hk regulates sleep --> an alcohol + NADP+ = an aldehyde or a ketone + NADPH + H

La: transcription termination by RNA polymerase III

CG6870 : enable heme binding activity, involved in sterol biosynthetic process

hgo - Homogentisate 1,2-dioxygenase activity. It is involved in the biological process described with: tyrosine catabolic process; L-phenylalanine catabolic process; oxidation-reduction process; tyrosine metabolic process: homogentisate + O2 <=> 4-maleylacetoacetate + H(+)

**Up in 42, down in Arc**

*These are FDR-significant*
Lsd-1 - Lipid storage droplet-1 (Lsd-1) is a protein associated with lipid droplets. It protects lipid droplets from lipase mediated remobilization and facilitates lipolysis by serving as an anchoring point for lipases such as Hsl. Lsd-1 is involved in lipid storage amount regulation and energy homeostasis in concert with Lsd-2
<br>

Tsp42Ee - Tetraspanin 42ee, isoform b; Tetraspanin 42Ee, ortholgous to TSPAN3


<br>

CG8665 - 10-formyltetrahydrofolate dehydrogenase; Hydroxymethyl-, formyl- and related transferase activity; methyltransferase activity; formyltetrahydrofolate dehydrogenase activity; oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor. It is involved in the biological process described with: oxidation-reduction process; *one-carbon metabolic process*; biosynthetic process; 10-formyltetrahydrofolate catabolic process [a.k.a. FBgn0032945, CG8665-RA, DmelCG8665]

Action: 
(6S)-10-formyltetrahydrofolate + H2O + NADP+ = (6S)-5,6,7,8-tetrahydrofolate + CO2 + H+ + NADPH

<br>
Sytbeta - Synaptotagmin beta, isoform D; Calcium-dependent phospholipid binding; transporter activity. It is involved in the biological process described with: neurotransmitter secretion; synaptic vesicle exocytosis

*These are not FDR-significant*

CG13071: Is expressed in adult head

Abcd3: enable ATP binding activity; *ATPase-coupled transmembrane transporter activity*; and long-chain fatty acid transporter activity. 
human: ABCD3 --> transport long chain fatty acid across membrane? 

CT19169: orthologous to SLC22A16, enable transmembrane transporter activity, organic zwitterion transporter protein family which transports carnitine

CG33120: long-chain-alcohol O-fatty-acyltransferase, involved in triglyceride biosynthetic process


Tig: extracellular matrix protein and integrin ligand that accumulates at the embryonic and larval muscle attachment site, involved_in axon guidance.
Tiggrin is a secreted glycoprotein that contains an RGD motif and is considered to be a ligand of the PS2 integrin (Fogerty et al., 1994;Bunch et al., 1998). Embryos homozygous for a loss of function allele of Tiggrin have a subtle Fas II phenotype reminiscent of integrin mutants.



Double check with the calculations by cat

[46] "Abundance.Ratio.log2.ARC..42"                   
 [47] "Abundance.Ratio.log2.PLUS..42"                  
 [48] "Abundance.Ratio.log2.PLUS..ARC"                 
 [49] "Abundance.Ratio.P.Value.ARC..42"                
 [50] "Abundance.Ratio.P.Value.PLUS..42"               
 [51] "Abundance.Ratio.P.Value.PLUS..ARC"              
 [52] "Abundance.Ratio.Adj.P.Value.ARC..42"            
 [53] "Abundance.Ratio.Adj.P.Value.PLUS..42"           
 [54] "Abundance.Ratio.Adj.P.Value.PLUS..ARC"         
 

```r
dt_check = dt[,c(7,33,46:54)]
dt_check = subset(na.omit(dt_check), Accession != "sp")

dt_check_sig = subset(dt_check,
                                    Abundance.Ratio.Adj.P.Value.PLUS..42 < 0.05 & Abundance.Ratio.Adj.P.Value.PLUS..ARC<0.05)

dt_check_sig$ArcUp_42down = ifelse(
  dt_check_sig$Abundance.Ratio.log2.PLUS..42 < 0 & 
    dt_check_sig$Abundance.Ratio.log2.PLUS..ARC>0,
  1,0
)
dt_check_sig$ArcDown_42up = ifelse(
  dt_check_sig$Abundance.Ratio.log2.PLUS..ARC < 0 & 
    dt_check_sig$Abundance.Ratio.log2.PLUS..42>0,
  1,0
)
```



```r
subset(dt_check_sig, ArcDown_42up == 1 | ArcUp_42down == 1)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Accession"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Gene.Symbol"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Abundance.Ratio.log2.ARC..42"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Abundance.Ratio.log2.PLUS..42"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Abundance.Ratio.log2.PLUS..ARC"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Abundance.Ratio.P.Value.ARC..42"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Abundance.Ratio.P.Value.PLUS..42"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["Abundance.Ratio.P.Value.PLUS..ARC"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["Abundance.Ratio.Adj.P.Value.ARC..42"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["Abundance.Ratio.Adj.P.Value.PLUS..42"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["Abundance.Ratio.Adj.P.Value.PLUS..ARC"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["ArcUp_42down"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["ArcDown_42up"],"name":[13],"type":["dbl"],"align":["right"]}],"data":[{"1":"Q9VUJ6","2":"Sytbeta","3":"-0.29","4":"-0.18","5":"0.12","6":"1.123014e-06","7":"0.0008441019","8":"0.0009841609","9":"0.0004245483","10":"0.01592592","11":"0.02246005","12":"1","13":"0","_rn_":"3398"},{"1":"Q9VUJ3","2":"Sytbeta","3":"-0.29","4":"-0.18","5":"0.12","6":"1.123014e-06","7":"0.0008441019","8":"0.0009841609","9":"0.0004245483","10":"0.01592592","11":"0.02246005","12":"1","13":"0","_rn_":"9897"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>



