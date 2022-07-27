Final Project Figures
================
Alice Tolley
31/03/2022

# Figures 1A and 2 setup

## Load packages

The following packages were loaded in order to create Figures 3A and 5:

``` r
library(karyoploteR)
library(GenomeInfoDb)
library(GenomeInfoDbData)
library(BSgenome)
library(BiocManager)
```

## Custom genome

A custom genome, created using the bostau9 (ARS-UCD1.2) assembly, was
required as a base graph to plot SNP density, gene density and gene
location.

First a file stating the length of each chromosome within the *Bos
taurus* genome was created:

    chr start end
    chr1 1 158534110
    chr2 1 136231102
    chr3 1 121005158
    chr4 1 120000601
    chr5 1 120089316
    chr6 1 117806340
    chr7 1 110682743
    chr8 1 113319770
    chr9 1 105454467
    chr10 1 103308737
    chr11 1 106982474
    chr12 1 87216183
    chr13 1 83472345
    chr14 1 82403003
    chr15 1 85007780
    chr16 1 81013979
    chr17 1 73167244
    chr18 1 65820629
    chr19 1 63449741
    chr20 1 71974595
    chr21 1 69862954
    chr22 1 60773035
    chr23 1 52498615
    chr24 1 62317253
    chr25 1 42350435
    chr26 1 51992305
    chr27 1 45612108
    chr28 1 45940150
    chr29 1 51098607
    chrX 1 139009144

Data within “mybostau9.txt” was converted to GRanges as required by
KaryoploteR.

``` r
custom.genome <- toGRanges("mybostau9.txt")
```

# Figure 1A: SNP & Gene density

## Create input data

Two input files were created highlighting the position of identified
SNPs and Genes along the genome.

File one: SNP position within each chromosome (Displaying first four
lines of the file)

    chr pos
    chr1    29018
    chr1    38959
    chr1    38960
    ...

File two: Gene position within each chromosome (Displaying first four
lines of the file)

    chr start   end
    chr1    339070  350389
    chr1    475398  475516
    chr1    477378  477504
    ...

Both files were loaded into R and assigned names:

``` r
SNPpositions <- read.table("SNPdata.txt", header = TRUE)
GENEpositions <- read.table("GENEdata.txt", header = TRUE)
```

Data was converted to GRanges, as required by KaryoploteR.

``` r
SNPs <- toGRanges(SNPpositions)
GENEs <- toGRanges(GENEpositions)
```

## Plot the data

``` r
# Define plot parameters
pp <- getDefaultPlotParams(plot.type = 4)
pp$data1inmargin <- 0
pp$bottommargin <- 20

# Build the karyoplot
kp <- plotKaryotype(genome=custom.genome, plot.type=4, labels.plotter = NULL, plot.params = pp, main="Gene vs SNP Density")
kpAddChromosomeNames(kp, srt=45)

# Plot the data
kp <- kpPlotDensity(kp, data=GENEs, window.size = 10e6, col = transparent("orange", amount = 0.5), border=darker("orange", amount = 150))
kp <- kpPlotDensity(kp, data=SNPs, window.size = 10e6, col = transparent("red", amount = 0.8), border=darker("red", amount = 150))
```

<img src="Figure 1A.png" width="1536" />

# Figure 2: Gene position

## Create input data

The chromosomal position of genes associated with each pathway was
stated in a text file. Example below: genes associated with glutamate
signalling:

    chr start end name
    chr15 1889884 2497441 GRIA4
    chr1 6112447 6577795 GRIK1
    chr9 47889245 48618507 GRIK2
    chr25 8496962 8687352 GRIN2A

Each file containing genes within the same pathway was loaded into R and
given an assigned name:

``` r
Glutamate <- read.table("Genepositions_glutamate.txt", header = TRUE)
Gprotein <- read.table("Genepositions_Gprotein.txt", header = TRUE)
Immunity <- read.table("Genepositions_Immunity.txt", header = TRUE)
Other <- read.table("Genepositions_other.txt", header = TRUE)
Proteins <- read.table("Genepositions_proteins.txt", header = TRUE)
RNA <- read.table("Genepositions_RNA.txt", header = TRUE)
```

Convert to GRanges as required by KaryoploteR:

``` r
Glutamatelabel <- toGRanges(Glutamate)
Gproteinlabel <- toGRanges(Gprotein)
Immunitylabel <- toGRanges(Immunity)
Otherlabel <- toGRanges(Other)
Proteinslabel <- toGRanges(Proteins)
RNAlabel <- toGRanges(RNA)
```

## Plot the data:

``` r
# Plot the data along custom genome
kp <- plotKaryotype(genome=custom.genome, main = "Gene position")

kpPlotMarkers(kp, data=Glutamatelabel, labels=Glutamatelabel$name, adjust.label.position = FALSE,
              text.orientation = "horizontal", r1=0.2, cex=0.8, line.color = "white", label.color = "orange")

kpPlotMarkers(kp, data=Gproteinlabel, labels=Gproteinlabel$name, adjust.label.position = FALSE,
              text.orientation = "horizontal", r1=0.2, cex=0.8, line.color = "white", label.color = "red")

kpPlotMarkers(kp, data=Immunitylabel, labels=Immunitylabel$name, adjust.label.position = FALSE,
              text.orientation = "horizontal", r1=0.2, cex=0.8, line.color = "white", label.color = "green")

kpPlotMarkers(kp, data=Otherlabel, labels=Otherlabel$name, adjust.label.position = FALSE,
              text.orientation = "horizontal", r1=0.2, cex=0.8, line.color = "white", label.color = "black")

kpPlotMarkers(kp, data=Proteinslabel, labels=Proteinslabel$name, adjust.label.position = FALSE,
              text.orientation = "horizontal", r1=0.2, cex=0.8, line.color = "white", label.color = "purple")

kpPlotMarkers(kp, data=RNAlabel, labels=RNAlabel$name, adjust.label.position = FALSE,
              text.orientation = "horizontal", r1=0.2, cex=0.8, line.color = "white", label.color = "blue")

# Plot gene length on the chromosome
kpPlotRegions(kp, data=Glutamatelabel, data.panel = "ideogram")
kpPlotRegions(kp, data=Gproteinlabel, data.panel = "ideogram")
kpPlotRegions(kp, data=Immunitylabel, data.panel = "ideogram")
kpPlotRegions(kp, data=Otherlabel, data.panel = "ideogram")
kpPlotRegions(kp, data=Proteinslabel, data.panel = "ideogram")
kpPlotRegions(kp, data=RNAlabel, data.panel = "ideogram")

# Add a figure legend
legend(x = "bottomright", fill = c("orange", "red", "green", "purple", "blue", "black"), 
       legend = c("Glutamate Signalling", "G-protein Signalling", "Immunity pathways",
                  "Uncharacterised proteins", "RNA function / DNA transcription", "Various pathways"))
```

<img src="Figure 2.png" width="1536" />

Final formatting to separate overlapping gene names was carried out in
PowerPoint.

# Figure 1B: Number of SNPs vs Chromosome length

## Load packages

``` r
library(ggplot2)
library(ggpmisc)
```

## Create input data

A file stating chromosome length and number of SNPs identified was
created, according to the output of SNPEff eff:

    chr length noSNPs
    chr1 158534110 573
    chr2 136231102 363
    chr3 121005158 427
    chr4 120000601 483
    chr5 120089316 456
    chr6 117806340 475
    chr7 110682743 353
    chr8 113319770 340
    chr9 105454467 427
    chr10 103308737 447
    chr11 106982474 372
    chr12 87216183 610
    chr13 83472345 761
    chr14 82403003 307
    chr15 85007780 732
    chr16 81013979 548
    chr17 73167244 259
    chr18 65820629 473
    chr19 63449741 114
    chr20 71974595 191
    chr21 69862954 403
    chr22 60773035 116
    chr23 52498615 233
    chr24 62317253 205
    chr25 42350435 149
    chr26 51992305 179
    chr27 45612108 158
    chr28 45940150 206
    chr29 51098607 207
    chrX 139009144 903

Data was loaded and assigned a name:

``` r
Chromscatter <- read.table("scatterplot_Chromosomes.txt", header = TRUE)
```

## Plot the data

``` r
# Calculate regression line intercept and slope
require(stats)
reg<-lm(noSNPs ~ length, data = Chromscatter)
reg
```

    ## 
    ## Call:
    ## lm(formula = noSNPs ~ length, data = Chromscatter)
    ## 
    ## Coefficients:
    ## (Intercept)       length  
    ##   4.165e+01    3.889e-06

``` r
# Plot the data
my.formula = y ~ x
sp <- ggplot(Chromscatter, aes(x=length, y=noSNPs)) +
  labs(title="Number of SNPs against chromosome length",
       x ="Chromosome length (bp)", y = "Number of SNPs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point() +
  scale_x_continuous(labels = scales::comma) +
  
# Add regression line
  geom_abline(intercept = 4.165e+01, slope = 3.889e-06) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point()
print(sp)
```

![](Final_Project_Figures_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
<img src="Figure 1B.png" width="672" />

Final formatting to enlarge numbers on axis for ease of viewing was
carried out in PowerPoint.
