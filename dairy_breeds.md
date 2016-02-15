#Introduction

This is a R Markdown file to accompany the mansucript: “Rumen microbial composition in Holstein and Jersey cows is different under same dietary condition and is not affected by sampling method” by Paz et al. 2016 in AEM. It was written in [R markdown](http://rmarkdown.rstudio.com) and converted to html using the R knitr package. This enables us to embed the results of our analyses directly into the text to allow for a reproducible data analysis pipeline. A [github repository is available](https://github.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds).

##BEFORE YOU RENDER:

The analysis from Paz et al., 2016 can be recreated by using the associated R Markdown file. This is setup to be ran on Mac OS X and was initially performed with 4 GB RAM. No root access is neeeded. Analyses should also work in a linux environment as well if the linux versions of USEARCH and [anaconda package manager](https://www.continuum.io/downloads) are used.

  1. Run the bash script to create a virtual enironment and download/install programs **LOCALLY** with the anaconda package manager. This will recreate the same enivronment used during the analyses of the data in the manuscript.
  2. Render the R Markdown file with knitR to recreate the workflow and outputs.

Due to licensing constraints, USEARCH could not be included in the setup. To obtain a download link, go to the USEARCH [download page](http://www.drive5.com/usearch/download.html) and select version USEARCH v8.1.1756 for Mac OSX. **A link (expires after 30 days) will be sent to the provided email. Use the link as an argument for shell script below.**

Simply download the bash script from the github repository and run it (provide the link to download your licensed USEARCH version as an argument for setup.sh):

  1. wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/setup.sh
  2. chmod 775 setup.sh
  3. ./setup.sh usearch_link
  
**Anaconda is downloaded and prompts you during installataion of the packages above. The prompts are as follows:**

  1. Press enter to view the license agreement
  2. Press enter to read the license and q to exit
  3. Accept the terms
  4. Prompts you where to install anaconda. Simply type anaconda to create a directory within the current directory. Should be: [/Users/user/anaconda2] >>> anaconda
  5. No to prepend anaconda to your path. Choosing yes does not impact the installation.
  6. Will be asked a few times if you wish to proceed with installing the packages, agree to it.
  7. After installation, enter ‘source anaconda/bin/activate projectEnv’ on the command line to activate the virtual enviornment with all dependencies.

To convert the R markdown to html use the command: **render(“dairy_breeds.Rmd”)**. To start a R session and run the workflow, use these commands from within the direcotry you initiated installation:

 1. source anaconda/bin/activate projectEnv
 2. R
 3. install.packages("rmarkdown", repos='http://cran.us.r-project.org')
 4. install.packages("knitr", repos='http://cran.us.r-project.org')
 5. source("http://bioconductor.org/biocLite.R")
 6. biocLite("Heatplus", ask=FALSE, suppressUpdates=TRUE)
 7. library(rmarkdown)
 8. library(knitr)
 9. render("dairy_breeds/dairy_breeds.Rmd")
 
The following packages and knitR settings were used to compile this Markdown:

```{r}
install.packages("ggplot2", repos='http://cran.us.r-project.org')
install.packages("matrixStats", repos='http://cran.us.r-project.org')
install.packages("plyr", repos='http://cran.us.r-project.org')
install.packages("tidyr", repos='http://cran.us.r-project.org')
install.packages("biom", repos='http://cran.us.r-project.org')
install.packages("gplots", repos='http://cran.us.r-project.org')
install.packages("RColorBrewer", repos='http://cran.us.r-project.org')
install.packages("vegan", repos='http://cran.us.r-project.org')
install.packages("mvtnorm", repos='http://cran.us.r-project.org')
install.packages("modeltools", repos='http://cran.us.r-project.org')
install.packages("coin", repos='http://cran.us.r-project.org')

sessionInfo()

perm = 1e4
```

```{r knitr.settings, eval=TRUE, echo=TRUE, cache=TRUE}
opts_chunk$set("tidy" = TRUE)
opts_chunk$set("echo" = TRUE)
opts_chunk$set("eval" = TRUE)
opts_chunk$set("warning" = FALSE)
opts_chunk$set("cache" = TRUE)
```

#Production measures of the cows

Download the file containing the dry matter intake (kg/d) and milk yield (kg/d) across cows.

```{r, engine='bash'}
wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/production_measures.txt
```

Generate bar chart of lactation measures for Holstein and Jersey breeds

```{r}
library(tidyr)
library(plyr)
library(ggplot2)

production_wide <- read.table("production_measures.txt", header = T, sep = "\t")
production_long <- gather(production_wide, Measure, kg_d, DMI:MY)

production_summary <- ddply(production_long, c("Breed", "Measure"), summarise, N = length(kg_d), mean = mean(kg_d), sd = sd(kg_d), se = sd / sqrt(N))

# Error bars represent standard error of the mean
production_plot <- ggplot(production_summary, aes(x=Measure, y=mean, fill=Breed)) + ylab("kg/d") + scale_fill_brewer("Paired") + geom_bar(position=position_dodge(), stat="identity", color = "black", size = 0.5, width = 0.9) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), size = 0.5, width = 0.3, position=position_dodge(.9)) + theme(legend.title = element_blank(), axis.text = element_text(color = "black", size = 12, face = "bold"), axis.title = element_text(color = "black", size = 12, face = "bold"), legend.text = element_text(color = "black", size = 10, face = "bold")) + scale_x_discrete(breaks=c("DMI", "MY"), labels=c("Dry Matter Intake", "Millk Yield"))  

png("FigS1.png", units = "in", height = 5, width = 6, res = 300)
production_plot
dev.off()

pdf("FigS1.pdf", height = 6, width = 6)
production_plot
dev.off()
```

```{r}
library(car)

cow_measures <- read.table ("production_measures.txt", header = T, sep="\t")

#Levene test
#Ho: Variances of breeds dry matter intake are equal
leveneTest (cow_measures$DMI~cow_measures$Breed)

#T test
#Ho: mean DMI Holstein = mean DMI Jersey 
#Two-sided t test
#default options used: mu = 0, alt = "two.sided", conf = 0.95, paired = F

t.test(cow_measures$DMI~cow_measures$Breed, var.eq=T)

#Levene test
#Ho: Variances of breeds milk yield are equal
leveneTest (cow_measures$MY~cow_measures$Breed)

#T test
#Ho: mean MY Holstein = mean MY Jersey 
#Two-sided t test
#default options used: mu = 0, alt = "two.sided", conf = 0.95, paired = F

t.test(cow_measures$MY~cow_measures$Breed, var.eq=T)
```

#Data Curation

Sequences were downloaded from the Torrent server in a fastq file (seqs.fastq) which was used to generate a fasta file (seqs.fna). Initial quality control steps were performed using the Torrent Suite Software (refer to manuscript) 

```{r, engine='bash', eval=FALSE}
convert_fastaqual_fastq.py  -c fastq_to_fastaqual –f seqs.fastq –o fastaqual
```

Because of file size limit in GitHub, the previously described seqs.fastq and seqs.fna files were not uploaded.

#Demulitplex and Quality Control

Demultiplex the sequences in the library using the mapping file.

```{r, engine='bash', eval=FALSE}

split_libraries.py -f seqs.fna -b variable_length -l 0 -L 1000 -x -M 1  -m mappingfile.txt -o seqs_demultiplexed
```

The seqs_demultiplexed.fna generated is provided. Remove primer and subsequent sequence.

```{r, engine='bash'}
wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/seqs_demultiplexed.fna.tar.gz 

wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/mappingfile.txt

tar -zxvf seqs_demultiplexed.fna.tar.gz

truncate_reverse_primer.py -f seqs_demultiplexed.fna -m mappingfile.txt -z truncate_only -M 2 -o seqs_revprimer_truncated
```

Trim the sequences to a fixed length (130 basepairs for this study) using 2 custom perl scripts to improve OTU selection in UPARSE downstream.  

```{r, engine='bash'}
wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/Scripts/min_max_length.pl 

chmod 775 min_max_length.pl

wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/Scripts/second_min_max_length.pl

chmod 775 second_min_max_length.pl

./min_max_length.pl -min=131 -max=131 -fasta=seqs_revprimer_truncated/seqs_demultiplexed_rev_primer_truncated.fna

./second_min_max_length.pl -min=130 -max=130 -fasta=seqs_trimmed.fasta
```

Reverse complement the sequences in mothur.

```{r, engine='bash'}
mothur "#reverse.seqs(fasta=seqs_finaltrimmed.fasta)"
```

#Pick OTUs, assign taxonomy, align sequences, and generate phylogenetic tree

Use a custom perl script to convert the fasta file from QIIME format to UPARSE format to generate the OTU table.

```{r, engine='bash', results='hide'}
wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/Scripts/qiime_to_usearch.pl

chmod 775 qiime_to_usearch.pl

./qiime_to_usearch.pl -fasta=seqs_finaltrimmed.rc.fasta -prefix=cow

mv format.fasta dairybreeds.format.fasta
```

Run the sequences through the UPARSE pipeline to pick OTUs.

```{r, engine='bash'}
svn export https://github.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/trunk/Scripts/usearch_python_scripts --non-interactive --trust-server-cert

chmod -R 775 usearch_python_scripts

chmod 775 usearch_python_scripts/uc2otutab.py

chmod 775 usearch_python_scripts/fasta_number.py

wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/gold.fasta.gz 
gzip -d gold.fasta.gz

chmod 775 gold.fasta

mkdir usearch_results

usearch -derep_fulllength dairybreeds.format.fasta -sizeout -fastaout usearch_results/dairybreeds.derep.fasta

usearch -sortbysize usearch_results/dairybreeds.derep.fasta -minsize 2 -fastaout usearch_results/dairybreeds.derep.sort.fasta

usearch -cluster_otus usearch_results/dairybreeds.derep.sort.fasta -otus usearch_results/dairybreeds.otus1.fasta

usearch -uchime_ref usearch_results/dairybreeds.otus1.fasta -db gold.fasta -strand plus -nonchimeras usearch_results/dairybreeds.otus1.nonchimera.fasta

python usearch_python_scripts/fasta_number.py usearch_results/dairybreeds.otus1.nonchimera.fasta > usearch_results/dairybreeds.otus2.fasta

usearch -usearch_global dairybreeds.format.fasta -db usearch_results/dairybreeds.otus2.fasta -strand plus -id 0.97 -uc usearch_results/dairybreeds.otu_map.uc

python usearch_python_scripts/uc2otutab.py usearch_results/dairybreeds.otu_map.uc > usearch_results/dairybreeds.otu_table.txt

cp usearch_results/dairybreeds.otu_table.txt ./
```

Assign taxonomy to the OTU representative sequences:

```{r, engine='bash'}
wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/gg_13_5_otus.tar.gz  

tar -zxvf gg_13_5_otus.tar.gz

assign_taxonomy.py -i usearch_results/dairybreeds.otus2.fasta -t gg_13_5_otus/97_otu_taxonomy.txt -r gg_13_5_otus/97_otus.fasta -o assigned_gg_taxa
```

Add the taxa outputted to the OTU table with the column header "taxonomy" and output the resulting file to biom format:

```{r, engine='bash'}
awk 'NR==1; NR > 1 {print $0 | "sort"}' dairybreeds.otu_table.txt > dairybreeds.otu_table.sort.txt 
sort assigned_gg_taxa/dairybreeds.otus2_tax_assignments.txt > assigned_gg_taxa/dairybreeds.otus2_tax_assignments.sort.txt
{ printf '\ttaxonomy\t\t\n'; cat assigned_gg_taxa/dairybreeds.otus2_tax_assignments.sort.txt ; }  > assigned_gg_taxa/dairybreeds.otus2_tax_assignments.sort.label.txt
awk '!/2496/{print}' assigned_gg_taxa/dairybreeds.otus2_tax_assignments.sort.label.txt > assigned_gg_taxa/dairybreeds.otus2_tax_assignments.sort.label2.txt

paste dairybreeds.otu_table.sort.txt <(cut -f 2 assigned_gg_taxa/dairybreeds.otus2_tax_assignments.sort.label2.txt) > dairybreeds.otu_table.tax.txt

rm dairybreeds.otu_table.sort.txt

biom convert --table-type "OTU table" -i dairybreeds.otu_table.tax.txt -o dairybreeds.otu_table.tax.biom --process-obs-metadata taxonomy --to-json
```

Remove sample x106(Jersey cow). Sample was not obtained from this cow. 

```{r, engine='bash'}
printf "x106" > remove_sample.txt

filter_samples_from_otu_table.py -i dairybreeds.otu_table.tax.biom -o dairybreeds.otu_table.tax.filtered.biom --sample_id_fp remove_sample.txt --negate_sample_id_fp
```

Sequences were aligned using the [RDP aligner](https://pyro.cme.msu.edu/login.spr). The aligned_dairybreeds.otus2.fasta_alignment_summary.txt obtained from the alignment is provided. This file was used to assess OTUs that did not aligned properly based on the following criteria:

  1. The minimum starting position was 324. As previously described, sequences were trimmed to be 130 basepairs in length and since the RDP aligner sets the end position to be 454, sequences must start at 324 to be 130 basepairs (i.e. 454-324=130).
  2. The maximun starting position was 354. The minimum number of basepairs accepted was 100, thus these sequences should start at 354 to contain 100 basepairs (i.e. 454-354=100).

Output files from the RDP aligner are provided.

```{r, engine='bash', results='hide'}
wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/aligned_dairybreeds.otus2.fasta

wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/aligned_dairybreeds.otus2.fasta_alignment_summary.txt 
```

Make a file containing all the OTUs that aligned properly

```{r}
alignedOTUs <- read.table("aligned_dairybreeds.otus2.fasta_alignment_summary.txt", header=T, sep="\t")

properOTUtable <- subset(alignedOTUs, (Start >= 324 & Start <= 354) & End == 454, select = SequenceId)

write.table (properOTUtable, file = "proper_aligned_otu.txt", col.names = F, row.names = F)
```

Remove those OTUs that did not align well from the OTU table. Then remove OTUs with Cyanobacteria classification. UPARSE should have removed sinlgeton OTUs, but double check  with the (-n 2 parameter):

```{r, engine='bash'}
filter_otus_from_otu_table.py -i dairybreeds.otu_table.tax.filtered.biom -o dairybreeds.otu_table.tax.filtered.filtered.biom --negate_ids_to_exclude -e proper_aligned_otu.txt -n 2

filter_taxa_from_otu_table.py -i dairybreeds.otu_table.tax.filtered.filtered.biom -o dairybreeds.otu_table.tax.final.biom -n p__Cyanobacteria
```

Generate quality filtered OTU table with assigend taxonomy 

```{r, engine='bash'}
biom convert -i dairybreeds.otu_table.tax.final.biom -o Table_OTUs.txt --to-tsv --header-key taxonomy
```

```{r}
otus_data <- read.table("Table_OTUs.txt", header=F, sep="\t")
colnames(otus_data) <- c("OTU ID", "Holstein5", "Holstein3", "Jersey3", "Holstein1", "Jersey4", "Holstein3C", "Holstein2", "Holstein4C", "Holstein4", "Holstein5C", "Jersey1", "Holstein2C", "Jersey2", "Holstein1C", "Taxonomy")
otus_data_arranged <- otus_data[c("OTU ID", "Holstein1", "Holstein2", "Holstein3", "Holstein4", "Holstein5", "Holstein1C", "Holstein2C", "Holstein3C", "Holstein4C", "Holstein5C", "Jersey1", "Jersey2", "Jersey3", "Jersey4", "Taxonomy")]
write.table(otus_data_arranged, sep = "\t", file = "TableS3.txt", col.names = T, row.names = F)
```

Use the aligned file to generate a phylogenetic tree using clearcut in mothur. Note that using the unfiltered aligned file does not affect downstream results. Clearcut requires ID lengths greater than ~10 characters, thus add 10 'A's to the front of all sequence names. Then remove the 'A's after the tree is formed. 

```{r, engine='bash'}
sed -i -e 's/>/>AAAAAAAAAA/g' aligned_dairybreeds.otus2.fasta

mothur "#dist.seqs(fasta=aligned_dairybreeds.otus2.fasta, output=lt)"

mothur "#clearcut(phylip=aligned_dairybreeds.otus2.phylip.dist)"

sed -i -e 's/AAAAAAAAAA//g' aligned_dairybreeds.otus2.phylip.tre
```

#Alpha diversity and rarefactions curves

Rarefy quality filtered OTU table to the lowest depth across samples and compute alpha diversity metrics and generate rarefaction plots. 
Note from QIIME: "If the lines for some categories do not extend all the way to the right end of the x-axis, that means that at least one of the samples in that category does not have that many sequences."

```{r, engine='bash'}
biom summarize-table -i dairybreeds.otu_table.tax.final.biom -o dairybreeds.otu_table.tax.final.summary.txt
```

```{r, engine='bash', eval=FALSE}
single_rarefaction.py -i dairybreeds.otu_table.tax.final.biom -d 12141 -o dairybreeds.otu_table_rarefied.biom 
```

```{r, engine='bash'}
wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/dairybreeds.otu_table_rarefied.biom
```

```{r, engine='bash'}
alpha_diversity.py -i dairybreeds.otu_table_rarefied.biom -m chao1,observed_otus,goods_coverage -o alpha_rarefaction_even.txt

wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/qiime_parameters.txt

alpha_rarefaction.py -i dairybreeds.otu_table_rarefied.biom -n 10 --min_rare_depth 1 -m mappingfile.txt -p qiime_parameters.txt -o alpha_rarefaction_plots_even
```

Generate alpha diversity charts

```{r}
library(XML)
library(matrixStats)

rare_table <- readHTMLTable("alpha_rarefaction_plots_even/alpha_rarefaction_plots/rarefaction_plots.html")
rare_table$rare_data[rare_table$rare_data == "nan"] <- NA
alpha_rare <- na.omit(rare_table$rare_data)
colnames(alpha_rare)[2] <- "Seqs.Sample"
colnames(alpha_rare)[3] <- "chao1.Ave."
colnames(alpha_rare)[4] <- "chao1.Err."
colnames(alpha_rare)[5] <- "observed_otus.Ave."
colnames(alpha_rare)[6] <- "observed_otus.Err."

cols = c(2, 3, 4, 5, 6)
alpha_rare[, cols] <- lapply(cols, function(x) as.numeric(as.character(alpha_rare[, x])))
breed_rare <- subset(alpha_rare, Description %in% c("Holstein", "HolsteinC", "Jersey"))
method_rare <- subset(alpha_rare, Description %in% c("Cannula", "TubingH"))
breed_rare$label <- "Breed"
method_rare$label <- "Sampling Method"
alpha_rare <- rbind(breed_rare, method_rare)
pd <- position_dodge(width = 275)

rare_otu_plot <- ggplot(alpha_rare, aes(x = Seqs.Sample, y = observed_otus.Ave., colour = Description, group = Description, ymin = observed_otus.Ave. - observed_otus.Err., ymax = observed_otus.Ave. + observed_otus.Err.)) + geom_line(position = pd) + geom_pointrange(position = pd) + scale_colour_manual(values = c("#FF0000", "#0000FF", "#FFA500", "#006400", "#FF00FF"), limits=c("Holstein", "HolsteinC", "Jersey", "Cannula", "TubingH"), labels=c("Holstein", "HolsteinC", "Jersey", "Cannula", "Tubing")) + labs(x = "Sequences per Sample", y = "Mean Observed OTUs") + theme(legend.title = element_blank(), axis.text = element_text(color = "black", size = 8, face = "bold"), axis.title = element_text(color = "black", size = 12, face = "bold"), legend.text = element_text(face = "bold"), strip.text = element_text(size=12, face = "bold")) + facet_grid(~label)

rare_chao1_plot <- ggplot(alpha_rare, aes(x = Seqs.Sample, y = chao1.Ave., colour = Description, group = Description, ymin = chao1.Ave. - chao1.Err., ymax = chao1.Ave. + chao1.Err.)) + geom_line(position = pd) + geom_pointrange(position = pd) + scale_colour_manual(values = c("#FF0000", "#0000FF", "#FFA500", "#006400", "#FF00FF"), limits=c("Holstein", "HolsteinC", "Jersey", "Cannula", "TubingH"), labels=c("Holstein", "HolsteinC", "Jersey", "Cannula", "Tubing")) + labs(x = "Sequences per Sample", y = "Mean Chao1") + theme(legend.title = element_blank(), axis.text = element_text(color = "black", size = 8, face = "bold"), axis.title = element_text(color = "black", size = 12, face = "bold"), legend.text = element_text(face = "bold"), strip.text = element_text(size=12, face = "bold")) + facet_grid(~label)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

png("FigS2.png", units = "in", height = 6, width = 12, res = 300)
multiplot(rare_chao1_plot, rare_otu_plot, cols=2)
dev.off()

pdf("FigS2.pdf", height = 6, width = 12)
multiplot(rare_chao1_plot, rare_otu_plot, cols=2)
dev.off()
```

Perform a two-sided t-test to evaluate alpha diversity metrics

```{r}
alpha <- read.table ("alpha_rarefaction_even.txt", header = T, sep="\t")
alpha_sub <- alpha[,c(2,3)]
breed_data <- c("Holstein", "Holstein", "Jersey", "Holstein", "Jersey", "Holstein", "HolsteinC", "HolsteinC", "Holstein", "HolsteinC", "Jersey", "HolsteinC", "Jersey", "HolsteinC") 
alpha_metrics <- data.frame(alpha_sub, breed_data)
breed_final <- subset(alpha_metrics, breed_data != "HolsteinC")
method_final <- subset(alpha_metrics, breed_data != "Jersey")

#Levene test
#Ho: Variances of breeds are equal
leveneTest (breed_final$chao1~breed_final$breed_data)
leveneTest (breed_final$observed_otus~breed_final$breed_data)

#T test
#Ho: alpha metrics Holstein = alpha metrics Jersey 
#Two-sided t test
#default options used: mu = 0, alt = "two.sided", conf = 0.95, paired = F

t.test(breed_final$chao1~breed_final$breed_data, var.eq=T)
t.test(breed_final$observed_otus~breed_final$breed_data, var.eq=T)

#Levene test
#Ho: Variances of methods are equal
leveneTest (method_final$chao1~method_final$breed_data)
leveneTest (method_final$observed_otus~method_final$breed_data)

#T test
#Ho: alpha metrics Tubing = alpha metrics Cannula 
#Two-sided t test
#default options used: mu = 0, alt = "two.sided", conf = 0.95, paired = F

t.test(method_final$chao1~method_final$breed_data, var.eq=T)
t.test(method_final$observed_otus~method_final$breed_data, var.eq=T)
```

#Taxonomy

Assign taxonomy and generate plots for desired taxa. 

```{r, engine='bash'}
summarize_taxa.py -i dairybreeds.otu_table_rarefied.biom -o summarize_taxa -L 2,3,4,5,6,7
```

Generate a stacked bar chart for the phylum taxonomic rank across samples.

```{r}
taxa_data <- read.table("summarize_taxa/dairybreeds.otu_table_rarefied_L2.txt", header=F, sep="\t")
colnames(taxa_data) <- c("Phylum", "Holstein5", "Holstein3", "Jersey3", "Holstein1", "Jersey4", "Holstein2", "Holstein3C", "Holstein4C", "Holstein4", "Holstein5C", "Jersey1", "Holstein2C", "Jersey2", "Holstein1C") 
taxa_data2 <- taxa_data[c("Phylum", "Holstein1", "Holstein1C", "Holstein2", "Holstein2C", "Holstein3", "Holstein3C", "Holstein4", "Holstein4C", "Holstein5", "Holstein5C", "Jersey1", "Jersey2", "Jersey3", "Jersey4")]
taxa_data2$Phylum <- sub("k__", "", taxa_data2$Phylum)
taxa_data2$Phylum <- sub("Bacteria;", "", taxa_data2$Phylum)
taxa_data2$Phylum <- sub("p__", "", taxa_data2$Phylum)
taxa_data2$Phylum <- sub("Unassigned;Other", "Unassigned", taxa_data2$Phylum)
taxa_data2$Phylum <- sub("Other", "No Assigned Phylum", taxa_data2$Phylum)
plot_taxa_long <- gather(taxa_data2, Samples, Proportion, Holstein1:Jersey4)

png("FigS3.png", units = "in", height = 12, width = 17, res = 300)
ggplot(plot_taxa_long, aes(x=Samples, y=Proportion, fill=Phylum)) + geom_bar(stat="identity") + scale_fill_manual(values= c("#FF0000", "#0000FF", "#FFA500", "#556B2F", "#800080", "#FFFF00", "#00FFFF", "#FFC0CB", "#5F9EA0", "#A52A2A", "#808080", "#00FF00", "#000080", "#ADD8E6", "#E9967A", "#00FA9A", "#DA70D6", "#FFD700")) + theme(axis.text = element_text(color = "black", size = 11, face = "bold"), axis.title = element_text(color = "black", size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12))
dev.off()                   

pdf("FigS3.pdf", height = 12, width = 17)
ggplot(plot_taxa_long, aes(x=Samples, y=Proportion, fill=Phylum)) + geom_bar(stat="identity") + scale_fill_manual(values= c("#FF0000", "#0000FF", "#FFA500", "#556B2F", "#800080", "#FFFF00", "#00FFFF", "#FFC0CB", "#5F9EA0", "#A52A2A", "#808080", "#00FF00", "#000080", "#ADD8E6", "#E9967A", "#00FA9A", "#DA70D6", "#FFD700")) + theme(axis.text = element_text(color = "black", size = 11, face = "bold"), axis.title = element_text(color = "black", size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12))
dev.off()
```

#Shared OTUs and Sequences

Generate venn diagrams with the number of shared OTUs 

```{r, engine='bash'}
mkdir venns

split_otu_table.py -i dairybreeds.otu_table_rarefied.biom -o split_by_breed -m mappingfile.txt -f Breed

biom convert -i split_by_breed/dairybreeds.otu_table_rarefied__Breed_Holstein__.biom -o split_by_breed/dairybreeds.otu_table_rarefied__Breed_Holstein__json.biom --to-json --table-type="OTU table"

biom convert -i split_by_breed/dairybreeds.otu_table_rarefied__Breed_Jersey__.biom -o split_by_breed/dairybreeds.otu_table_rarefied__Breed_Jersey__json.biom --to-json --table-type="OTU table"

collapse_samples.py -b dairybreeds.otu_table_rarefied.biom -m mappingfile.txt --output_biom_fp venns/dairybreeds.otu_table_rarefied_breed.biom --output_mapping_fp venns/mappingfile_breeds.txt --collapse_mode mean --collapse_fields Breed

biom convert -i venns/dairybreeds.otu_table_rarefied_breed.biom -o venns/dairybreeds.otu_table_rarefied_breed_json.biom --to-json --table-type="OTU table"
```

```{r}
library(biom)
library(gplots)

png("FigS4a.png", units = "in", height = 6, width = 8, res = 300)
holstein_biom <- read_biom("split_by_breed/dairybreeds.otu_table_rarefied__Breed_Holstein__json.biom")
holstein_df <- as.data.frame(as(biom_data(holstein_biom), "matrix"))
colnames(holstein_df) <- c("Holstein5", "Holstein3", "Holstein1", "Holstein2", "Holstein4")
holstein_boolean_df <- as.data.frame(holstein_df > 0 + 0)
holstein_venn <- venn(holstein_boolean_df)
dev.off()

png("FigS4b.png", units = "in", height = 9, width = 10, res = 300)
jersey_biom <- read_biom("split_by_breed/dairybreeds.otu_table_rarefied__Breed_Jersey__json.biom")
jersey_df <- as.data.frame(as(biom_data(jersey_biom), "matrix"))
colnames(jersey_df) <- c("Jersey3", "Jersey4", "Jersey1", "Jersey2")
jersey_boolean_df <- as.data.frame(jersey_df > 0 + 0)
jersey_venn <- venn(jersey_boolean_df)
dev.off()

png("FigS4c.png", units = "in", height = 6, width = 6, res = 300)
breed_biom <- read_biom("venns/dairybreeds.otu_table_rarefied_breed_json.biom")
breed_df <- as.data.frame(as(biom_data(breed_biom), "matrix"))
breed_sub <- breed_df[,-2]
breed_boolean_df <- as.data.frame(breed_sub > 0 + 0)
breed_venn <- venn(breed_boolean_df)
dev.off()

png("FigS4d.png", units = "in", height = 6, width = 6, res = 300)
method_biom <- read_biom("venns/dairybreeds.otu_table_rarefied_breed_json.biom")
method_df <- as.data.frame(as(biom_data(method_biom), "matrix"))
method_sub <- method_df[,-3]
colnames(method_sub) <- c("Tubing", "Cannula")
method_boolean_df <- as.data.frame(method_sub > 0 + 0)
method_venn <- venn(method_boolean_df)
dev.off()
```

Calculate the number of shared sequences

```{r, engine='bash'}
mkdir sharedseqs
```

```{r}
#Holstein shared sequences
collapse_biom <- read_biom("split_by_breed/dairybreeds.otu_table_rarefied__Breed_Holstein__json.biom")
collapse <- as.matrix(as(biom_data(collapse_biom), "matrix"))
collapse_df <- as.data.frame(collapse)

seq_shared_func <- function(x) {
  single_combo <- unlist(x)
  collapse_sub <- collapse_df[, names(collapse_df) %in% single_combo]
  collapse_sub[, 3] <- collapse_sub[, 1] + collapse_sub[, 2]
  sub_sum <- colSums(collapse_sub)
  collapse_sub[collapse_sub == 0] <- NA
  collapse_sub2 <- na.omit(collapse_sub)
  collapse_sub2[, 3] <- collapse_sub2[, 1] + collapse_sub2[, 2]
  sub_sum2 <- colSums(collapse_sub2)
  per <- (sub_sum2["V3"]/sub_sum["V3"]) * 100
  collapse_out <- c(names(collapse_sub2)[1], names(collapse_sub2)[2], toString(sub_sum2["V3"]), 
                    toString(sub_sum["V3"]), toString(per))
  write(collapse_out, file = "sharedseqs/Holsteinsharedseqs.txt", sep = "\t", ncolumns = 5, append = TRUE)
}

combn(colnames(collapse), 2, simplify = FALSE, FUN = seq_shared_func)

#Jersey shared sequences
collapse_biom <- read_biom("split_by_breed/dairybreeds.otu_table_rarefied__Breed_Jersey__json.biom")
collapse <- as.matrix(as(biom_data(collapse_biom), "matrix"))
collapse_df <- as.data.frame(collapse)

seq_shared_func <- function(x) {
  single_combo <- unlist(x)
  collapse_sub <- collapse_df[, names(collapse_df) %in% single_combo]
  collapse_sub[, 3] <- collapse_sub[, 1] + collapse_sub[, 2]
  sub_sum <- colSums(collapse_sub)
  collapse_sub[collapse_sub == 0] <- NA
  collapse_sub2 <- na.omit(collapse_sub)
  collapse_sub2[, 3] <- collapse_sub2[, 1] + collapse_sub2[, 2]
  sub_sum2 <- colSums(collapse_sub2)
  per <- (sub_sum2["V3"]/sub_sum["V3"]) * 100
  collapse_out <- c(names(collapse_sub2)[1], names(collapse_sub2)[2], toString(sub_sum2["V3"]), 
                    toString(sub_sum["V3"]), toString(per))
  write(collapse_out, file = "sharedseqs/Jerseysharedseqs.txt", sep = "\t", ncolumns = 5, append = TRUE)
}

combn(colnames(collapse), 2, simplify = FALSE, FUN = seq_shared_func)

#Breeds and Methods shared sequences
collapse_biom <- read_biom("venns/dairybreeds.otu_table_rarefied_breed_json.biom")
collapse <- as.matrix(as(biom_data(collapse_biom), "matrix"))
collapse_df <- as.data.frame(collapse)

seq_shared_func <- function(x) {
  single_combo <- unlist(x)
  collapse_sub <- collapse_df[, names(collapse_df) %in% single_combo]
  collapse_sub[, 3] <- collapse_sub[, 1] + collapse_sub[, 2]
  sub_sum <- colSums(collapse_sub)
  collapse_sub[collapse_sub == 0] <- NA
  collapse_sub2 <- na.omit(collapse_sub)
  collapse_sub2[, 3] <- collapse_sub2[, 1] + collapse_sub2[, 2]
  sub_sum2 <- colSums(collapse_sub2)
  per <- (sub_sum2["V3"]/sub_sum["V3"]) * 100
  collapse_out <- c(names(collapse_sub2)[1], names(collapse_sub2)[2], toString(sub_sum2["V3"]), 
                    toString(sub_sum["V3"]), toString(per))
  write(collapse_out, file = "sharedseqs/breedssharedseqs.txt", sep = "\t", ncolumns = 5, append = TRUE)
}

combn(colnames(collapse), 2, simplify = FALSE, FUN = seq_shared_func)
```

#Beta Diversity

Generate core file of OTUS. To include an OTU as a member of the core, a cut off of 80% (4/5 cows) was assigned for Holstein cows and 75% (3/4 cows) for Jersey cows.

```{r, engine='bash'}
mkdir cores 

filter_otus_from_otu_table.py -i split_by_breed/dairybreeds.otu_table_rarefied__Breed_Holstein__.biom -o cores/core_Holstein.biom -s 4 

filter_otus_from_otu_table.py -i split_by_breed/dairybreeds.otu_table_rarefied__Breed_HolsteinC__.biom -o cores/core_HolsteinC.biom -s 4 

filter_otus_from_otu_table.py -i split_by_breed/dairybreeds.otu_table_rarefied__Breed_Jersey__.biom -o cores/core_Jersey.biom -s 3

merge_otu_tables.py -i cores/core_Holstein.biom,cores/core_HolsteinC.biom,cores/core_Jersey.biom -o cores/merged_cores.biom

biom convert -i cores/merged_cores.biom -o cores/merged_cores.txt --to-tsv
```

When cores files are merged, values are set to 0 for an OTU that is not a core in any other of the merged files. However, this does not mean that the OTU was not present but that it did not meet the cutoff criteria in a particular core file. Reads across OTUs have to be adjusted in the merged core file. From the adjusted core file, filter the OTU table into by breed and sampling method to generate respective beta diversity outputs. Additionally, from a file with all the data generate beta diversity outputs.

```{r}
cores <- read.table("cores/merged_cores.txt", sep = "\t", header = F)
cores_sub <- cores [,1]
write.table(cores_sub, file = "cores/core_otus.txt", col.names = F, row.names = F)
```

```{r, engine='bash'}
filter_otus_from_otu_table.py -i dairybreeds.otu_table_rarefied.biom -o final_core.biom --negate_ids_to_exclude -e cores/core_otus.txt 

filter_samples_from_otu_table.py -i final_core.biom -o filter_breed.biom -m mappingfile.txt -s 'Breed:*,!HolsteinC'

filter_samples_from_otu_table.py -i final_core.biom -o filter_method.biom -m mappingfile.txt -s 'Breed:*,!Jersey'

beta_diversity_through_plots.py -i final_core.biom  -o beta_div_all -t aligned_dairybreeds.otus2.phylip.tre -m mappingfile.txt -p qiime_parameters.txt

beta_diversity_through_plots.py -i filter_breed.biom -o beta_div_breed -t aligned_dairybreeds.otus2.phylip.tre -m mappingfile.txt -p qiime_parameters.txt 

beta_diversity_through_plots.py -i  filter_method.biom -o beta_div_method -m mappingfile.txt -t aligned_dairybreeds.otus2.phylip.tre -p qiime_parameters.txt 
```

Generate the principal coordinate analyses plots

```{r}
#breed and sampling method
all_data <- read.table("beta_div_all/weighted_unifrac_pc.txt", skip = 9, nrows = 14, sep = "\t")
all_pc <- all_data[, c("V2", "V3")]
samples_all <- c("H5", "H3", "J3", "H1", "J4", "H2", "H3C", "H4C", "H4", "H5C", "J1", "H2C", "J2", "H1C")
breed_all <- c("Holstein", "Holstein", "Jersey", "Holstein", "Jersey", "Holstein", "HolsteinC", "HolsteinC", "Holstein", "HolsteinC", "Jersey", "HolsteinC", "Jersey", "HolsteinC")
pc_all <- data.frame (samples_all, all_pc, breed_all)
colnames(pc_all) <- c("SampleID", "PC1", "PC2", "Breed")

all_pc_plot <- ggplot(pc_all, aes(x = PC1, y = PC2, shape = Breed)) + geom_point(size=3) + labs(x = "PC1 (31.8%)", y = "PC2 (16.5%)") + ggtitle("Breed and Sampling Method") + theme(title = element_text(color = "black", size = 14, face = "bold"), axis.text = element_text(color = "black", size = 14, face = "bold"), axis.title = element_text(color = "black", size = 14, face = "bold"), legend.text = element_text(size = 14, face = "bold"), legend.title = element_blank()) + geom_text(aes(label = SampleID), size = 4, vjust = -0.5, hjust = -0.3)

#breed
breed_data <- read.table("beta_div_breed/weighted_unifrac_pc.txt", skip = 9, nrows = 9, sep = "\t")
breed_pc <- breed_data[, c("V2", "V3")]
samples <- c("H5", "H3", "J3", "H1", "J4", "H2", "H4", "J1", "J2")
breed <- c("Holstein", "Holstein", "Jersey", "Holstein", "Jersey", "Holstein", "Holstein", "Jersey", "Jersey")
pc_breed <- data.frame (samples, breed_pc, breed)
colnames(pc_breed) <- c("SampleID", "PC1", "PC2", "Breed")

breed_pc_plot <- ggplot(pc_breed, aes(x = PC1, y = PC2, shape = Breed)) + geom_point(size=3) + labs(x = "PC1 (36.9%)", y = "PC2 (20.3%)") + ggtitle("Breed") + theme(title = element_text(color = "black", size = 14, face = "bold"), axis.text = element_text(color = "black", size = 14, face = "bold"), axis.title = element_text(color = "black", size = 14, face = "bold"), legend.text = element_text(size = 14, face = "bold"), legend.title = element_blank()) + geom_text(aes(label = SampleID), size = 4, vjust = -0.5, hjust = -0.3)

#sampling method
method_data <- read.table("beta_div_method/weighted_unifrac_pc.txt", skip = 9, nrows = 10, sep = "\t")
method_pc <- method_data[, c("V2", "V3")]
samples_method <- c("H5", "H3", "H1", "H2", "H3C", "H4C", "H4", "H5C", "H2C", "H1C")
method <- c("Tubing", "Tubing", "Tubing", "Tubing", "Cannula", "Cannula", "Tubing", "Cannula", "Cannula", "Cannula")
pc_method <- data.frame (samples_method, method_pc, method)
colnames(pc_method) <- c("SampleID", "PC1", "PC2", "Method")

method_pc_plot <- ggplot(pc_method, aes(x = PC1, y = PC2, shape = Method)) + geom_point(size=3) + labs(x = "PC1 (32.4%)", y = "PC2 (24.3%)") + ggtitle("Sampling Method") + theme(title = element_text(color = "black", size = 14, face = "bold"), axis.text = element_text(color = "black", size = 14, face = "bold"), axis.title = element_text(color = "black", size = 14, face = "bold"), legend.text = element_text(size = 14, face = "bold"), legend.title = element_blank()) + geom_text(aes(label = SampleID), size = 4, vjust = -0.5, hjust = -0.3)

png("Fig2.png", height = 8, width = 26, units = "in", res = 300)
multiplot(all_pc_plot, breed_pc_plot, method_pc_plot, cols = 3)
dev.off()

pdf("Fig2.pdf", height = 8, width = 26)
multiplot(all_pc_plot, breed_pc_plot, method_pc_plot, cols = 3)
dev.off()
```

#Community Statistics

Evaluate the effect of breed and method on microbial community composition using the permutational multivariate analysis of variance (permanova) test (implemented in R as adonis function in the vegan package) using the unweighted unifrac distance matrix

```{r}
library(vegan)

#Mapping files
mappingfile <- read.table("mappingfile.txt", sep = "\t", header = F)
colnames(mappingfile) <- c("SampleID", "BarcodeSequence", "LinkerPrimerSequence", "ReversePrimer", "Method", "Breed", "BCcount", "Well", "Description")
mappingcorrect <- mappingfile[-3,]
breed_map <- subset(mappingcorrect, Breed != "HolsteinC")
method_map <- subset(mappingfile, Breed != "Jersey") 

breed_weighted <- read.table("beta_div_breed/weighted_unifrac_dm.txt", sep = "\t", header = TRUE)
#Order breed mapping file Sample IDs like matirx Sample IDs
breed_map = breed_map[match(breed_weighted$X, breed_map$SampleID),]
row.names(breed_weighted) <- breed_weighted$X
breed_weighted <- breed_weighted[, -1]
breed_weighted <- as.dist(breed_weighted)

method_weighted <- read.table("beta_div_method/weighted_unifrac_dm.txt", sep = "\t", header = TRUE)
#Order method mapping file Sample IDs like matirx Sample IDs
method_map = method_map[match(method_weighted$X, method_map$SampleID),]
row.names(method_weighted) <- method_weighted$X
method_weighted <- method_weighted[, -1]
method_weighted <- as.dist(method_weighted)

all_weighted <- read.table("beta_div_all/weighted_unifrac_dm.txt", sep = "\t", header = TRUE)
#Order method mapping file Sample IDs like matirx Sample IDs
mappingcorrect = mappingcorrect[match(all_weighted$X, mappingcorrect$SampleID),]
row.names(all_weighted) <- all_weighted$X
all_weighted <- all_weighted[, -1]
all_weighted <- as.dist(all_weighted)

adonis(breed_weighted ~ Breed, permutations = 999, data = breed_map)
adonis(method_weighted ~ Method, permutations = 999, data = method_map)
adonis(all_weighted ~ Breed, permutations = 999, data = mappingcorrect)
```

HOMOVA

```{r, engine='bash'}
wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/breed.dist.txt

wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/breed.design.txt

wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/method.dist.txt

wget https://raw.githubusercontent.com/enriquepaz/2016_Paz_et_al_Dairy_Breeds/master/method.design.txt
```

```{r, engine='bash'}
mothur "#homova(phylip=breed.dist.txt, design=breed.design.txt)"

mothur "#homova(phylip=method.dist.txt, design=method.design.txt)"
```

Make heat map of all core OTUs

```{r, engine='bash', results='hide'}
mkdir heatmap

biom convert -i final_core.biom -o heatmap/final_core_json.biom --to-json --table-type="OTU table"
```

```{r}
library(Heatplus)
library(RColorBrewer) 

otu_core <- read_biom("heatmap/final_core_json.biom")
core_table <- as.data.frame(as(biom_data(otu_core), "matrix"))
colnames(core_table) <- c("Holstein5", "Holstein3", "Jersey3", "Holstein1", "Jersey4", "Holstein2", "Holstein3C", "Holstein4C", "Holstein4", "Holstein5C", "Jersey1", "Holstein2C", "Jersey2", "Holstein1C")
core_table_trans <- as.data.frame(t(core_table))
data_propor <- core_table_trans/rowSums(core_table_trans)
scalewhiteblack <- colorRampPalette(c("white", "black"), space = "rgb")(100)
maxab <- apply(data_propor, 2, max)
n1 <- names(which(maxab < 0.01))
data_abun <- data_propor[, -which(names(data_propor) %in% n1)]
data.dist <- vegdist(data_propor, method = "bray")
row.clus <- hclust(data.dist, "aver")
data.dist.g <- vegdist(t(data_abun), method = "bray")
col.clus <- hclust(data.dist.g, "aver")

png("FigS5.png", height = 8, width = 8, units = "in", res = 300)
heatmap.2(as.matrix(data_abun), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scalewhiteblack, margins = c(2, 8), trace = "none", density.info = "none", labCol = "", xlab = "OTUs", ylab = "Samples", lhei = c(2, 8))
dev.off()

pdf("FigS5.pdf", height = 8, width = 8)
heatmap.2(as.matrix(data_abun), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scalewhiteblack, margins = c(2, 8), trace = "none", density.info = "none", labCol = "", xlab = "OTUs", ylab = "Samples", lhei = c(2, 8))
dev.off()
```

LEfSe was used to determine different OTUs between breeds and sampling methods

```{r, engine='bash', results='hide'}
biom convert -i filter_breed.biom -o filter_breed_json.biom --to-json --table-type="OTU table"

biom convert -i filter_method.biom -o filter_method_json.biom --to-json --table-type="OTU table"
```

Create files with LEfSe format

```{r}
breed_biom <- read_biom("filter_breed_json.biom")
breed_table <- as.data.frame(as(biom_data(breed_biom), "matrix"))
breed_rel <- sweep(breed_table, 2, colSums(breed_table), FUN = "/")
rownames(breed_rel) <- paste("OTU", rownames(breed_rel), sep = "")
breed_rel2 <- breed_rel[0, ]
breed_rel2[nrow(breed_rel2) + 1, ] <- c("Holstein", "Holstein", "Jersey", "Holstein", "Jersey", "Holstein", "Holstein", "Jersey", "Jersey")
breed_rel2[nrow(breed_rel2) + 1, ] <- c("Holstein5", "Holstein3", "Jersey3", "Holstein1", "Jersey4", "Holstein2", "Holstein4", "Jersey1", "Jersey2")
row.names(breed_rel2) <- c("breed", "ID")
breed_abundance <- rbind(breed_rel2, breed_rel)
write.table(breed_abundance, sep = "\t", file = "breed_abundance.txt", row.names = TRUE, col.names = FALSE, quote = FALSE)

method_biom <- read_biom("filter_method_json.biom")
method_table <- as.data.frame(as(biom_data(method_biom), "matrix"))
method_rel <- sweep(method_table, 2, colSums(method_table), FUN = "/")
rownames(method_rel) <- paste("OTU", rownames(method_rel), sep = "")
method_rel2 <- method_rel[0, ]
method_rel2[nrow(method_rel2) + 1, ] <- c("Tubing", "Tubing", "Tubing", "Tubing", "Cannula", "Cannula", "Tubing", "Cannula", "Cannula", "Cannula")
method_rel2[nrow(method_rel2) + 1, ] <- c("Holstein5", "Holstein3", "Holstein1", "Holstein2", "Holstein3C", "Holstein4C", "Holstein4", "Holstein5C", "Holstein2C", "Holstein1C")
row.names(method_rel2) <- c("method", "ID")
method_abundance <- rbind(method_rel2, method_rel)
write.table(method_abundance, sep = "\t", file = "method_abundance.txt", row.names = TRUE, col.names = FALSE, quote = FALSE)
```

LEfSe was used to determine OTUs that differ

```{r, engine='bash'}
wget https://bitbucket.org/nsegata/lefse/get/default.zip -O lefse.zip
unzip lefse.zip
mv nsegata* lefse

python lefse/format_input.py breed_abundance.txt breed_lefse_format.txt -c 1 -u 2 -o 1000000
python lefse/run_lefse.py breed_lefse_format.txt breed_lefse_result.txt

python lefse/format_input.py method_abundance.txt method_lefse_format.txt -c 1 -u 2 -o 1000000 
python lefse/run_lefse.py method_lefse_format.txt method_lefse_result.txt
```

```{r}
breed_lefse <- read.table("breed_lefse_result.txt", header = FALSE, sep = "\t")
breed_lefse <- breed_lefse[complete.cases(breed_lefse), ]
breed_lefse$V1 <- gsub("OTU", "", breed_lefse$V1)
write.table(breed_lefse$V1, file = "breed_lefse_otus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

method_lefse <- read.table("method_lefse_result.txt", header = FALSE, sep = "\t")
method_lefse <- method_lefse[complete.cases(method_lefse), ]
method_lefse$V1 <- gsub("OTU", "", method_lefse$V1)
write.table(method_lefse$V1, file = "method_lefse_otus.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

Create input files to generate heat maps

```{r, engine='bash'}
biom convert -i filter_breed.biom -o heatmap/filter_breed.txt --to-tsv --header-key taxonomy

biom convert -i filter_method.biom -o heatmap/filter_method.txt --to-tsv --header-key taxonomy

awk '{gsub("; ","\t",$0); print;}' heatmap/filter_breed.txt | awk  '{gsub("#OTU","OTU",$0); print;}' | cut -f-10,15 | tail -n +2 | awk '{if(NR==1){print $0,"\ttaxonomy"}else{print }}' > heatmap/breed_tax.txt

awk '{gsub("; ","\t",$0); print;}' heatmap/filter_method.txt | awk  '{gsub("#OTU","OTU",$0); print;}' | cut -f-11,16 | tail -n +2 | awk '{if(NR==1){print $0,"\ttaxonomy"}else{print }}' > heatmap/method_tax.txt
```

```{r}
#breed
breed_otu <- read.table("heatmap/breed_tax.txt", header = T, sep = "\t", fill = TRUE)
breed_otu$taxonomy <- sub("f__", "", breed_otu$taxonomy)
breed_otu$taxonomy <- sub("\\]", "", breed_otu$taxonomy)
breed_otu$taxonomy <- sub("\\[", "", breed_otu$taxonomy)
breed_otu$taxonomy <- sub("^$", "No Assigned Family", breed_otu$taxonomy)
colnames(breed_otu) <- c("OTUs", "Holstein5", "Holstein3", "Jersey3", "Holstein1", "Jersey4", "Holstein2", "Holstein4", "Jersey1", "Jersey2", "taxonomy")

OTU.ID <- breed_otu[,1]
taxonomy <- breed_otu[,11]
breed_samples <- breed_otu[,-11]

row.names(breed_samples) <- breed_samples$OTUs
breed_samples <- breed_samples[, -1]
breed_trans <- as.data.frame(t(breed_samples))
breed_propor <- breed_trans/rowSums(breed_trans)
breed_propor_trans <- as.data.frame(t(breed_propor))

breed_proportion <- data.frame(OTU.ID, breed_propor_trans, taxonomy)
write.table(breed_proportion, file = "heatmap/breed_propor.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#method
method_otu <- read.table("heatmap/method_tax.txt", header = T, sep = "\t", fill = TRUE)
method_otu$taxonomy <- sub("f__", "", method_otu$taxonomy)
method_otu$taxonomy <- sub("\\]", "", method_otu$taxonomy)
method_otu$taxonomy <- sub("\\[", "", method_otu$taxonomy)
method_otu$taxonomy <- sub("^$", "No Assigned Family", method_otu$taxonomy)
colnames(method_otu) <- c("OTUs", "Holstein5", "Holstein3", "Holstein1", "Holstein2", "Holstein3C", "Holstein4C", "Holstein4", "Holstein5C", "Holstein2C", "Holstein1C", "taxonomy")

OTU.ID <- method_otu[,1]
taxonomy <- method_otu[,12]
method_samples <- method_otu[,-12]

row.names(method_samples) <- method_samples$OTUs
method_samples <- method_samples[, -1]
method_trans <- as.data.frame(t(method_samples))
method_propor <- method_trans/rowSums(method_trans)
method_propor_trans <- as.data.frame(t(method_propor))

method_proportion <- data.frame(OTU.ID, method_propor_trans, taxonomy)
write.table(method_proportion, file = "heatmap/method_propor.txt", sep = "\t", row.names = F, col.names = T, quote = F)
```

```{r, engine='bash'}
awk '{gsub("OTU.ID","#OTU ID",$0); print;}' heatmap/breed_propor.txt> heatmap/breed_propor_final.txt

awk '{gsub("OTU.ID","#OTU ID",$0); print;}' heatmap/method_propor.txt> heatmap/method_propor_final.txt

biom convert -i heatmap/breed_propor_final.txt -o heatmap/breed_propor_final.biom --to-json --table-type="OTU table" --process-obs-metadata taxonomy

biom convert -i heatmap/method_propor_final.txt -o heatmap/method_propor_final.biom --to-json --table-type="OTU table" --process-obs-metadata taxonomy

filter_otus_from_otu_table.py -i heatmap/breed_propor_final.biom -o heatmap/breed_lefse.biom -e breed_lefse_otus.txt --negate_ids_to_exclude

filter_otus_from_otu_table.py -i heatmap/method_propor_final.biom -o heatmap/method_lefse.biom -e method_lefse_otus.txt --negate_ids_to_exclude

biom convert -i heatmap/breed_lefse.biom -o heatmap/breed_lefse.txt --to-tsv --header-key taxonomy

biom convert -i heatmap/method_lefse.biom -o heatmap/method_lefse.txt --to-tsv --header-key taxonomy

awk '{gsub("#OTU ID","OTUs",$0); print;}' heatmap/breed_lefse.txt> heatmap/breed_lefse_final.txt

awk '{gsub("#OTU ID","OTUs",$0); print;}' heatmap/method_lefse.txt> heatmap/method_lefse_final.txt
```

Make heat maps

```{r}

#breed all 
breed_lefse_all <- read.table("heatmap/breed_lefse_final.txt", header = T, sep = "\t")
breed_lefse_data <- breed_lefse_all[,-11]
row.names(breed_lefse_data) <- breed_lefse_data$OTUs
breed_lefse_data <- breed_lefse_data[,-1]
breed_all_trans <- as.data.frame(t(breed_lefse_data))
scalewhiteblack <- colorRampPalette(c("white", "black"), space = "rgb")(100)
data.dist <- vegdist(breed_all_trans, method = "bray")
row.clus <- hclust(data.dist, "aver")
data.dist.g <- vegdist(t(breed_all_trans), method = "bray")
col.clus <- hclust(data.dist.g, "aver")

png("Fig3a.png", height = 8, width = 8, units = "in", res = 300)
heatmap.2(as.matrix(breed_all_trans), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scalewhiteblack, margins = c(2, 8), trace = "none", density.info = "none", labCol = "", xlab = "OTUs", ylab = "Samples", lhei = c(2, 8))
dev.off()

pdf("Fig3a.pdf", height = 8, width = 8)
heatmap.2(as.matrix(breed_all_trans), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scalewhiteblack, margins = c(2, 8), trace = "none", density.info = "none", labCol = "", xlab = "OTUs", ylab = "Samples", lhei = c(2, 8))
dev.off()

#breed most abundant
breed_lefse <- read.table("heatmap/breed_lefse_final.txt", header = T, sep = "\t")
row.names(breed_lefse) <- breed_lefse$OTUs
breed_lefse <- breed_lefse[, -1]
breed_lefse_tax <- subset(breed_lefse, select = c(taxonomy))

breed_lefse_samples <- breed_lefse[,-10]
breed_lefse_trans <- as.data.frame(t(breed_lefse_samples))
scalewhiteblack <- colorRampPalette(c("white", "black"), space = "rgb")(100)
maxab <- apply(breed_lefse_trans, 2, max)
n1 <- names(which(maxab < 0.01))
breed_abun <- breed_lefse_trans[, -which(names(breed_lefse_trans) %in% n1)]
data.dist <- vegdist(breed_lefse_trans, method = "bray")
row.clus <- hclust(data.dist, "aver")
breeds_abun_dist <- vegdist(t(breed_abun), method = "bray")
col.clus <- hclust(breeds_abun_dist, "aver")

breed_abun_trans <-  as.data.frame(t(breed_abun))
merge_data <- merge (breed_abun_trans, breed_lefse_tax, by = "row.names")
row.names(merge_data) <- merge_data$Row.names
merge_data <- merge_data[,-1]                                   
breed_abun_fam <- subset(merge_data, select = c(taxonomy))

png("Fig3b.png", height = 8, width = 8, units = "in", res = 300)
heatmap.2(as.matrix(breed_abun), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scalewhiteblack, margins = c(10, 7), trace = "none", density.info = "none", labCol = breed_abun_fam$taxonomy, xlab = "Family", ylab = "Samples", lhei = c(2, 8))
dev.off()

pdf("Fig3b.pdf", height = 8, width = 8)
heatmap.2(as.matrix(breed_abun), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scalewhiteblack, margins = c(10, 7), trace = "none", density.info = "none", labCol = breed_abun_fam$taxonomy, xlab = "Family", ylab = "Samples", lhei = c(2, 8))
dev.off()
```

