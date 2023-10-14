# Functional analysis 

The output of RNA-seq differential expression analysis is a list of significant differentially expressed genes (DEGs). To gain greater biological insight on the differentially expressed genes there are various analyses that can be done. While the tools for these functional analyses span a wide variety of techniques, they can loosely be categorized into three main types: over-representation analysis, functional class scoring, and pathway topology [[1](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/raw/master/resources/pathway_tools.pdf)]. 

![Pathway analysis tools](img/pathway_analysis.png)

## Over-representation analysis

There are a plethora of functional enrichment tools that perform some type of "over-representation" analysis by querying databases containing information about gene function and interactions.

These databases typically categorize genes into groups based on shared function, or involvement in a pathway, or presence in a specific cellular location, or other categorizations, e.g. functional pathways, etc. These categorries have been consistently named (controlled vocabulary) based on how the gene has been annotated functionally. 

### Gene Ontology project

One of the most widely-used categorizations is the **Gene Ontology (GO)** established by the Gene Ontology project.

"The Gene Ontology project is a collaborative effort to address the need for consistent descriptions of gene products across databases" [[2](http://geneontology.org/page/introduction-go-resource)]. The [Gene Ontology Consortium](http://geneontology.org/page/go-consortium-contributors-list) maintains the GO terms, and these GO terms are incorporated into gene annotations in many of the popular repositories for animal, plant, and microbial genomes. 

Tools that investigate enrichment of biological functions or interactions often use the Gene Ontology (GO) categorizations, i.e. the GO terms to determine whether any have significantly modified representation in a given list of genes.

To determine whether any GO terms are over-represented, you can determine the probability of having the observed proportion of genes associated with a specific GO term in your gene list based on the proportion of genes associated with the same GO term in the background set.

<img src="img/go_proportions.png" width="600">

<img src="img/go_proportions_table3.png" width="600">

The statistical test that will determine whether something is actually over-represented is the *Hypergeometric test*.

### Hypergeometric testing

Using the example of the first functional category above, hypergeometric distribution is a probability distribution that describes the probability of 25 genes (k) being associated with "Functional category 1", for all genes in our gene list (n=1000), from a population of all of the genes in entire genome (N=13,000) which contains 35 genes (K) associated with "Functional category 1" [[4](https://en.wikipedia.org/wiki/Hypergeometric_distribution)].

The calculation of probability of k successes follows the formula:

![hypergeo](img/hypergeo.png) 

This test will result in an adjusted p-value (after multiple test correction) for each category tested.


#### GO Ontologies

To best use and interpret the results from these functional analysis tools, it is helpful to have a good understanding of the GO terms themselves and their organization. To describe the roles of genes and gene products, GO terms are organized into three independent controlled vocabularies (ontologies) in a species-independent manner: 

- **Biological process:** refers to the biological role involving the gene or gene product, and could include "transcription", "signal transduction", and "apoptosis". A biological process generally involves a chemical or physical change of the starting material or input.

- **Molecular function:** represents the biochemical activity of the gene product, such activities could include "ligand", "GTPase", and "transporter". 

- **Cellular component:** refers to the location in the cell of the gene product. Cellular components could include "nucleus", "lysosome", and "plasma membrane".

Each GO term has a term name (e.g. **DNA repair**) and a unique term accession number (**GO:0005125**), and a single gene product can be associated with many GO terms, since a single gene product "may function in several processes, contain domains that carry out diverse molecular functions, and participate in multiple alternative interactions with other proteins, organelles or locations in the cell" [[3](https://github.com/hbctraining/In-depth-NGS-Data-Analysis-Course/raw/master/resources/go_terms.pdf)]. 

#### GO term hierarchy

The GO ontologies were developed to describe and query biological knowledge with differing levels of information available. To do this, GO ontologies are loosely hierarchical, ranging from general, 'parent', terms to more specific, 'child' terms. The GO ontologies are "loosely" hierarchical since 'child' terms can have multiple 'parent' terms.

Some genes with less information may only be associated with general 'parent' terms or no terms at all, while other genes with a lot of information be associated with many terms.

![Nature Reviews Cancer 7, 23-34 (January 2007)](img/go_heirarchy.jpg)

[Tips for working with GO terms](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003343)


## clusterProfiler

We will be using [clusterProfiler](http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) to perform over-representation analysis on GO terms associated with our list of significant genes. The tool takes as input a significant gene list and a background gene list and performs statistical enrichment analysis using hypergeometric testing. The basic arguments allow the user to select the appropriate organism and GO ontology (BP, CC, MF) to test. 

### Running clusterProfiler

Load the following libraries:

```r

# you may have to install some of these libraries; use 
# BiocManager::install(c("org.Hs.eg.db","clusterProfiler","enrichplot","fgsea"))

library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)
library(tidyverse)
library(fgsea)
```

```r

res_tableOE = readRDS("data/res_tableOE.rds")

res_tableOE_tb <- res_tableOE %>%
data.frame() %>%
rownames_to_column(var="gene") %>%
dplyr::filter(!is.na(log2FoldChange))  %>% as_tibble()

```

To perform the over-representation analysis, we need a list of background genes and a list of significant genes. For our background dataset we will use all genes tested for differential expression (all genes in our results table). For our significant gene list we will use genes with p-adjusted values less than 0.05 (we could include a fold change threshold too if we have many DE genes).

```r
## background set of genes
allOE_genes <- res_tableOE_tb$gene

sigOE = dplyr::filter(res_tableOE_tb, padj < 0.05)
sigOE_genes = sigOE$gene

```
Now we can perform the GO enrichment analysis and save the results:

```r
## Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db, 
                minGSSize = 20,
                maxGSSize = 300,
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
                
## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.csv(cluster_summary, "results/clusterProfiler_Mov10oe.csv")
```
>**NOTE:** The different organisms with annotation databases available to use with clusterProfiler can be found [here](img/orgdb_annotation_databases.png)


![cluster_summary](img/cluster_summary.png)                

### Visualizing clusterProfiler results

clusterProfiler has a variety of options for viewing the over-represented GO terms. We will explore the dotplot, enrichment plot, and the category netplot.

The **dotplot** shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.

```r
## Dotplot 
dotplot(ego, showCategory = 50)
```

**To save the figure,** click on the `Export` button in the RStudio `Plots` tab and `Save as PDF...`. In the pop-up window, change:
- `PDF size` to `8 x 14` to give a figure of appropriate size for the text labels

<img src="img/mov10oe_dotplot.png" width="600">

The next plot is the **enrichment GO plot**, which shows the relationship between the top 50 most significantly enriched GO terms (padj.), by grouping similar terms together. The color represents the p-values relative to the other displayed terms (brighter red is more significant) and the size of the terms represents the number of genes that are significant from our list.

```r
## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms

pwt <- pairwise_termsim(
ego,
method = "JC",
semData = NULL,
showCategory = 50
)

emapplot(pwt, showCategory = 50)
```

**To save the figure,** click on the `Export` button in the RStudio `Plots` tab and `Save as PDF...`. In the pop-up window, change the `PDF size` to `24 x 32` to give a figure of appropriate size for the text labels.

<img src="img/mov10oe_enrichmap.png" width="800">

Finally, the **category netplot** shows the relationships between the genes associated with the top five most significant GO terms and the fold changes of the significant genes associated with these terms (color). The size of the GO terms reflects the pvalues of the terms, with the more significant terms being larger. This plot is particularly useful for hypothesis generation in identifying genes that may be important to several of the most affected processes. 

```r
## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector

OE_foldchanges <- sigOE$log2FoldChange

names(OE_foldchanges) <- sigOE$gene

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange = OE_foldchanges, 
         vertex.label.font=6)
         
```

**Again, to save the figure,** click on the `Export` button in the RStudio `Plots` tab and `Save as PDF...`. Change the `PDF size` to `24 x 32` to give a figure of appropriate size for the text labels.

<img src="img/mov10oe_cnetplot.png" width="800">

Over-representation analysis is only a single type of functional analysis method. Let's look GeneSet Enrichment Analysis (GSEA), a functional class scoring method. 

![Pathway analysis tools](img/pathway_analysis.png)

## Gene set enrichment analysis (GSEA)

[GSEA](http://software.broadinstitute.org/gsea/index.jsp), most often uses the gene-level statistics or log2 fold changes for all genes from the differential expression results, then look to see whether gene sets for particular biological pathways are enriched among the large positive or negative fold changes. 

<img src="img/gsea_theory.png" width="600">

The hypothesis of GSEA is that although large changes in individual genes can have significant effects on pathways (and will be detected via ORA methods), weaker but coordinated changes in sets of functionally related genes (i.e., pathways) can also have significant effects.  Thus, rather than setting an arbitrary threshold to identify 'significant genes', **all genes are considered** in the analysis. The gene-level statistics from the dataset are aggregated to generate a single pathway-level statistic and statistical significance of each pathway is reported. This type of analysis can be particularly helpful if the differential expression analysis only outputs a small list of significant DE genes. 

### GSEA using clusterProfiler

Using the log2 fold changes obtained from the differential expression analysis for every gene, gene set enrichment analysis and pathway analysis can be performed using clusterProfiler tools.

For gene set or pathway analysis using clusterProfiler, coordinated differential expression over gene sets is tested instead of changes of individual genes. "Gene sets are pre-defined groups of genes, which are functionally related. Commonly used gene sets include those derived from KEGG pathways, Gene Ontology terms, MSigDB, Reactome, or gene groups that share some other functional annotations, etc. Consistent perturbations over such gene sets frequently suggest mechanistic changes" [[1](pathway_tools.pdf)]. 


Extract and name the fold changes:

```r
## Extract the foldchanges
foldchanges <- res_tableOE_tb$log2FoldChange

## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_tableOE_tb$gene
```

Next we need to order the fold changes in decreasing order. To do this we'll use the `sort()` function, which takes a vector as input. This is in contrast to Tidyverse's `arrange()`, which requires a data frame.

```r
## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)

head(foldchanges)
```


We can explore the enrichment of BP Gene Ontology terms using gene set enrichment analysis: 

```r

# GSEA using gene sets associated with BP Gene Ontology terms

gseaGO <- clusterProfiler::gseGO(
  geneList = foldchanges,
  ont = "BP",
  keyType = "SYMBOL",
  eps = 0,
  minGSSize = 20,
  maxGSSize = 300,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = TRUE,
  OrgDb = "org.Hs.eg.db",
  by = "fgsea"
)

gseaGO_results <- gseaGO@result

goplot(gseaGO)

gseaplot2(gseaGO, geneSetID = 1:3)

```


We can also use our homemade GO enrichment analysis. To do this we need to load the library GOenrichment:

```

# Uncomment the following if you haven't yet installed GOenrichment.

# devtools::install_github("gurinina/GOenrichment")

library(GOenrichment)

ls("package:GOenrichment")
```

One of the problems with GO enrichment analysis is that the GO annotations are in constant flux. 

Here we can use the GO annotations in `hGOBP.gmt` (downloaded recently) to run GSEA using the `fgsea` package to run GSEA:

```


fgseaRes <-  fgsea::fgseaSimple(pathways = hGOBP.gmt,stats=foldchanges,nperm=1000,maxSize = 300,minSize = 20)

fgsea <- data.frame(fgseaRes,stringsAsFactors = F)

w = which(fgsea$ES > 0)

fposgsea <- fgsea[w,]

fposgsea <- fposgsea %>% arrange(padj)

```


We are going to compare these results to runing the GO enrichment function `runGORESP`. `runGORESP` uses over-representation analysis to identify enriched GO terms, so we need to define a significance cutoff for the `querySet`.
```
args(runGORESP)

?runGORESP

# we'll define our significance cutoff as 0.58, corresponding to 1.5x change.

# `runGORESP` requires a matrix, so we can turn foldchanges into a matrix using `cbind`:
matx <- cbind(foldchanges,foldchanges)

hresp = GOenrichment::runGORESP(fdrThresh = 0.2,mat=matx,
coln=1,curr_exp = colnames(matx)[1], sig = 0.58,
bp_input = hGOBP.gmt,go_input = NULL,minSetSize = 20,
maxSetSize = 300)

names(hresp$edgeMat)
names(hresp$enrichInfo)
View(hresp$enrichInfo[,c(2,3,4,5,10)])
```

Let's check the overlap between the enriched terms found using `runGORESP` and those found using `fgseaSimple` as they used the same GO term libraries:

```
w = which(fposgsea$padj <= 0.2)

lens <- length(intersect(fposgsea$pathway[w],hresp$enrichInfo$term))

length(w)
dim(hresp$enrichInfo)

percent_overlap <- lens/nrow(hresp$enrichInfo)*100

percent_overlap
```

80%, that's very good, especially because we are using two different GO enrichment methods, over-representation analysis and GSEA. The overlap between these enrichment and the ones using the other GO enrichment tools will be very small because of the differences in the GO annotation libraries.

Now to set up the results for viewing in a network, we use the function `visSetup`, which creates a set of nodes and edges in the network, where nodes are GO terms (node size proportional to FDR score) and edges represent the overlap between GO terms (proportional to edge width). This network analysis is based on [Cytoscape](https://cytoscape.org/), an open source bioinformatics software platform for visualizing molecular interaction networks.

```
vis = visSetup(hresp$enrichInfo,hresp$edgeMat)
names(vis)
```

Now we use runNetwork to view the map: 

```

runNetwork(vis$nodes,vis$edges)
```

[Cytoscape](https://cytoscape.org/) is one of the best visualizations available out of all the GO packages.
There are other gene sets available for GSEA analysis in clusterProfiler (Disease Ontology, Reactome pathways, etc.). In addition, it is possible to supply your own gene set GMT file, such as a GMT for MSigDB from the Broad Institute called c2.

** The C2 subcollection CGP: Chemical and genetic perturbations
Gene sets that represent expression signatures of genetic and chemical perturbations.**

## Other Tools

### GeneMANIA

[GeneMANIA](http://genemania.org/) is a tool for predicting the function of your genes. Rather than looking for enrichment, the query gene set is evaluated in the context of curated functional association data and results are displayed in the form of a network. Association data include protein and genetic interactions, pathways, co-expression, co-localization and protein domain similarity. Genes are represented as the nodes of the network and edges are formed by known association evidence. The query gene set is highlighted and so you can find other genes that are related based on the toplogy in the network. This tool is more useful for smaller gene sets (< 400 genes), as you can see in the figure below our input results in a bit of a hairball that is hard to interpret.

![genemania](img/genemania.png)

> Use the significant gene list generated from the analysis we performed in class as input to GeneMANIA. Using only pathway and coexpression data as evidence, take a look at the network that results. Can you predict anything functionally from this set of genes? 


## Resources for functional analysis

* g:Profiler - http://biit.cs.ut.ee/gprofiler/index.cgi 
* DAVID - http://david.abcc.ncifcrf.gov/tools.jsp 
* clusterProfiler - http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
* GeneMANIA - http://www.genemania.org/
* GenePattern -  http://www.broadinstitute.org/cancer/software/genepattern/ (need to register)
* WebGestalt - http://bioinfo.vanderbilt.edu/webgestalt/ (need to register)
* AmiGO - http://amigo.geneontology.org/amigo
* ReviGO (visualizing GO analysis, input is GO terms) - http://revigo.irb.hr/ 
* WGCNA - http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork
* GSEA - http://software.broadinstitute.org/gsea/index.jsp
* SPIA - https://www.bioconductor.org/packages/release/bioc/html/SPIA.html
* GAGE/Pathview - http://www.bioconductor.org/packages/release/bioc/html/gage.html

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
