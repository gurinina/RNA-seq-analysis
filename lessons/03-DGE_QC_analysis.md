# Quality Control

The next step in the DESeq2 workflow is QC, which includes sample-level and gene-level steps to perform QC checks on the count data to help us ensure that the samples/replicates look good. 

<img src="img/deseq_workflow_qc_2018.png" width="400">

## Sample-level QC

A useful initial step in an RNA-seq analysis is often to assess overall similarity between samples: 

- Which samples are similar to each other, which are different? 
- Does this fit to the expectation from the experiment’s design? 
- What are the major sources of variation in the dataset?

To explore the similarity of our samples, we will be performing sample-level QC using Principal Component Analysis (PCA) and hierarchical clustering methods. Our sample-level QC allows us to see how well our replicates cluster together, as well as, observe whether our experimental condition represents the major source of variation in the data. Performing sample-level QC can also identify any sample outliers, which may need to be explored further to determine whether they need to be removed prior to DE analysis. 

<img src="img/sample_qc.png" width="700">

When using these unsupervised clustering methods, log2-transformation of the normalized counts improves the distances/clustering for visualization. DESeq2 uses a **regularized log transform** (rlog) of the normalized counts for sample-level QC as it moderates the variance across the mean, improving the clustering.

<img src="img/rlog_transformation.png" width="500">

### [Principal Component Analysis (PCA)](https://hbctraining.github.io/DGE_workshop/lessons/principal_component_analysis.html)

Principal Component Analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset (dimensionality reduction). Details regarding PCA are given below but if you would like a more thorough description, you can explore [StatQuest's video](https://www.youtube.com/watch?v=_UVHneBUBW0). 

Suppose we had a dataset with two samples and four genes. Based on this expression data we want to evaluate the relationship between these samples. We could plot the counts of one sample versus another, with Sample 1 on the x-axis and Sample 2 on the y-axis as shown below:

<img src="img/PCA_2sample_genes.png" width="600">

For PCA analysis, the first step is taking this plot and drawing a line through the data in the direction representing the most variation (best fit). In this example, the most variation is along the diagonal. That is, the **largest spread in the data** is between the two endpoints of this line. **This is called the first principal component, or PC1.**  The genes at the endpoints of this line (Gene B and Gene C) have the **greatest influence** on the direction of this line. 

<img src="img/pca_with_PC1_line.png" width="300">

After drawing this line and establishing the amount of influence per gene, **PCA will compute a per sample score**. The per sample PC1 score is computed by taking the product of the influence and the normalized read count and summing across all genes. 

Calculating the influence of each gene is a bit complicated, but to give you an idea,the first step is to calculate a z-score for each gene:

<img src="img/zscore.png" style= "width:150px; align:center"/>

The z-score is a measure of variability. So it's easy to see how in our plot, the influence of the two endpoints will be greater than the other points as they are further away from the mean and their z-scores will be larger. Therefore these points will have a greater impact on PC1.

```
Sample1 PC1 score = (read count Gene A * influence Gene A) + (read count Gene B * influence Gene B) + .. for all genes
```

We could draw another line through the data representing the second most amount of variation in the data (PC2) and compute scores, followed by a third line and so on until you hit the total number of samples in your dataset.

In reality, your dataset will have larger dimensions (more samples, and many, many more genes). The initial sample-to-sample plot, will therefore be in *n*-dimensional space with *n* axes representing the total number of samples you have. The end result is a 2-dimensional matrix with rows representing samples and columns reflecting scores for each of the principal components. To evaluate the results of a PCA, we usually plot principal components against each other, starting with PCs that explain the most amount of variation in your data.

<img src="img/PCA_samples.png" width="600">

**If two samples have similar levels of expression for the genes that contribute significantly to the variation represented by PC1, they will be plotted close together on the PC1 axis.** Therefore, we would expect that biological replicates to have similar scores (since the same genes are changing) and cluster together on PC1 and/or PC2, and the samples from different treatment groups to have different score. This is easiest to understand by visualizing example PCA plots.

#### Interpreting PCA plots

We have an example dataset and a few associated PCA plots below to get a feel for how to interpret them. The metadata for the experiment is displayed below. The main condition of interest is `treatment`.

<img src="img/example_metadata.png" width="600">

When visualizing on PC1 and PC2, we don't see the samples separate by `treatment`, so we decide to explore other sources of variation present in the data. We hope that we have included all possible known sources of variation in our metadata table, and we can use these factors to color the PCA plot. 

<img src="img/example_PCA_treatmentPC1.png" width="600">

We start with the factor `cage`, but the `cage` factor does not seem to explain the variation on PC1 or PC2.

<img src="img/example_PCA_cage.png" width="600">

Then we color by the `strain` factor and find that it explains the variation on PC1. 

<img src="img/example_PCA_treatmentPC3.png" width="600">

It's great that we have been able to identify the sources of variation for both PC1 and PC2. By accounting for it in our model, we should be able to detect more genes differentially expressed due to `treatment`.

***

**Exercise**

The figure below was generated from a time course experiment with sample groups 'Ctrl' and 'Sci' and the following timepoints: 0h, 2h, 8h, and 16h. 

* Determine the sources explaining the variation represented by PC1 and PC2.

* Do the sample groups separate well?

* Do the replicates cluster together for each sample group?

* Are there any outliers in the data?

* Should we have any other concerns regarding the samples in the dataset?

<img src="img/PCA_example3.png" width="600">

***

### Hierarchical Clustering Heatmap

Similar to PCA, hierarchical clustering is another, complementary method for identifying strong patterns in a dataset and potential outliers. The heatmap displays **the correlation of gene expression for all pairwise combinations of samples** in the dataset. Since the majority of genes are not differentially expressed, samples generally have high correlations with each other (values higher than 0.80). Samples below 0.80 may indicate an outlier in your data and/or sample contamination.  

The hierarchical tree can indicate which samples are more similar to each other based on the normalized gene expression values. The color blocks indicate substructure in the data, and you would expect to see your replicates cluster together as a block for each sample group. Additionally, we expect to see samples clustered similar to the groupings observed in a PCA plot. 

**In the plot below, we would be quite concerned about 'Wt_3' and 'KO_3' samples not clustering with the other replicates. We would want to explore the PCA to see if we see the same clustering of samples.**

<img src="img/heatmap_example.png" width="500">


## Gene-level QC

In addition to examining how well the samples/replicates cluster together, there are a few more QC steps. Prior to differential expression analysis it is beneficial to omit genes that have little or no chance of being detected as differentially expressed. The genes omitted fall into three categories:

- Genes with zero counts in all samples
- Genes with an extreme count outlier
- Genes with a low mean normalized counts

<img src="img/gene_filtering.png" width="600">

**DESeq2 will perform this filtering by default**.  Filtering is a necessary step as it will increase the power to detect differentially expressed genes. 

## Mov10 quality assessment and exploratory analysis using DESeq2	

Now that we have a good understanding of the QC steps normally employed for RNA-seq, let's implement them for the Mov10 dataset we are going to be working with.

### Transform normalized counts using the rlog transformation

**To improve the distances/clustering for the PCA and heirarchical clustering visualization methods**, we need to moderate the variance across the mean by applying the rlog transformation to the normalized counts. 

> The rlog transformation of the normalized counts is only necessary for these visualization methods during this quality assessment. We will not be using these tranformed counts downstream.

```r
### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
```
The `blind=TRUE` argument results in a transformation unbiased to sample condition information. When performing quality assessment, it is important to include this option. 

The `rlog` function returns a `DESeqTransform` object. We use this object to plot the PCA and heirarchical clustering figures for quality assessment.

### Principal components analysis (PCA)

DESeq2 has a built-in function for plotting PCA plots, that uses `ggplot2` under the hood. It takes the `rlog` object as an input directly, hence saving us the trouble of extracting the relevant information from it.

The function `plotPCA()` requires two arguments as input: an `rlog` object and the `intgroup` (the column in our metadata that we are interested in). 

```r
### Plot PCA 
plotPCA(rld, intgroup="sampletype")
```

![pca](img/pca_500.png)

**What does this plot tell you about the similarity of samples? Does it fit the expectation from the experimental design?** By default the function uses the *top 500 most variable genes*. You can change this by adding the `ntop` argument and specifying how many genes you want to use to draw the plot.

### Hierarchical Clustering

Since there is no built-in function for heatmaps in DESeq2 we will be using the `pheatmap()` function from the `pheatmap` package. This function requires a matrix/dataframe of numeric values as input, and so the first thing we need to is retrieve that information from the `rld` object:

The rlog is a "DESeqTransform" object and similar to the dds DESeqDataSet object has slotNames containing other data.We can see this by:

```r
slotNames(rld)

To extract the count matrix we use:
```r

rld_mat <- assay(rld)   
```

Then we need to compute the pairwise correlation values for samples. We can do this using the `cor()` function:

```r
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames
```

And now to plot the correlation values as a heatmap:

```r
### Plot heatmap
pheatmap(rld_cor)
```

![heatmap1](img/pheatmap-1.png)

Overall, we observe pretty high correlations across the board ( > 0.999) suggesting no outlying sample(s). Also, similar to the PCA plot you see the samples clustering together by sample group. Together, these plots suggest to us that the data are of good quality and we have the green light to proceed to differential expression analysis.

> NOTE: The `pheatmap` function has a number of different arguments that we can alter from default values to enhance the aesthetics of the plot. If you are curious try running the code below. *How does your plot change?* Take a look through the help pages (`?pheatmap`) and identify what each of the added arguments is contributing to the plot.
>
> ```r
> heat.colors <- brewer.pal(6, "Blues")
> pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
>			fontsize_row = 10, height=20)
> ```       
> Curious on all of the available [color palettes offered by the RColorBrewer package](http://www.r-graph-gallery.com/38-rcolorbrewers-palettes/)? Try typing in your console `display.brewer.all()` and see what happens!
> 

***

