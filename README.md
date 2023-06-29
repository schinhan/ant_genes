# Molecular Signatures of Attention Networks

Attention network theory proposes three distinct types of attention - alerting, orienting, and control - that are supported by 
separate brain networks and modulated by different neurotransmitters, i.e., noradrenaline, acetylcholine, and dopamine.  This repository contains code, data, and resources related to our study exploring the extent of cortical, genetic, and molecular 
dissociation of these three attention systems using multimodal neuroimaging. We evaluated the spatial overlap between fMRI activation maps from the attention network test (ANT) and cortex-wide 
gene expression data from the Allen Human Brain Atlas. The goal was to identify genes associated with each of the attention networks 
in order to determine whether specific groups of genes are co-expressed with the corresponding attention networks. Furthermore, 
we analysed publicly available PET-maps of neurotransmitter receptors and transporters to investigate their spatial overlap with 
the attention networks. 

## Code
### Parcellation

Brain Parcellation: Creates the Lausanne-219 parcellation with 219 cortical brain regions. We created our own group version of 
this parcellation based on T1-weighted structural images, aligning individual parcellations with a group atlas for subsequent 
analyses.

Null Model: Employs a spatial permutation approach that randomly swaps parcellation labels while accounting for the 
intrinsic geometry of the cortex.

### fMRI

Overlap of activation maps: Quantifies the degree of spatial overlap between the phenotypic activation maps (brain images)
by calculating product-moment correlations between activation vectors of different attention networks. To test signficance, 
1 million random rotations of each activation vector are created - these are then  correlated with the obsvered activation
vectors of the other conditions.

Gene-expression-similarities of attention networks: correlates each of the three regions-by-activation vectors from the 
ANT conditions with the regions-by-genes expression matrix obtained from abagen. Corresponding p-values are computed by 
randomly rotating the ANT activation maps relative to the gene expression maps.

Gene Set Enrichment Analysis: 
Genes are ranked according to their z-values - these are based on the previously determined p-values and the empirical 
correlations’ signs. No thresholds are applied, and all genes, along with their respective rank indices, 
are preserved without exclusions to perform the gene set enrichment analysis with Panther.

### PET

Relationships between attention networks and 19 PET publicly available maps are assessed through Pearson correlations. 
Corresponding p-values are obtained through permutation testing with the spin test (5000 permutations each).

## Data

### Templates
Contains various templates from the Human Connectome Project which are necessary to execute all code related to the Lausanne parcellation and the creation of the null model.

### Parcellation
Contains a group version of the Lausanne Atlas with all related files and the genes-by-regions expression matrix with a total of 15632 incorporated genes that was created with abagen.


### fMRI
Contains the parcellated ANT images, the observed and permutated associations between activation vectors of different attention networks, the observed and permutated associations between activation vectors and gene expression data, the input and output files of needed for the gene set enrichment analysis.

### PET
Contains all files related to the analysis of PET-maps representing neurotransmitter receptors and transporters to evaluate their spatial overlap with the attention networks.
## Data Visualization
Contains all the code that was used to generate the figures in the study.

### GoFigure
Summarized visualisation of gene set enrichment analysis results. The GO-terms of all 
significant Molecular Functions (p<.05, FDR corrected) with the corresponding p-value serve as input. 

### Heatmaps
Heatmaps show how attention network images are decomposed into activation vectors using Lausanne parcellation. 
Activation vectors and abagen gene expression matrix are then correlated with each other; correlation coefficients 
are displayed.

### Manhattan_QQ_Plots
Manhattan Plots show association results of all expressed genes and the different ANT conditions. Note that Manhattan 
plots do not highlight specific genomic regions that contain accumulations of strong association p-values as typically 
shown in genome-wide association studies. In this expression-based analysis, clusters of genes may span the whole genome 
and still overlap in their cortex-wide expression patterns, producing ‘horizontal band of associations’ with similar 
test-statistics in Manhattan plots.

QQ-plots also refer to association results of all expressed genes and the different ANT conditions. They compare the 
distribution of observed p-values (y-axes) against the distribution of expected p-values under the null hypothesis.

### ViolinePlots
Violine plots of correlations between brain activation maps show the observed correlations of ANT maps vs. the 
distribution of expected correlations after one-million bidirectional random rotations of ANT maps.
Violine plots of correlations between gene associations results showing the observed correlations of gene association
results vs. the distribution of expected correlations after one-million random rotations of ANT maps.

