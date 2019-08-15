# SlideSeq_implementation_ST
Attempts to implement and adapt the SlideSeq method to identify spatially variable genes to Spatial Transcriptomics data.


Direct_efficient_MHT.R contains the most efficient implementation of the method described by Rodriques et al (Science 363(6434), 
p 1463–1467, 2019. doi 10.1126/science.aaw1219) for identifying spatially variable genes, with Benjamini-Hochberg multiple 
hypothesis testing correction added. 

Direct_slow.R, Direct_slow_2.R and Direct_efficient.R are various attempts towards direct implementation. All are functional, but 
the two first are extremely time-consuming, and the third doesn't include any MHT correction.

Transcript_sampling.R attempts to adapt the method to ST data by sampling a number of transcripts rather than a number of spots, 
but gives enormous numbers of hits, and isn't practically useful.

Thresholding.R attempts to adapt the method to ST data by defining a threshold of a certain number of transcripts that must be
present in order for a gene to be said to be expressed at a spot. This approach is functional, but there seems to be some kind of 
bias, as different thresholds lead to the identification of different spatially variable genes. Unattractive in that it doesn't 
in itself utilize the full range of the ST data. 

Oversampling.R uses conditional probabilities to generate the true distribution for each gene, sampling a larger number of spots 
than those with expression of the gene in question. Attractive since it uses the full range of expression levels, but gives 
enormous numbers of hits, and isn't practically useful.

Oversampling_limited.R combines the direct implementation with a certain extent of oversampling. Gives quite reasonable results, 
but this is not the most useful version.

Oversampling_limited_genes_as_input.R and SlideSeq_pvalues_function.R give p values for genes that are suspected to be spatially 
variable. The former only uses one sample, while the latter is a function capable of integrating the data from several samples.
SlideSeq_pvalues_function.R is useful, but questions regarding its foundation, and thereby how trustworthy its results are, 
remain to be answered.

Oversampling_limited_multiple_samples.R and Oversampling_limited_multiple_samples_attempt_parallel.R analyze multiple samples 
to identify genes with spatial patterns, using one or multiple cores, respectively. The former is quite slow, the latter
faster.

SlideSeq_SVid_function.R and SlideSeq_SVid_function_parallel.R are similar to Oversampling_limited_multiple_samples.R and 
Oversampling_limited_multiple_samples_attempt_parallel.R, but are formulated as functions. NB that multiple samples greatly
increases the number of genes that are identified, something which needs to be investigated further.
