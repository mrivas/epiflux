# Data 
Contents:
- epigenetics: Files in .csv format that contain the epigenetic information (ChIP-seq) of the H3K9ac and H3K4me3 marks. This information was obtained in the study by Kuang *et al.* (2014) under the access code [GSE52339](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52339)
- gems: Metabolic network of S. cerevisiae (iMM904) that contains acetylation, methylation and glycogen reactions. This metabolic network was downloaded from the [BiGG Models](http://bigg.ucsd.edu/) repository.
- knownFluxes: Known fluxes oxygen which was obtained by the percentage of dissolved oxygen reported by Kuang *et al.* (2014)
- mediums: Medium proposed by Kuang *et al.* (2014) where the sugars trehalose and glycogen were added.
- transcriptomes: Genetic expression (RNA-seq) across the 16 times of YMC. This information was obtained in the study by Kuang *et al.* (2014) under the access code [GSE52339](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52339)
- InputData_yeast_example.csv : Required input file in ../../algorithm/example.sh

*The RNA-seq and ChIP-seq data changed from 16 times to 15, this is due to the specifications of the study by Sánchez *et al.* (2018) where times 10 and 11 are averaged for RNA-seq and times 13 and 14 for ChIP-seq data.

# References

Kuang, Z., Cai, L., Zhang, X., Ji, H., Tu, B.P. and Boeke, J.D., 2014. High-temporal-resolution view of transcription and chromatin states across distinct metabolic states in budding yeast. *Nature structural & molecular biology*, 21(10), pp.854-863.

King, Z. A., Lu, J. S., Dräger, A., Miller, P., Federowicz, S., Lerman, J. A., Ebrahim, A., Palsson, B. Ø., & Lewis, N. E. (2015). BIGG Models: a platform for integrating, standardizing and sharing genome-scale models. *Nucleic Acids Research*, 44(D1), D515-D522. https://doi.org/10.1093/nar/gkv1049

Sánchez-Gaya, V., Casaní-Galdón, S., Ugidos, M., Kuang, Z., Mellor, J., Conesa, A., & Tarazona, S. (2018). Elucidating the role of chromatin state and transcription factors on the regulation of the yeast metabolic Cycle: a Multi-OmiC Integrative Approach. *Frontiers in Genetics*, 9. https://doi.org/10.3389/fgene.2018.00578
