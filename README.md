# CSKS
Coalescent simulations the generate k-mer spectra

## The motivation
Simulating genetic data is not trivial. In the case of one diploid (i.e. two genomes), one can just throw in heterozygous sites according to the desired heterozygosity level (θ). If more than two genomes are simulated, the task becomes more complicated. This is because related genomes have a common ancestor and thus share mutations. A summary of the number of times different mutations are shared between genomes in a sample is called the [site-frequency spectrum](https://en.wikipedia.org/wiki/Allele_frequency_spectrum).
