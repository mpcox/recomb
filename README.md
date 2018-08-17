# recomb

## The recomb package simulates and summarizes genetic datasets with recombination rates that change over time.

The program *ms_recomb*, a modified version of [Richard Hudson](http://home.uchicago.edu/~rhudson1/)'s [*ms*](http://home.uchicago.edu/%7Erhudson1/source/mksamples.html), uses the coalescent to simulate genetic datasets with recombination rates that change through time. Richard Hudson kindly suggested changes to the code.

The program *msstats_recomb*, an extensively modified version of [Kevin Thornton](http://www.molpopgen.org/markdown/krthornt)'s [*msstats*](http://www.molpopgen.org/markdown/software.html), implements *n*-tuple subsampling and calculates a suite of summary statistics with specific focus on recombination events. The code includes a new summary, min(d<sub>ij</sub>), which shows some potential for detecting recombination rates through time. *msstats_recomb* requires a working installation of the [*libsequence*](https://molpopgen.github.io/libsequence/) library.

For details, see:

Cox MP, BR Holland, MC Wilkins, J Schmid. 2013. [Reconstructing past changes in locus-specific recombination rates](https://doi.org/10.1186/1471-2156-14-11). *BMC Genetics* 14:11.

Author: [Murray Cox](https://www.genomicus.com)
