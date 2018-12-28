# recomb

The *recomb* package simulates and summarizes genetic datasets with recombination rates that change over time.

For details, see:

Cox MP, BR Holland, MC Wilkins and J Schmid. 2013. [Reconstructing past changes in locus-specific recombination rates](https://doi.org/10.1186/1471-2156-14-11). *BMC Genetics* 14:11.

The program *ms_recomb*, a modified version of [Richard Hudson](http://home.uchicago.edu/~rhudson1/)'s [*ms*](http://home.uchicago.edu/%7Erhudson1/source/mksamples.html), uses the coalescent to simulate genetic datasets with recombination rates that change through time. Richard Hudson kindly suggested changes to the code. The resulting code, *ms_recomb*, can simulate changing recombination rates for a single population only.

The program *msstats_recomb*, an extensively modified version of [Kevin Thornton](http://www.molpopgen.org/markdown/krthornt)'s [*msstats*](http://www.molpopgen.org/markdown/software.html), implements *n*-tuple subsampling and calculates a suite of summary statistics with specific focus on recombination events. The code includes a new summary, min(d<sub>ij</sub>), which shows potential for detecting recombination rates through time. The program *msstats_recomb* requires a working installation of the [*libsequence*](https://molpopgen.github.io/libsequence/) library.

## *ms_recomb*

Full instructions on how to use most functions in *ms_recomb* are available in the original [*ms*](http://home.uchicago.edu/%7Erhudson1/source/mksamples.html) documentation.

This sister program, *ms_recomb*, allows you to change recombination rates over time, but only for a single population.  

> Note: Do not use this code to model multiple populations and change recombination at the same time. The changing recombination function is not expected to behave correctly when also simulating multiple populations.

Compilation:
```
gcc -O3 -o ms_recomb ms_recomb.c streec.c rand1.c -lm
```

Worked example:
```
ms_recomb 5 2 -t 6.0 -r 1.0 1000 -eR 0.5 3.0 > ms_recomb.out
```

This command line changes the recombination parameter from 1.0 to 3.0 at time 0.5.  All other parameters are exactly as for *ms*.
   
For more detailed descriptions of these time and recombination rate parameters, see the original [*ms*](http://home.uchicago.edu/%7Erhudson1/source/mksamples.html) documentation.


## *msstats_recomb*

Full instructions on how to use most functions in msstats_recomb are available in the original [*msstats*](http://www.molpopgen.org/markdown/software.html) documentation.

This sister program, *msstats_recomb*, calculates *n*-tuples and specifically focuses on summaries of recombination, including a new summary statistic, min(d<sub>ij</sub>).  Note that this code requires a working installation of the [*libsequence*](https://molpopgen.github.io/libsequence/) library.

Compilation:
```
g++ -O3 -o msstats_recomb msstats_recomb.cc -lsequence
```

Worked examples:<br>

To calculate recombination summaries on an entire *ms* dataset (i.e., with no *n*-tuple subsampling):

```
ms_recomb 5 2 -t 6.0 -r 1.0 1000 -eR 0.5 3.0 | msstats_recomb
```

To calculate recombination summaries on 5 *n*-tuples of size 4:

```
ms_recomb 5 1 -t 6.0 -r 1.0 1000 -eR 0.5 3.0 | msstats_recomb -q 4 -p 5
```

When calculating summaries on *n*-tuples, it is best practice to feed only a single *ms_recomb* dataset to *msstats_recomb*.  If multiple datasets are provided, the first *p* lines contain *n*-tuple permutations of the first dataset, the second *p* lines contain *n*-tuple permutations of the second dataset, and so forth.
