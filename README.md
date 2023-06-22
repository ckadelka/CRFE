# Concise Ranked Functional Enrichment (CRFE) User Manual
Welcome to CRFE! This user manual provides a comprehensive guide to using the software effectively.

CRFE provides a high-level interpretation of high-throughput gene expression data. For this type of analysis, the user possesses a ranked list of genes, ranked f.e. by differential expression between treatment and control samples. The genes above a certain threshold are considered differentially expressed or perturbed, all other genes form the background or the unperturbed genes. Given categories that describe the biological functions of different sets of genes (e.g., biological process terms and annotations from the Gene Ontology), we can ask the question of which functional categories best explain the observed gene list. A good collection of functional categories explains as many as possible perturbed genes, and as few as possible unperturbed genes. In particular, particularly many highly perturbed genes should be explained and the collection should consist of only a small number of specific categories. CRFE returns a list of functional categories that satisfies all these requirements. Investigating this list is a good first step to interpreting the results of any high-throughput gene expression experiment. CRFE can also be used in other omics functional enrichment analyses, including genome, epigenome, proteome, and metabolome.

If you use CRFE in your work, please cite the CRFE paper: [insert bioarxiv link]. 

And feel free to link to CRFE in your Methods: https://github.com/ckadelka/CRFE

----

## Download / Installation
A working Python distribution and the program `crfe.py` are required to use CRFE. The program has been developed and tested extensively under Python 2.7 and Python 3.4 but should work for Python 2.6+. It requires the packages, `sys`, `math`, `random`, `os`, `datetime`, `csv`, `cPickle` (Python 2.x) / `pickle` (Python 3), `optparse`. The correct setup can be tested by executing the example below and comparing the output to the provided output.

----

## Using CRFE
CRFE accepts different types of command line options. See the [table](#complete-option-list) below or type `crfe.py -h` for details about all available command line options.

A typical invocation of CRFE will look like this:
```
python crfe.py -a <annotations file> -g <gene file>
```

We recommend using CRFE with user-defined parameters based on the data. Specifically, user should set `-T THRESHOLD_TYPE` and `-t THRESHOLD` to match the input gene list, `-b BELIEF` for the emphasis of the top-ranked genes, and `-n REPEATS` for the number of repeated CRFE runs.

----

## Example
To test CRFE on your local machine, please download the annotation file `GOhuman_ft_named.txt` and sample data `GSE87340_tumor_normal_log2fc-overexp.txt` from the data folder.

Run GSEA by default using the following command line: 
```
python crfe.py -a 'GOhuman_propagate_combined.txt' -g 'GSE87340_tumor_normal_log2fc-overexp.txt'
```
CRFE will create a `saved_data\` directory that stores information about the annotation file and the input gene list. The result text file CRFE generates is stored in the directory `output\`. 

----

## Result
There are two main parts of the result: The first part records all the parameters for this run. The second part is a table consisting of a ranked set of returned functional categories (row) and their statistics (column).

| Statistic      | Description                          |
|------------|--------------------------------------|
| Order    | The rank of the functional category based on their Mean Posterior Probability|
| #(Annotations)     | The size of the category, indicating the number of elements (genes) annotated to it|
| #(Annotations of perturbed genes)    | The count of perturbed elements (genes) annotated to the category|
| Mean Posterior Probability     | The average posterior probability calculated from all repeats|
|Standard Deviation Posterior Probability  | The standard deviation of posterior probability from all repeats|
|p-value|The p-value obtained from the Hypergeometric test|
|Average Expression Level| The average expression level (second column in the gene list) of elements (genes) in this category|
|Average position perturbed genes| The average position of elements (genes) within the full gene list of the category|
|Term Name|Identifier or name of the category|

----

## Complete option list

```
Options:
  -h, --help            show this help message and exit
  -a ANNOTATIONS_FILE, --annotations_file=ANNOTATIONS_FILE
                        Annotation data (csv or txt). Each row contains the
                        term name followed by a tab, followed by all annotated
                        genes separated by space. No header!
  -g GENE_FILE, --gene_file=GENE_FILE
                        Expression data (csv or txt). Each row contains the
                        gene name, optionally followed by tab plus expression
                        level. No header!
  -t THRESHOLD, --threshold=THRESHOLD
                        If threshold_type==proportion, this describes the
                        proportion of genes considered perturbed, otherwise it
                        describes the threshold in the expression values to be
                        used (default 0.3)
  -T THRESHOLD_TYPE, --threshold_type=THRESHOLD_TYPE
                        Whether threshold should be interpreted as 1. a value,
                        2. proportion ('proportion') of perturbed genes, or 3.
                        different levels of genes already provided in the
                        expression data ('levels'), (default proportion)
  -c LOWER_CUTOFF, --lower_cutoff=LOWER_CUTOFF
                        Only terms with at least this many annotations are
                        considered (default 20)
  -C UPPER_CUTOFF, --upper_cutoff=UPPER_CUTOFF
                        Only terms with at most this many annotations are
                        considered. Use 0 to excluded an upper cutoff (default
                        500).
  -b BELIEF, --belief=BELIEF
                        Belief for how much more active the highest gene set
                        is compared to the least active gene set (default 5)
  -n REPEATS, --repeats=REPEATS
                        Number of independent repeats of CRFE (default 1).
  -s BURNIN_STEPS, --burnin_steps=BURNIN_STEPS
                        Length of (unrecorded) burnin period of MCMC
                        simulation. Used to initialize the Markov chain
                        (default 50000)
  -S MCMC_STEPS, --MCMC_steps=MCMC_STEPS
                        Length of (recorded) steps in the MCMC simulation
                        after the (unrecorded) burnin period (default 50000)
  -x ALPHA_BETA_MAX, --alpha_beta_max=ALPHA_BETA_MAX
                        Maximal value for the false positive and false
                        negative rate that can be learned (default 0.5).
                        Should not be changed. CRFE already automatically
                        chooses the false poositive and negative rates such
                        that (i) an explained perturbed gene is never
                        penalized more than an unexplained perturbed gene, and
                        that (ii) an unexplained unperturbed gene is never
                        penalized more than an explained unperturbed gene.
  -X NUMBER_DIFFERENT_FALSE_RATES, --number_different_false_rates=NUMBER_DIFFERENT_FALSE_RATES
                        Number of different possible values for the false
                        positive and false negative rate (default 20).
  -Q NUMBER_DIFFERENT_PENALIZATIONS, --number_different_penalizations=NUMBER_DIFFERENT_PENALIZATIONS
                        Number of different possible values for the term size
                        penalization parameter q (default 20)
  -p PROBABILITY_PARAMETER_CHANGE, --probability_parameter_change=PROBABILITY_PARAMETER_CHANGE
                        Proportion of time parameters are changed in MCMC
                        (default 0.2), if 0 no parameters are being learnt.
  -r PARAMETER_INITIAL_MCMC_SET, --parameter_initial_MCMC_set=PARAMETER_INITIAL_MCMC_SET
                        When nonzero, the algorithm does not start with an
                        empty set. If r>=1, r random processes are chosen, if
                        0<r<1, each process is in the initial set with
                        probability r (default 0, i.e., empty set)
  -o OUTPUT_FOLDER, --output_folder=OUTPUT_FOLDER
                        Local reference to folder where output is saved
                        (default: output/)
  -i IDENTIFIER, --identifier=IDENTIFIER
                        Identifying string used to distinguish different
                        experiments (empty default)
  -v                    If True, more information is printed during the CRFE
                        run (slows down CRFE).
  -R SEED, --seed=SEED  If nonnegative, this seed is used too initialize the
                        random number generators (default: -1, i.e. random
                        seed)
```

----

## Contact
Please contact Xinglin Jia (xjia@iastate.edu), An Phan (ahphan@iastate.edu), or Claus Kadelka (ckadelka@iastate.edu) with any questions about this program.
