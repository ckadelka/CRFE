# Concise Ranked Functional Enrichment (CRFE) User Manual
Welcome to CRFE! This user manual provides a comprehensive guide to using this software effectively.

CRFE provides a high-level interpretation of high-throughput gene expression data. For this type of analysis, the user possesses a ranked list of genes, ranked f.e. by differential expression between treatment and control samples. The genes above a certain threshold are considered differentially expressed or perturbed, all other genes form the background or the unperturbed genes. Given categories that describe the biological functions of different sets of genes (e.g., biological process terms and annotations from the Gene Ontology), we can ask the question of which functional categories best explain the observed gene list. A good collection of functional categories explains as many as possible perturbed genes, and as few as possible unperturbed genes. In particular, particularly many highly perturbed genes should be explained and the collection should consist of only a small number of specific categories. CRFE returns a list of functional categories that satisfies all these requirements. Investigating this list is a good first step to interpreting the results of any high-throughput gene expression experiment. CRFE can also be used in other omics functional enrichment analyses, including genome, epigenome, proteome, and metabolome.

<!--If you use CRFE in your work, please cite the CRFE paper: [insert bioarxiv link].-->

Feel free to link to CRFE in your Methods: https://github.com/ckadelka/CRFE

----

## Download / Installation
A working Python distribution and the program `crfe.py` are required to use CRFE. The program has been developed and tested extensively under Python 3.10 but should work with lower Python 3.x versions as well. It requires the packages, `sys`, `math`, `random`, `os`, `datetime`, `time`, `numpy`, `scipy`, `pandas`, `pickle`, and `optparse`. The correct setup can be tested by executing the example below and comparing the output to the provided output.

----

## Implementation
### Data preprocessing
CRFE requires two input file: an annotation file (describing categories with gene annotations such as Gene Ontology Annotation File or KEGG pathway) and a list of genes to be enriched (containing gene names and their corresponding levels or ranking).

- Annotation file is a tab-separated file with two columns: the first column contains gene categories and the second column contains all gene annotations in that category, whitespace-separated. Below is an example of a valid annotation file:

| | |
|--|--|
| reproduction | A1CF A2M AAAS ABAT ABCC8 ABHD2 ACE |
|  developmental growth | ABL1 ACVR1C ACVR2B ADAM15 ADARB1 ADM |
| | |
- Gene file is a tab-separated file with two columns: the first column contains gene names (using the same naming system as the gene names in the annotation file) and the second column contains numeric values upon which the list of genes will be ranked (in descending order by default). If the values of the genes are not provided, the gene order in the first column will be used.

#### Getting the annotation file from a GAF:
We also provide a Python script to get the required annotation file using a Gene Ontology Annotation File (GAF) as the input. If users wish to parse a GAF to a valid annotation file for CRFE, an OBO file is needed and can be retrieved from [GO Consortium](http://release.geneontology.org/). 

To run the `parse_DAG.py` script, use the following command:
```
python parse_DAG.py -i <path to obo file> -ont <ontology to use> -g <path to GAF> -o <output directory>
```

List of arguments used in `parse_DAG.py`:

```
  -h, --help            show this help message and exit
  --inputobo INPUTOBO, -i INPUTOBO
                        An .obo file of the Gene Ontology, preferably obtained 
                        from the same release as the Gene Annotation file .gaf to avoid obsolete terms.
  --ontology ONTOLOGY, -ont ONTOLOGY
                        Which Gene Ontology namespace to use, one of 'F','P', or 'C' 
                        (corresponds to molecular_function, biological_process, or cellular_component, respectively),
                        default to 'P'
  --gafinput GAFINPUT, -g GAFINPUT
                        A Gene Ontology annotation file (.gaf), preferably obtained from 
                        the same release as the Gene Ontology structure .obo file to avoid obsolete annotations
  --odir ODIR, -o ODIR  Annotation file output directory
```

The script creates a directed acyclic graph using the provided `obo` file in the chosen ontology (Molecular Function, Biological Process, or Cellular Component). Any gene annotated to a Gene Ontology (GO) term is also annotated to the ancestors of that term; therefore, the gene annotations in GAF are propagated using the constructed graph. In addition, GO terms with exactly similar gene annotations are combined. As a result, the script outputs two annotation files in the specified output directory which differ in the format of GO terms, one with GO term ID (e.g., GO:0000003), one with GO term name (e.g., reproduction). The gene annotations (the second column) are similar in these two files. 


### Running CRFE
CRFE accepts different types of command line options. See the [table](#complete-option-list) below or type `crfe.py -h` for details about all available command line options.

A typical invocation of CRFE will look like this:
```
python crfe.py -a <annotations file> -g <gene file> -i test
```

We recommend using CRFE with user-defined parameters based on the data. Specifically, users should set the following parameters:

- `-T THRESHOLD_TYPE` and `-t THRESHOLD`: These parameters should be adjusted based on which genes are considered perturbed.
- `-b BELIEF`: This parameter determines the emphasis CRFE places on explaining the top-ranked genes.
- `-n REPEATS`: The MCMC optimization of CRFE is stochastic. This parameter specifies the number of repeated CRFE runs.

By customizing these parameters, users can optimize the performance of CRFE for their particular use case.


----

## Example
To test CRFE on your local machine, please download the annotation file `GOhuman_ft_named.txt` and the sample data `GSE87340_tumor_normal_log2fc-overexp.txt` from the **data** folder.

To run GSEA with the default settings, use the following command line: 
```
python crfe.py -a GOhuman_propagate_combined.txt -g GSE87340_tumor_normal_log2fc-overexp.txt
```
Once executed, CRFE will create a **saved_data** directory where it will store information about the annotation file and the input gene list for faster reference. The result text file generated by CRFE will be saved in the **output** directory.

----

## Result
The result text file consists of two parts: The first part records all the parameters of the current run for reproducibility. The second part is a table that presents all enriched functional categories (rows, ranked by mean posterior probability), as well as several statistics for each category (columns).

| Statistic      | Description                          |
|------------|--------------------------------------|
| Order    | the rank of the functional category based on its mean posterior probability|
| #(Annotations)     | the number of elements (genes) annotated to the functional category|
| #(Annotations of perturbed genes)    | the number of perturbed elements (genes) annotated to the functional category|
| Mean Posterior Probability     | the average posterior probability across all repeats|
|Standard Deviation Posterior Probability  | The standard deviation of the posterior probability across all repeats|
|p-value|the one-tailed p-value from the hypergeometric distribution (values < 10^{-16} are 0) |
|Average Expression Level| the average expression level (if provided) of all annotated elements (genes)|
|Average position perturbed genes| the average position of all annotated perturbed elements (genes), compared to all perturbed genes|
|Term Name|identifier or name of the category|

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
