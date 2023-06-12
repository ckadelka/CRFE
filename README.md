# Concise Ranked Functional Enrichment (CRFE) User Manual
CRFE is a functional enrichment method that takes a full list of genes and functional category annotations as input and outputs a small set of non-redundant functional categories that are significantly over-represented in highly ranked genes. It can also be used in other omics functional enrichment analyses, including genome, epigenome, proteome, and metabolome.

If you use CRFE in your work, please cite the CRFE paper: [insert bioarxiv link]. 
And feel free to link to CRFE in your Methods: [insert git link]

## Requirements
  <list package versions>

## Workflow with sample data
Please download the main python script <crfe.py>, the annotation file <xxx> and sample data <xxx>

'''
crfe.py -annotation_file 'data/GOhuman_ft_named.txt' -gene_file 'data/GSE87340_tumor_normal_log2fc-overexp.txt' -identifier 'real_an' -threshold
'''
## Result

## Complete option list
'''
usage: crfe [-h]
            [-x]
            [-x]
            [-x]
            [-x]
optional arguments:
  -h, --help show this help message and exit\\
  -x  --gene_file\\
              direction of the input gene list\\
  -x  --annotations_file\\
              direction of the annotation file\\
  -x  --repeats <10>\\
              number of repeats\\
              [DEFAULT: 10]
  -x  --nr_categories <0>
              number of categories in result
              [DEFAULT: 0]
  -x  --lower_cutoff <20>
              minimal number of genes in category
              [DEFAULT: 20]
  -x  --upper_cutoff <500>
              maximual number of genes in category
              [DEFAULT: 500]
  -x  --belief <5>
              belief parameter, the emphysiasis ratio of the top ranked gene (most perturbed) and the least bottom ranked gene (least perturbed) in the P set
              [DEFAULT: 5]
  -x  --threshold_type {levels, proportion}
              set the threshold type for input gene list. "levels" requires a second column with numeric values in the input gene list 
  -x  --burnin_steps <50000>
              set the burin steps in MCMC process
              [DEFAULT: 50000]
  -x  --MCMC_steps <50000>
              set the sampling steps in MCMC process
              [DEFAULT: 50000]
  -x  --alpha_beta_max <1>
              set the max value of alpha and beta
              [DEFAULT: 1]
  -x  --proportion_parameter_change <0.2>
              set the proportion of making a proposal of paramters alpha, beta, and penalization parameter in the MCMC process
              [DEFAULT: 0.2]
  -x  --A
  -x  --P
  -x  --output_folder
  -x  --out_str
  -x  --LEARN_PENALIZATION_PARAMETER
  -x  --penalization_parameter
'''

## Contributions
