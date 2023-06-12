import GeneSetRankedEnrichment31q_web as tools
import pylab as pl
import numpy as np

cutoff=5
top_cutoff=100
hist_steps=top_cutoff-cutoff+1

m = tools.GeneSetRankedEnrichment(2, 0, cutoff, top_cutoff, 0.05, 
                                   0.25, 0.001, 2, 0, 0.05016865, 10, 25, 2000, 10000, 0.05157523,0.2,19,20, '3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt',0,#'3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt',
                                   'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                   'output/', 'LIVER')

(genes_list, level_list) = m.get_sets_and_genes()
(T, unique_genes, large_goids, term_names) = m.load_all_new(genes_list)
lug=len(unique_genes)
cc=[len(T[i]) for i in xrange(len(T))]
#a=pl.hist(cc,hist_steps)
#pl.title('Occurence of GO terms that annotate a given number of genes')


#gg=[0 for g in unique_genes]
#for t in T:
#    for gene in t:
#        gg[gene]+=1
#hist_steps=max(gg)+1
#a=pl.hist(gg,hist_steps/5)
#pl.title('Occurence of genes that are annotated by a given number of GO terms')

glist=[genes_list.index(unique_genes[i]) for i in xrange(len(unique_genes))]
levels=[level_list[glist[i]] for i in xrange(len(unique_genes))]
mm=m.mean_gene_level(range(len(T)),[],T,levels,[])
a=pl.hist(mm,20)
pl.title('Histogram of the average expression level of all terms')

#h,x=tools.hist(cc,hist_steps,0,top_cutoff)
#pl.plot([el+0.5/hist_steps*top_cutoff for el in x],h)



pl.show()