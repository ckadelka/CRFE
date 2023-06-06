import GeneSetRankedEnrichment31b as tools

cutoff=10
top_cutoff=500

m = tools.GeneSetRankedEnrichment(2, 1, 1, 0, 0.05, 
                                0.25, 1e-80, 5, 0.559889, 1000, 1000, 20000, int(1e6), 19,20,'bkz_gene_list.csv',
                                'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                'output/', 'HDF')
                                
unique_genes=m.loadvar('unique_genes')
lug=len(unique_genes)
T=m.loadvar('T')
large_goids=m.loadvar('large_goids')
Tset=set(range(len(T)))
term_names=m.loadvar('term_names')
cc=[len(T[i]) for i in xrange(len(T))]

#remove terms that have too many or too little annotations
for i in xrange(len(T)-1,-1,-1):
    if cc[i]>top_cutoff or cc[i]<cutoff:
        cc.pop(i)
        T.pop(i)
        large_goids.pop(i)
        term_names.pop(i)
        
#remove genes that are no longer annotated by any term
(unique_genes,T) = m.remove_genes_at_most_annotated_to_root(unique_genes,T,{},kind=0)