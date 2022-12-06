import sys
import GeneSetRankedEnrichment32e_web as tools
#
gene_file=sys.argv[1]
onto_file=sys.argv[2]
#
#print gene_file+","+onto_file+","+str(4) #for testing purposes only

#gene_file="uploaded_data/bkz_gene_list.txt"
#onto_file='human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv'

m = tools.GeneSetRankedEnrichment(1, 1, 0, 0.1, 0.25, 0.01,
                            5, 0, 0, 0.5, 'proportion', 20000, int(1e5), 
                            0.5,0.2,19,20,
                            gene_file,0,onto_file,'output/','FULL')

m.get_sets_and_genes()

#print str(sum([l<0 for l in level_list]))+',0,-2'

unique_genes=m.loadvar('unique_genes')
if unique_genes==None:
    (_,unique_genes,_,_)=m.load_all_new()
    m.savevar('unique_genes', unique_genes)
    
print str(len(m.genes_list))+","+str(len(unique_genes))+','+str(len(set(m.genes_list)&set(unique_genes)))
