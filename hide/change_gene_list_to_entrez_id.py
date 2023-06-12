import GeneSetRankedEnrichment31t_web as tools
s=0
cutoff=20
top_cutoff=200
belief=5
weight_kind=0
category_kind=0
burnin=100
steps=10000

m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
            0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12.txt',0,
            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
            'output1/', 'DISEASE')
gene_file=m.gene_file

import creates_files_for_funcassociate as trans

(dict_gene,_)=dict_transform(gene_file)

reader=open(m.gene_file,'r')
writer=open(m.gene_file[:-4]+'_ID'+'.txt','w')
c=0
for row in reader:
    row=row.split('\t')
    try:
        writer.write(dict_gene[row[0]]+'\t'+row[1])
    except Exception:
        c+=1
        continue
writer.close()
reader.close()