import GeneSetRankedEnrichment31q_web as tools

#leave these parameters fixed
cutoff=20
top_cutoff=200

#really change these parameters
which=3 #0: HDF, 1: 3-cell vs 2-cell - GO, 2: 3-cell vs 2-cell - MSigDB

#code, don't change
if which==0:
    m = tools.GeneSetRankedEnrichment(2, 0, cutoff, top_cutoff, 0.1, 
                                0.25, 1e-50, 5,0, 0.559, 1000, 1000, 0,0, 0.5,0.2,20,20,'bkz_gene_list.csv',0,
                                'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                'output1/', 'HDF')
elif which==1:
    m = tools.GeneSetRankedEnrichment(2, 0, cutoff, top_cutoff, 0.1, 
                            0.25, 1e-50, 5,0, 0.08016865, 1000, 1000, 5,0, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt',0,
                                'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                'output/', 'LIVER')
elif which==2:
    m = tools.GeneSetRankedEnrichment(2, 0, cutoff, top_cutoff, 0.1, 
                            0.25, 1e-50, 5,0, 0.08016865, 1000, 1000, 0,0, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12.txt',0,
                                'msigdb.txt',
                                'output/', 'LIVER')
elif which==3:
    m = tools.GeneSetRankedEnrichment(2, 0, cutoff, top_cutoff, 0.1, 
                0.25, 1e-50, 5,0, 0.559, 1000, 1000, 0,0, 0.5,0.2,20,20,'datasets/disease_adenocarcinoma_of_lung.txt',0,
                'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                'output1/', 'DISEASE')  
                    
(genes_list, level_list) = m.get_sets_and_genes()
(T, unique_genes, large_goids, term_names) = m.load_all_new(genes_list)

if m.out_str=='DISEASE':
    m.out_str=m.gene_file.split('_')[-1][:-4]  

#TRANSFORM AND SAVE TO FILE
filename='ontology_for_funcass'+m.out_str+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)

f_out = open(filename+'.txt', 'w')
for i in range(len(T)):
    text='MaxNala'+str(i+1)+'\t'+term_names[i]+'\t'
    
    for ind in T[i]:
        text+=str(unique_genes[ind])+' '
    if i==16:
        print text
    text=text[:-1]+'\n'
    f_out.write(text)
f_out.flush()
f_out.close()
