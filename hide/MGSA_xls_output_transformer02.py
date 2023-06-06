#Run Gorilla, check option to get output in xls file,
#open xls file and save as comma delimited csv file, named GO.csv
import graphics_gene_enrichment as gr
import GeneSetRankedEnrichment31s_web as tools
import csv, math

method_name="MGSA"
disease_name='Obesity'

disease_name='Cardiomyopathy'
disease_name="Liver"
#disease_name='Lung'
HDF=0
s=0
cutoff=20
top_cutoff=200
belief=2
weight_kind=0
category_kind=0
burnin=100
steps=10000

for min_level in [30]:

    if disease_name=='Liver':
        m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                    0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt',0,
                                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                        'output/', 'LIVER')
    elif disease_name.lower()=='glaucoma':
        m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                    0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'datasets/disease_glaucoma.txt',0,
                    'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                    'output1/', 'DISEASE')
    elif disease_name.lower()=='obesity':
        m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                    0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'datasets/disease_obesity.txt',0,
                    'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                    'output1/', 'DISEASE')
    elif disease_name.lower()=='lung':
        m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                    0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'datasets/disease_adenocarcinoma_of_lung.txt',0,
                    'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                    'output1/', 'DISEASE')
    elif disease_name.lower()=='cardiomyopathy':
        m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                    0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'datasets/disease_ischemic_cardiomyopathy.txt',0,
                    'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                    'output1/', 'DISEASE')                
                    
    (genes_list, level_list) = m.get_sets_and_genes()
    (T, unique_genes, large_goids, term_names) = m.load_all_new(genes_list)
    lug=len(unique_genes)
    Tset=set(range(len(T)))
    m.min_level_for_activation=gr.correct_cutoff(min_level/100.,unique_genes,level_list,genes_list,m.kind_of_list)
    cc=[len(T[i]) for i in range(len(T))]
    (G,_, dict_G)=m.getG(unique_genes,genes_list,T,Tset,level_list)
    glist=[genes_list.index(unique_genes[i]) for i in range(len(unique_genes))]
    levels=[level_list[glist[i]] for i in range(len(unique_genes))]
    dummy=glist[:]
    dummy.sort()
    levels_ordinal = [dummy.index(g) for g in glist]
                
    gf = open('MGSA_files/results/correct_MGSA_'+disease_name+'_Full_results.csv','rU')    
    reader = csv.reader(gf, delimiter=',')
    
    terms=[]
    means=[]
    i=0
    gg=[cc[t]-len(set(T[t])&G[-1]) for t in range(len(T))]
    for row in reader:
        if i==0:
            i+=1
            continue
        try:
            if int(row[1])==cc[int(row[0])] and int(row[2])==gg[int(row[0])]:
                terms.append(int(row[0]))
                means.append(float(row[3]))
            else:
                for j in [-1,1,-2,2,-3,3]:
                    if i+j>0 and i+j<len(T):
                        if int(row[1])==cc[int(row[0])+j] and int(row[2])==gg[int(row[0])+j]:
                            terms.append(int(row[0])+j)
                            means.append(float(row[3]))
                            break
        except:
            continue
        i+=1
