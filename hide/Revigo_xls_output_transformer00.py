#Run Gorilla, check option to get output in xls file,
#open xls file and save as comma delimited csv file, named GO.csv
import graphics_gene_enrichment as gr
import GeneSetRankedEnrichment31s_web as tools
import csv, math

method_name="GOrilla"
disease_name="Liver"
#disease_name='Obesity'
disease_name='Lung'
#disease_name='GDS3216'

#disease_name='Ischemic_cardiomyopathy'
HDF=0
s=0
cutoff=20
top_cutoff=200
belief=2
weight_kind=0
category_kind=0
burnin=100
steps=10000

identifier = "liver_revigo" if disease_name=="Liver" else (disease_name[3:]+"_revigo" if 'gds'==disease_name.lower()[:3] else "lung_cancer_revigo")

for min_level in [30]:
    
    try:
        gf = open('../GSEA-P-R/REVIGO/REVIGO_'+disease_name+'.csv','rU')    
        reader = csv.reader(gf, delimiter=',')
        
        C_help=[]
        for row in reader:
            try:
                GOid=int(row[0][3:])
            except:
                continue
            C_help.append(GOid)
        gf.close()
    except IOError:
        gf = open('GSEA-P-R/REVIGO/REVIGO_'+disease_name+'.txt','rU')    
        reader = gf.read().splitlines()
        
        C_help=[]
        for row in reader:
            row=row.split('\t')
            try:
                GOid=int(row[0][3:])
            except:
                continue
            C_help.append(GOid)
        gf.close()        


    if HDF:
        m = tools.GeneSetRankedEnrichment(2, 0, 20, 200, 0.1, 
                                    0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'bkz_gene_list.csv',0,
                                    'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                    'output1/', 'HDF')
    elif disease_name=='Liver':
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
    elif disease_name.lower()=='ischemic_cardiomyopathy':
        m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                    0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'datasets/disease_ischemic_cardiomyopathy.txt',0,
                    'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                    'output1/', 'DISEASE')
    elif disease_name.lower()[:3]=='gds':
        m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                    0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'datasets/gene_file_'+disease_name+'.txt',0,
                    'ONTOLOGY_PARSER/annotation_file_human_association_human_biological_process.txt',
                    'output1/', 'NCBI_'+disease_name)                

    if 'NCBI' in m.out_str:
        if int(m.gene_file.split('/')[-1][13:17]) in [3004,4419,4610,4974]:
            m.go_file='ONTOLOGY_PARSER/annotation_file_human_combined_association_human_biological_process.txt'
        elif int(m.gene_file.split('/')[-1][13:17]) in [2969,3866]:
            m.go_file='ONTOLOGY_PARSER/annotation_file_yeast_combined_association_yeast_biological_process.txt'
        elif int(m.gene_file.split('/')[-1][13:17]) in [1686]:
            m.go_file='ONTOLOGY_PARSER/annotation_file_fly_combined_association_fly_biological_process.txt'
        elif int(m.gene_file.split('/')[-1][13:17]) in [3928]:
            m.go_file='ONTOLOGY_PARSER/annotation_file_rat_combined_association_rat_biological_process.txt'
        elif int(m.gene_file.split('/')[-1][13:17]) in [3246,4423,4929,5092]:
            m.go_file='ONTOLOGY_PARSER/annotation_file_mouse_combined_association_mouse_biological_process.txt'
        elif int(m.gene_file.split('/')[-1][13:17]) in [3216,3933]:
            m.go_file='ONTOLOGY_PARSER/annotation_file_arabidopsis_combined_association_arabidopsis_biological_process.txt'
        #m.out_str+='_'+m.gene_file.split('_')[-1][:-4]
        
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
                
    large_goids_int=map(int,large_goids)
    #if large_goids exists
    C=[]
    for c in C_help:
        try:
            C.append(large_goids_int.index(c))
        except:
            pass
    
    #To test whether conversion worked
    for c in C:
        print term_names[c]

    a=gr.plot_jaccard_term_by_term2a([C],T,[],G,m,min_level=30,cutoff=10,kind=1,colors=[gr.rgb('e')],title="GSEA + Revigo",identifier=identifier,SAVE=True,fontsize=54)
    a=gr.plot_jaccard_term_by_term2a([C],T,[],G,m,min_level=30,cutoff=20,kind=1,colors=[gr.rgb('e')],title="GSEA + Revigo",identifier=identifier,SAVE=True,fontsize=54)
    a=gr.plot_jaccard_term_by_term2a([C],T,[],G,m,min_level=30,cutoff=40,kind=1,colors=[gr.rgb('e')],title="GSEA + Revigo",identifier=identifier,SAVE=True,fontsize=54)

    
    #a=gr.plot_jaccard_term_by_term2a([C],T,[],G,m,min_level=20,cutoff=20,n_intervals=10,kind=1,colors=[gr.rgb('e')],title='GOrilla',identifier='lung_cancer_gorilla',SAVE=True,fontsize=35,STD=False)
    #
    #if HDF:
    #    a=gr.plot_jaccard_term_by_term2([C],T,[],G,m,min_level,cutoff=10,kind=1,title="HDF data, "+method_name+" ("+str(min_level)+"% perturbed genes)",identifier="hdf_gorilla",SAVE=True,fontsize=17)
    #else:
    #    #a=gr.plot_jaccard_term_by_term2([C],T,[],G,m,min_level,cutoff=20,kind=1,title="3-cell vs. 2-cell - GO, "+method_name+" ("+str(min_level)+"% perturbed genes)",identifier="liver_go_fa",SAVE=True,fontsize=17)
    #    a=gr.plot_jaccard_term_by_term2([C],T,[],G,m,min_level,cutoff=20,kind=1,title=method_name,identifier=disease_name+"_go_fa",SAVE=True,fontsize=17)
    
    #stri=gr.create_table_row([C],T,[],levels_ordinal,lug-len(G[-1]),lug,method_name,cutoff=20,kind=1)