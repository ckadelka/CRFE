#Run Gorilla, check option to get output in xls file,
#open xls file and save as comma delimited csv file, named GO.csv
import graphics_gene_enrichment as gr
import GeneSetRankedEnrichment33h as tools
#import GeneSetRankedEnrichment32d_web as tools
import csv, math

method_name="GSEA"
disease_name='GDS5092'
s=0
cutoff=20
top_cutoff=200
belief=2
weight_kind=0
category_kind=0
burnin=100
steps=10000


identifier = "gds"+disease_name+"_gsea"

for min_level in [30]:
    gf = open('../GSEA-P-R/GDS_data/GSEA_results_'+disease_name+'/gsea_report_for_na_pos.txt','rU') 
           
    reader = gf.read().splitlines()
    gf.close()
    
    C_help=[]
    cc=[]
    i=0
    for row in reader:
        row=row.split('\t')
        try:
            GOid=int(row[0][3:])
        except:
            continue
        C_help.append(GOid)
        i+=1

    try:
        if disease_name=='Liver':
            m = tools.CRFE(1, s, cutoff, top_cutoff, belief, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt',0,
                                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                            'output/', 'LIVER')
        elif disease_name.lower()=='glaucoma':
            m = tools.CRFE(1, s, cutoff, top_cutoff, belief, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'datasets/disease_glaucoma.txt',0,
                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                        'output1/', 'DISEASE')
        elif disease_name.lower()=='obesity':
            m = tools.CRFE(1, s, cutoff, top_cutoff, belief, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'datasets/disease_obesity.txt',0,
                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                        'output1/', 'DISEASE')
        elif disease_name.lower()=='adenocarcinoma_of_lung':
            m = tools.CRFE(1, s, cutoff, top_cutoff, belief, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'datasets/disease_adenocarcinoma_of_lung.txt',0,
                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                        'output1/', 'DISEASE')
        elif disease_name.lower()=='ischemic_cardiomyopathy':
            m = tools.CRFE(1, s, cutoff, top_cutoff, belief, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'datasets/disease_ischemic_cardiomyopathy.txt',0,
                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                        'output1/', 'DISEASE')
        elif disease_name.lower()[:3]=='gds':
            m = tools.CRFE(1, s, cutoff, top_cutoff, belief, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'datasets/gene_file_'+disease_name+'.txt',0,
                        'ONTOLOGY_PARSER/annotation_file_human_association_human_biological_process.txt',
                        'output1/', 'NCBI_'+disease_name)
    except:
        alpha=0.5
        beta=0.25
        prob=0.001
        if disease_name=='Liver':
            m = tools.GeneSetRankedEnrichment(s, cutoff, top_cutoff, alpha, beta, prob, belief, weight_kind, category_kind, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt',0,
                                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                            'output/', 'LIVER')
        elif disease_name.lower()=='glaucoma':
            m = tools.GeneSetRankedEnrichment(s, cutoff, top_cutoff, alpha, beta, prob, belief, weight_kind, category_kind, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'datasets/disease_glaucoma.txt',0,
                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                        'output1/', 'DISEASE')
        elif disease_name.lower()=='obesity':
            m = tools.GeneSetRankedEnrichment(s, cutoff, top_cutoff, alpha, beta, prob, belief, weight_kind, category_kind, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'datasets/disease_obesity.txt',0,
                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                        'output1/', 'DISEASE')
        elif disease_name.lower()=='adenocarcinoma_of_lung':
            m = tools.GeneSetRankedEnrichment(s, cutoff, top_cutoff, alpha, beta, prob, belief, weight_kind, category_kind, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'datasets/disease_adenocarcinoma_of_lung.txt',0,
                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                        'output1/', 'DISEASE')
        elif disease_name.lower()=='ischemic_cardiomyopathy':
            m = tools.GeneSetRankedEnrichment(s, cutoff, top_cutoff, alpha, beta, prob, belief, weight_kind, category_kind, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'datasets/disease_ischemic_cardiomyopathy.txt',0,
                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                        'output1/', 'DISEASE')
        elif disease_name.lower()[:3]=='gds':
            m = tools.GeneSetRankedEnrichment(s, cutoff, top_cutoff, alpha, beta, prob, belief, weight_kind, category_kind, 0.559, 'probs', burnin, steps, 0.5,0.2,20,20,'datasets/gene_file_'+disease_name+'.txt',0,
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
                            
    m.get_sets_and_genes()
    m.load_all()
    #lug=len(m.unique_genes)
    #Tset=set(range(len(m.T)))
    #try:
    #    m.min_level_for_activation=gr.correct_cutoff(min_level/100.,m.unique_genes,m.level_list,m.genes_list,m.gene_file_use)
    #except:
    #    m.min_level_for_activation=gr.correct_cutoff(min_level/100.,m.unique_genes,m.level_list,m.genes_list,m.kind_of_gene_file)
    #    
    #cc=[len(m.T[i]) for i in range(len(m.T))]
    #m.getG()
    #glist=[m.genes_list.index(m.unique_genes[i]) for i in range(len(m.unique_genes))]
    #levels=[m.level_list[glist[i]] for i in range(len(m.unique_genes))]
    #dummy=glist[:]
    #dummy.sort()
    #levels_ordinal = [dummy.index(g) for g in glist]
                
    large_goids_int=map(int,m.large_goids)
    #if large_goids exists
    C=[]
    for c in C_help:
        try:
            C.append(large_goids_int.index(c))
        except:
            pass
    
    #To test whether conversion worked
    for c in C:
        print m.term_names[c]
            
    ##To test whether conversion worked
    #for i,c in enumerate(C[:20]):
    #    print c,term_names[c],cc[c]
    #This list should be similar to the second column of the original csv file

    #a=gr.plot_jaccard_term_by_term2([C],T,[],G,m,min_level,cutoff=20,kind=1,title="FuncAssociate 2.0",identifier="lung_cancer_fa_suppl",SAVE=True,fontsize=17,colors=[gr.rgb('p'),gr.rgb(235,235,235)],STD=False)    
        
    #a=gr.plot_jaccard_term_by_term2a([C],m.T,[],m.G,m,min_level=30,cutoff=10,n_intervals=10,kind=1,colors=[gr.rgb('e')],title='GSEA',identifier=identifier,SAVE=True,fontsize=54,STD=False)
    #a=gr.plot_jaccard_term_by_term2a([C],m.T,[],m.G,m,min_level=30,cutoff=20,n_intervals=10,kind=1,colors=[gr.rgb('e')],title='GSEA',identifier=identifier,SAVE=True,fontsize=54,STD=False)
    #a=gr.plot_jaccard_term_by_term2a([C],m.T,[],m.G,m,min_level=30,cutoff=40,n_intervals=10,kind=1,colors=[gr.rgb('e')],title='GSEA',identifier=identifier,SAVE=True,fontsize=54,STD=False)

    #
    #if HDF:
    #    a=gr.plot_jaccard_term_by_term2([C],T,[],G,m,min_level,cutoff=10,kind=1,title="HDF data, "+method_name+" ("+str(min_level)+"% perturbed genes)",identifier="hdf_gorilla",SAVE=True,fontsize=17)
    #else:
    #    #a=gr.plot_jaccard_term_by_term2([C],T,[],G,m,min_level,cutoff=20,kind=1,title="3-cell vs. 2-cell - GO, "+method_name+" ("+str(min_level)+"% perturbed genes)",identifier="liver_go_fa",SAVE=True,fontsize=17)
    #    a=gr.plot_jaccard_term_by_term2([C],T,[],G,m,min_level,cutoff=20,kind=1,title=method_name,identifier=disease_name+"_go_fa",SAVE=True,fontsize=17)
    
    #stri=gr.create_table_row([C],T,[],levels_ordinal,lug-len(G[-1]),lug,method_name,cutoff=20,kind=1)