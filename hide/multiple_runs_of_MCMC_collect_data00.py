import GeneSetRankedEnrichment31t_web as tools
import numpy as np
import cPickle as pickle

def correct_cutoff(min_level,unique_genes,level_list,genes_list,kind_of_file):
    nr_active=int(round(len(unique_genes)*float(min_level)))
    counter=0
    BROKE=0
    for i in xrange(len(genes_list)):
        if genes_list[i] in unique_genes:
            counter+=1
            if counter>=nr_active:
                broke_at=i
                BROKE=1
                break
    if BROKE==0:
        broke_at=len(genes_list)-1
        nr_active=len(unique_genes)
    counter=len(set(unique_genes)&set(genes_list[:counter]))
    help=level_list[broke_at]
    if int(kind_of_file)==1:
        return -help
    else:
        return help

nsim=20

cutoff=20
top_cutoff=200
param_learning=1
MCMC=1
weight_kind=1
category_kind=0
burnin=int(1*1e3)
steps=int(1*1e4)
s=0
gene_files=['bkz_gene_list.csv','3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12.txt','datasets/disease_glaucoma.txt','datasets/disease_obesity.txt','3-cell-hep-d12_versus_CS-d12_ID.txt','datasets/disease_adenocarcinoma_of_lung.txt','datasets/disease_ischemic_cardiomyopathy.txt']
out_strs=['HDF','LIVER','LIVER','DISEASE','DISEASE','LIVER','DISEASE','DISEASE']
for hdhsfs in [6]:
    for belief in [100]:
        for min_level in [30]:
            m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                            0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,gene_files[hdhsfs],0,
                                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                            'output1/', out_strs[hdhsfs])
            if hdhsfs==2:
                m.go_file='msigdb.txt'
                                                                                        
            (genes_list, level_list) = m.get_sets_and_genes()
            (T, unique_genes, large_goids, term_names) = m.load_all_new(genes_list)
            m.min_level_for_activation=correct_cutoff(min_level/100.,unique_genes,level_list,genes_list,m.kind_of_list)
            lug=len(unique_genes)
            Tset=set(range(len(T)))
            cc=[len(T[i]) for i in range(len(T))]
            (G,_, dict_G)=m.getG(unique_genes,genes_list,T,Tset,level_list)
            glist=[genes_list.index(unique_genes[i]) for i in range(len(unique_genes))]
            levels=[level_list[glist[i]] for i in range(len(unique_genes))]
            dummy=glist[:]
            dummy.sort()
            levels_ordinal = [dummy.index(g) for g in glist]
            
            #data collection
            C=[[]]*nsim
            MCMC_distr=[[]]*nsim
            for i in range(nsim):
                print i
                (C[i],MCMC_distr[i])=m.runMe(0,0,param_learning,MCMC,0,0)
            
            if m.out_str=='DISEASE':
                m.out_str=m.gene_file.split('.')[0].split('datasets/disease_')[1]
        #save data
            if m.go_file=='msigdb.txt':
                f=open('saved_data/multiple_runs_new_msigdb_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','w+')
            else:
                f=open('saved_data/multiple_runs_new_go_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','w+')
            pickle.dump([C,MCMC_distr], f)
            f.close()