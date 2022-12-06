import GeneSetRankedEnrichment32d_web as tools
import numpy as np
import cPickle as pickle

nsim=20

cutoff=20
top_cutoff=200
param_learning=1
MCMC=1
weight_kind=1
category_kind=0
burnin=int(1*1e5)
steps=int(1*1e6)
alpha_beta_top=0.5
min_level=30

gene_files=['bkz_gene_list.csv','3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12.txt','datasets/disease_glaucoma.txt','datasets/disease_obesity.txt','3-cell-hep-d12_versus_CS-d12_ID.txt','datasets/disease_adenocarcinoma_of_lung.txt','datasets/disease_ischemic_cardiomyopathy.txt','3-cell-hep-d12_versus_2-cell-hep-d12_pos_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12_neg_ID.txt']
out_strs=['HDF','LIVER','LIVER','DISEASE','DISEASE','LIVER','DISEASE','DISEASE','Lpos','Lneg']
for hdhsfs in [1]:
    for (belief,s) in [(10,0),(5,0),(2,0),(1,1)]:
        for alpha_beta_top in [0.5]:
            m = tools.GeneSetRankedEnrichment(s, cutoff, top_cutoff, 0.1, 
                                            0.25, 1e-50, belief, weight_kind, category_kind, min_level/100., 'percentage', burnin, steps, alpha_beta_top,0.2,20,20,gene_files[hdhsfs],0,
                                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                            'output1/', out_strs[hdhsfs])
            if hdhsfs==2:
                m.go_file='msigdb.txt'
            
            #data collection
            C=[[]]*nsim
            MCMC_distr=[[]]*nsim
            for i in range(nsim):
                print i
                (C[i],MCMC_distr[i])=m.runMe(0,0,param_learning,MCMC,0,0)
            
            if m.out_str=='DISEASE':
                m.out_str=m.gene_file.split('.')[0].split('datasets/disease_')[1]
        #save data
            if alpha_beta_top==0.5:
                addon=""
            else:
                addon='_alphabetatop'+str(int(round(m.alpha_beta_top*1000)))
            if m.go_file=='msigdb.txt':
                f=open('saved_data/multiple_runs_new3_msigdb_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','w+')
            else:
                f=open('saved_data/multiple_runs_new3_go_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','w+')
            pickle.dump([C,MCMC_distr], f)
            f.close()