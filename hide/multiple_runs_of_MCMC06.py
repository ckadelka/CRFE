import GeneSetRankedEnrichment32b_web as tools
import cPickle as pickle
import graphics_gene_enrichment as gr

nsim=20

cutoff=20
top_cutoff=200
param_learning=1
MCMC=1
belief=2
burnin=int(1*1e5)
steps=int(1*1e6)
weight_kind=1
category_kind=0
data1=[[[[0 for iiii in range(4)] for iii in range(4)] for ii in range(3)] for i in range(10)]
Cs=[[[[0 for iiii in range(4)] for iii in range(4)] for ii in range(3)] for i in range(10)]
MCMC_distrs=[[[[0 for iiii in range(4)] for iii in range(4)] for ii in range(3)] for i in range(10)]
beliefs=[2,5,100]
alpha_beta_tops=[0.1,0.2,0.3,0.5]
min_levels=[10,20,30]

ss=[1]
which=[6]
MCMC_threshold=0.1
ymax=1

gene_files=['bkz_gene_list.csv','3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12.txt','datasets/disease_glaucoma.txt','datasets/disease_obesity.txt','3-cell-hep-d12_versus_CS-d12_ID.txt','datasets/disease_adenocarcinoma_of_lung.txt','datasets/disease_ischemic_cardiomyopathy.txt','3-cell-hep-d12_versus_2-cell-hep-d12_pos_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12_neg_ID.txt']
out_strs=['HDF','LIVER','LIVER','DISEASE','DISEASE','LIVER','DISEASE','DISEASE','Lpos','Lneg']

min_level=30
alpha_beta_top=0.5

for hdhsfs in which:
    for s,belief in zip([0,0,0,1],[10,5,2,1]): #s
        count=0
        for ind in [0]: #min_level
            for j in [0]: #belief
                #if s==0 and weight_kind==0 or weight_kind==0 and s==1 and hdhsfs==0:
                #    continue 
                #if jjj==0 and j>1 or jjj==1 and j<2:
                #    continue
                #belief=beliefs[ind]
                #alpha_beta_top=alpha_beta_tops[j]
                #min_level=min_levels[2]#ind]
                #s=ss[jjj]
                if s==1:
                    belief=100
                
                m = tools.GeneSetRankedEnrichment(s, cutoff, top_cutoff, 0.1, 
                                0.25, 1e-50, belief, weight_kind, category_kind, 0.559, 'percentage',burnin, steps, alpha_beta_top,0.2,20,20,gene_files[hdhsfs],0,
                                'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                'output1/', out_strs[hdhsfs])
                if which==2:
                    m.go_file='msigdb.txt'    
                    
                m.get_sets_and_genes()                
                m.load_all_new()                
                m.glist=[m.genes_list.index(m.unique_genes[i]) for i in xrange(len(m.unique_genes))]
                m.levels=[m.level_list[m.glist[i]] for i in xrange(len(m.unique_genes))]                
                if m.threshold_type=='percentage':
                    m.min_level_for_activation=m.find_min_level_of_activation(m.min_level)
                else:
                    m.min_level_for_activation=m.min_level                     
                dummy=m.glist[:]
                dummy.sort()
                m.levels_ordinal = [dummy.index(g) for g in m.glist]
                
                m.Tset=set(range(len(m.T)))
                m.number_of_genes = len(m.unique_genes)                
                m.getG()
                    
                #m.load_all_new()                                                                                                                                                                
                #m.get_sets_and_genes()
                ##(T, unique_genes, large_goids, term_names) = m.load_all_new(genes_list)
                #lug=len(m.unique_genes)
                #m.min_level_for_activation=gr.correct_cutoff(min_level/100.,m.unique_genes,m.level_list,m.genes_list,m.kind_of_gene_file)
                #Tset=set(range(len(m.T)))
                #cc=[len(m.T[i]) for i in range(len(m.T))]
                #(G,_, dict_G)=m.getG()
                #m.glist=[m.genes_list.index(m.unique_genes[i]) for i in range(len(m.unique_genes))]
                #m.levels=[m.level_list[m.glist[i]] for i in range(len(m.unique_genes))]
                #dummy=m.glist[:]
                #dummy.sort()
                #m.levels_ordinal = [dummy.index(g) for g in m.glist]
    
                if m.out_str=='DISEASE':
                    m.out_str=m.gene_file.split('.')[0].split('datasets/disease_')[1]           
                                        
                #load data
                if alpha_beta_top==0.5:
                    addon=""
                else:
                    addon='_alphabetatop'+str(int(round(m.alpha_beta_top*1000)))
                try:
                    nsim=20
                    if m.go_file=='msigdb.txt': #+'_s'+str(m.nr_categories)
                        f=open('saved_data/multiple_runs_new3_msigdb_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
                    else:
                        f=open('saved_data/multiple_runs_new3_go_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
                except IOError:
                    try:
                        nsim=20
                        if m.go_file=='msigdb.txt':
                            f=open('saved_data/multiple_runs_new3_msigdb_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
                        else:
                            f=open('saved_data/multiple_runs_new3_go_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
                    except IOError:
                        nsim=20
                        if m.go_file=='msigdb.txt':
                            f=open('saved_data/multiple_runs_new3_msigdb_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
                        else:
                            f=open('saved_data/multiple_runs_new3_go_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
                            
                [C,MCMC_distr]=pickle.load(f)
                f.close()
                
                
                if hdhsfs==0:
                    title='HDF data'
                    identifier='hdf'
                elif hdhsfs==1:
                    title='Liver'
                    identifier='liver_go'
                elif hdhsfs==2:
                    title='3-cell vs. 2-cell - MSigDB'
                    identifier='liver_msigdb'
                elif hdhsfs in [3,4] or hdhsfs>5:
                    title=m.out_str.capitalize().replace('_',' ')+' - GO'
                    identifier=m.out_str+'_go'
                elif hdhsfs==5:
                    title='3-cell vs. CS - GO'
                    identifier='liverCS_go'
                elif hdhsfs==6:
                    title='Adenocarcinoma of Lung'
                    identifier='lung_cancer'
                #if s==0:
                #    title+=", MCMC"
                #elif s==1:
                #    title+=", GOing Bayesian"
                #title+=" ("+str(min_level)+"% perturbed genes)"
                #identifier+='_s'+str(s)+'_pert'+str(min_level)
                
                #a=gr.plot_jaccard_term_by_term2(C,T,MCMC_distr,G,m,min_level=min_level,cutoff=20,n_intervals=10,kind=1,colors=[[1,0.7,0.7],[0.7,1,0.7]],title=title+", "+str(lug-len(G[-1]))+" categories ("+str(min_level)+"% perturbed genes)" if s==0 else title+", "+str(s)+" categories ("+str(min_level)+"% perturbed genes)",identifier=identifier+'_s'+str(s),SAVE=1)
                #normal:
                #a=gr.plot_jaccard_term_by_term2(C,T,MCMC_distr,G,m,min_level=min_level,cutoff=20,n_intervals=10,kind=1,colors=[[1,0.7,0.7],[0.7,1,0.7]],title=title+", CRFE ("+str(min_level)+"% perturbed genes)" if s==0 else title+", MGSA ("+str(min_level)+"% perturbed genes)",identifier=identifier+'_s'+str(s),SAVE=1)
                
                #for poster:
                #a=gr.plot_jaccard_term_by_term2(C,m.T,MCMC_distr,m.G,m,min_level=min_level,cutoff=20,n_intervals=10,kind=1,colors=[[1,0.7,0.7],[0.7,1,0.7]],title="CRFE" if s==0 else "MGSA",identifier=identifier+'_s'+str(s),SAVE=1,fontsize=17)
                
                if s==0:
                    nr_perturbed=m.number_of_genes-len(m.G[-1])
                
                #print "vsdfgdfgdfg",[len(T[t]) for t in C[0]]
                
                #print hdhsfs,s,min_level,belief
                data1[hdhsfs][0][ind][j]=gr.categories(C,m.T,MCMC_distr,m.levels_ordinal,m.number_of_genes-len(m.G[-1]),m.number_of_genes,10,cutoff=MCMC_threshold,kind=3)
                #gr.plotte_categories_explained2([data1[hdhsfs][0][count][j]],m,len(G[-1]),min_level,legend=[r"CRFE"],title="",identifier=identifier)#,MCMC_threshold=2)

                if j==3:#ss==[0,1] and s==1 and j==1:
                    m.belief=beliefs[ind] #set belief value back to normal for correct filename
                    #print "a"
                    #gr.plotte_categories_explained2([data1[hdhsfs][0][count][j],data1[hdhsfs][1][count][j],data1[hdhsfs][2][count][j]],m,len(G[-1]),min_level,legend=[r"%i categories ($\tau=%s$)" % (nr_perturbed,str(MCMC_threshold)), r"1 category ($\tau=%s$)" % str(MCMC_threshold), r"10 categories ($\tau=%s$)" % str(MCMC_threshold)],title=title+" ("+str(min_level)+"% perturbed genes)",identifier=identifier,MCMC_threshold=MCMC_threshold)
                    #gr.plotte_categories_explained2([data1[hdhsfs][0][count][j],data1[hdhsfs][1][count][j]],m,len(G[-1]),min_level,legend=[r"CRFE ($\tau=%s$)" % str(MCMC_threshold), r"MGSA ($\tau=%s$)" % str(MCMC_threshold)],title=title+" ("+str(min_level)+"% perturbed genes)",identifier=identifier,MCMC_threshold=MCMC_threshold)
                    
                    #gr.plotte_categories_explained2([data1[hdhsfs][0][count][j],data1[hdhsfs][1][count][j]],m,len(G[-1]),min_level,legend=[r"CRFE", r"MGSA"],title=title,identifier=identifier,ymax=ymax,MCMC_threshold=MCMC_threshold)
                    #gr.plotte_categories_explained2([data1[hdhsfs][0][1][j],data1[hdhsfs][0][0][0],data1[hdhsfs][1][count][j]],m,len(G[-1]),min_level,legend=[r"CRFE, belief=5", r"CRFE, belief=2", r"MGSA"],title=title,identifier=identifier,ymax=ymax,MCMC_threshold=MCMC_threshold)
                    
                    gr.plotte_categories_explained2([data1[hdhsfs][0][ind][0],data1[hdhsfs][0][ind][1],data1[hdhsfs][0][ind][2],data1[hdhsfs][0][ind][3]],m,len(G[-1]),min_level,legend=[r"CRFE, 0.1", r"CRFE, 0.2", r"CRFE, 0.3",r"CRFE, 0.5"],title=title,identifier=identifier,ymax=ymax,MCMC_threshold=MCMC_threshold,colors=['k','b','r','g'])

                count+=1