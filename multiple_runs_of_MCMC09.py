import GeneSetRankedEnrichment31t_web as tools
import numpy as np
import cPickle as pickle
import creates_files_for_funcassociate as dicter
import graphics_gene_enrichment as gr

def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

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

#leave these parameters fixed
cutoff=20
top_cutoff=200
burnin=int(1*1e5)
steps=int(1*1e6)

#you can change these parameters
weight_kind=1 #1 is what we used
category_kind=0

#really change these parameters
which=1 #0: HDF, 1: 3-cell vs 2-cell - GO, 2: 3-cell vs 2-cell - MSigDB
s=1#0: for our method, 1: for GOing Bayesian
min_level=30 #10: 10% cutoff, 20: 20% cutoff, 30: 30% cutoff
belief=5 #use only 2 and 5 and 10
alpha_beta_top=0.3
alpha_beta_tops=[0.1,0.2,0.3,0.5]


gene_files=['bkz_gene_list.csv','3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12.txt','datasets/disease_glaucoma.txt','datasets/disease_obesity.txt','3-cell-hep-d12_versus_CS-d12_ID.txt','datasets/disease_adenocarcinoma_of_lung.txt','datasets/disease_ischemic_cardiomyopathy.txt']
out_strs=['HDF','LIVER','LIVER','DISEASE','DISEASE','LIVER','DISEASE','DISEASE']

#for the heat plot:
n_intervals=10 #number of different groups of active categories

#code, don't change
termss=[[[0 for iii in range(4)] for ii in range(len(alpha_beta_tops))] for i in range(8)]
meanss=[[[0 for iii in range(4)] for ii in range(len(alpha_beta_tops))] for i in range(8)]

ii=0
for which in [6]:         
    jj=0
    for alpha_beta_top in alpha_beta_tops:
        kk=0
        for belief in [100]:       
                                                            
            m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                            0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, alpha_beta_top, 0.2,20,20,gene_files[which],0,
                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                            'output1/', out_strs[which])
            if which==2:
                m.go_file='msigdb.txt'                   
                                
            (genes_list, level_list) = m.get_sets_and_genes()
            (T, unique_genes, large_goids, term_names) = m.load_all_new(genes_list)
            lug=len(unique_genes)
            m.min_level_for_activation=correct_cutoff(min_level/100.,unique_genes,level_list,genes_list,m.kind_of_list)
            Tset=set(range(len(T)))
            cc=[len(T[i]) for i in range(len(T))]
            (G,_, dict_G)=m.getG(unique_genes,genes_list,T,Tset,level_list)
            glist=[genes_list.index(unique_genes[i]) for i in range(len(unique_genes))]
            levels=[level_list[glist[i]] for i in range(len(unique_genes))]
            dummy=glist[:]
            dummy.sort()
            levels_ordinal = [dummy.index(g) for g in glist]
            
            if m.out_str=='DISEASE':
                m.out_str=m.gene_file.split('.')[0].split('datasets/disease_')[1]
            
            #load data
            if alpha_beta_top==0.5:
                addon=""
            else:
                addon='_alphabetatop'+str(int(round(m.alpha_beta_top*1000)))
            try:
                nsim=5
                if m.go_file=='msigdb.txt':
                    f=open('saved_data/multiple_runs_new_msigdb_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
                else:
                    f=open('saved_data/multiple_runs_new_go_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
            except IOError:
                try:
                    nsim=20
                    if m.go_file=='msigdb.txt':
                        f=open('saved_data/multiple_runs_new_msigdb_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
                    else:
                        f=open('saved_data/multiple_runs_new_go_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
                except IOError:
                    print "There exists no saved dataset for this choice of parameters"
            [C,MCMC_distr]=pickle.load(f)
            f.close()
            
            #data analysis
            terms=[]
            distr=[]
            for i in range(nsim):
                for j in range(len(C[i])):
                    try:
                        ind=terms.index(C[i][j])  
                    except ValueError:
                        ind=len(terms)
                        terms.append(C[i][j])
                        distr.append([0]*nsim)
                    distr[ind][i]=MCMC_distr[i][j]
                        
            means,stds=[],[]
            for i in range(len(distr)):
                means.append(np.mean(distr[i]))
                stds.append(np.std(distr[i]))
                
            ind=sorted(range(len(distr)), reverse=True, key=lambda k: means[k])
            termss[ii][jj][kk]=[terms[i] for i in ind]
            meanss[ii][jj][kk]=[means[i] for i in ind]
            stds=[stds[i] for i in ind]
            
            kk+=1
        jj+=1
    ii+=1
