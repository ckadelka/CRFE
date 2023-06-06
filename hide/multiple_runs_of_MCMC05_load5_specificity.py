import GeneSetRankedEnrichment32d_web as tools
import numpy as np
import cPickle as pickle
import graphics_gene_enrichment as gr
from matplotlib.pyplot import *
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

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
which=6 #0: HDF, 1: 3-cell vs 2-cell - GO, 2: 3-cell vs 2-cell - MSigDB
min_level=0.3 #10: 10% cutoff, 20: 20% cutoff, 30: 30% cutoff
alpha_beta_top=0.5

nsim=20

gene_files=['bkz_gene_list.csv','3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12.txt','datasets/disease_glaucoma.txt','datasets/disease_obesity.txt','3-cell-hep-d12_versus_CS-d12_ID.txt','datasets/disease_adenocarcinoma_of_lung.txt','datasets/disease_ischemic_cardiomyopathy.txt','3-cell-hep-d12_versus_2-cell-hep-d12_pos_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12_neg_ID.txt']
out_strs=['HDF','LIVER','LIVER','DISEASE','DISEASE','LIVER','DISEASE','DISEASE','Lpos','Lneg']

#for the category plot:
n_intervals=10 #number of different groups of active categories
threshold=0.1

termss,meanss,stdss,data,cats,legend,cats,avg_lens=[],[],[],[],[],[],[],[]

colorlist=[gr.rgb(212,175,255),gr.rgb(176,232,249),gr.rgb(255,224,191),gr.rgb(255,176,176),gr.rgb(169,235,196),gr.rgb(153,76,0)]

iiii=-1
for belief,s in [(10,0),(5,0),(2,0),(100,1)]:
    iiii+=1
    min_level=0.3
    #code, don't change
    m = tools.GeneSetRankedEnrichment(s, cutoff, top_cutoff, 0.1, 
                                    0.25, 1e-50, belief, weight_kind, category_kind, min_level, 'percentage',burnin, steps, alpha_beta_top,0.2,20,20,gene_files[which],0,
                                    'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                    'output1/', out_strs[which])
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
    
    if belief==9:
        break
    
    min_level=int(min_level*100)
    
    if m.out_str=='DISEASE':
        m.out_str=m.gene_file.split('.')[0].split('datasets/disease_')[1]
    
    #load data
    if alpha_beta_top==0.5:
        addon=""
    else:
        addon='_alphabetatop'+str(int(round(m.alpha_beta_top*1000)))
    try:
        if m.go_file=='msigdb.txt':
            f=open('saved_data/multiple_runs_new_msigdb_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
        else:
            f=open('saved_data/multiple_runs_new_go_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
    except IOError:
        try:
            if m.go_file=='msigdb.txt':
                f=open('saved_data/multiple_runs_new_msigdb_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
            else:
                f=open('saved_data/multiple_runs_new_go_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
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
    terms=[terms[i] for i in ind]
    means=[means[i] for i in ind]
    stds=[stds[i] for i in ind]
    
    termss.append(terms)
    meanss.append(means)
    stdss.append(stds)
    data.append([C,MCMC_distr])
    #cats.append(gr.categories(C,T,MCMC_distr,levels_ordinal,lug-len(G[-1]),lug,n_intervals,cutoff=threshold,kind=3))
    legend.append('CRFE (belief='+str(belief)+')' if s==0 else 'MGSA')
    cats.append(gr.categories(C,m.T,MCMC_distr,m.levels_ordinal,m.number_of_genes-len(m.G[-1]),m.number_of_genes,10,cutoff=threshold,kind=3))

    cc=[len(m.T[i]) for i in terms]
    avg_len=[]
    sum_len=0
    count=0
    for el in cc:
        count+=1
        sum_len+=el
        avg_len.append(sum_len*1./count)
    avg_lens.append(avg_len)

#gr.pert_vs_unpert2a_ratio(termss,m.G,m.T,legend,meanss,fontsize=14,colors=[gr.rgb(212,175,255),gr.rgb(176,232,249),gr.rgb(255,224,191),gr.rgb(255,176,176)],coarsity=1,show_annotations=0)
#gr.plotte_categories_explained2(cats,m,len(m.G[-1]),min_level,legend=legend,title="",identifier="",MCMC_threshold=threshold,colors=[gr.rgb(222,189,255),gr.rgb(191,235,250),gr.rgb(224,237,199),gr.rgb(255,191,191)])
#a=gr.plotte_categories_explained2(cats,m,len(m.G[-1]),min_level,legend=legend,fontsize=15,title="",identifier="",MCMC_threshold=threshold,colors=[gr.rgb(212,175,255),gr.rgb(176,232,249),gr.rgb(255,224,191),gr.rgb(255,176,176)])
#plt.text(1.8,-0.1,'Perturbed genes grouped by rank')
#plt.title('Adenocarcinoma of Lung')
#gr.jaccard_similarity_tau(termss,m.G,m.T,legend,meanss,max_nr_terms=100,colors=['k','b','r','g',[0,0.2,0]])
#gr.jaccard_similarity_tau(termss,m.G,m.T,legend,meanss=[],max_nr_terms=100,colors=['k','b','r','g',[0,0.2,0]])

colorlist=[gr.rgb(212,175,255),gr.rgb(176,232,249),gr.rgb(255,224,191),gr.rgb(255,176,176),gr.rgb(255,105,180),gr.rgb(153,76,0),gr.rgb(55,240,40),gr.rgb(100,100,255)]


for el in [C_FA,C_GOrilla,C_GSEA,C_REVIGO]:
    cc=[len(m.T[i]) for i in el]
    avg_len=[]
    sum_len=0
    count=0
    for el in cc:
        count+=1
        sum_len+=el
        avg_len.append(sum_len*1./count)
    avg_lens.append(avg_len)

N_max=40
for i in range(len(avg_lens)):
    plt.plot(range(1,min(N_max,len(avg_lens[i]))+1),avg_lens[i][:N_max],color=colorlist[i],lw=2)
exp_value=sum([len(t) for t in m.T])/len(m.T)
plt.plot(np.array([0,N_max]),exp_value*np.ones(2),color='k',ls='--',lw=2)
font = {'size'   : 15}
matplotlib.rc('font', **font) 
#plt.legend(['CRFE (belief=10)', 'CRFE (belief=5)', 'CRFE (belief=2)', 'MGSA','FuncAssociate2.0','GOrilla','GSEA','GSEA+REVIGO'])
plt.ylabel('Average Number of Annotations')
plt.xlabel('Number of Considered Processes')
#plt.title('Liver')
plt.title('Adenocarcinoma of Lung')
