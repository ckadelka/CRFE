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
burnin=int(1*1e2)
steps=int(1*1e3)

#you can change these parameters
weight_kind=1 #1 is what we used
category_kind=0

#really change these parameters
which=1 #0: HDF, 1: 3-cell vs 2-cell - GO, 2: 3-cell vs 2-cell - MSigDB
s=0#0: for our method, 1: for GOing Bayesian
min_level=30 #10: 10% cutoff, 20: 20% cutoff, 30: 30% cutoff
belief=2 #use only 2 and 5 and 10
alpha_beta_top=0.3

gene_files=['bkz_gene_list.csv','3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12.txt','datasets/disease_glaucoma.txt','datasets/disease_obesity.txt','3-cell-hep-d12_versus_CS-d12_ID.txt','datasets/disease_adenocarcinoma_of_lung.txt','datasets/disease_ischemic_cardiomyopathy.txt']
out_strs=['HDF','LIVER','LIVER','DISEASE','DISEASE','LIVER','DISEASE','DISEASE']

#for the heat plot:
n_intervals=10 #number of different groups of active categories

#code, don't change
                
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
    nsim=20
    if m.go_file=='msigdb.txt':
        f=open('saved_data/multiple_runs_new_msigdb_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
    else:
        f=open('saved_data/multiple_runs_new_go_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','rb')
except IOError:
    try:
        nsim=2
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
terms=[terms[i] for i in ind]
means=[means[i] for i in ind]
stds=[stds[i] for i in ind]

##Plots
def pert_vs_unpert(terms,G,T,posteriors=[]):
    import matplotlib.pyplot as plt
    import numpy as np
    pep=[] #proportion explained perturbed genes
    peu=[] #proportion explained unperturbed genes
    MCMC_threshold=[] #at a level of this threshold
    P=set.union(*[set([])] + [G[i] for i in range(len(G)-1)])
    U=set(G[-1])
    sizeP=len(P)
    sizeU=len(U)
    E=set([])
    for i in range(len(terms)):
        E=set.union(E,set(T[terms[i]]))
        pep.append(len(E&P))
        peu.append(len(E)-pep[-1])
        if posteriors!=[]:
            MCMC_threshold.append(posteriors[i])
    plt.plot(np.array(peu)*1./sizeU,np.array(pep)*1./sizeP)
    #plt.plot(MCMC_threshold,np.array(pep)*1./sizeP,MCMC_threshold,np.array(peu)*1./sizeU)



#Generate heatplot
#width=(lug-len(G[-1]))*1./n_intervals
#G=[[] for i in range(n_intervals+1)]
#for gene in range(lug):
#    if levels[gene]<m.min_level_for_activation:
#        G[-1].append(gene)
#    else:
#        G[int(levels_ordinal[gene]/width)].append(gene)
#
#terms=terms[:20]        
#                        
#ccG=[len(g) for g in G]
#ccC=[cc[c] for c in terms]
#
#mat=[]
#ps=[]
#for i in range(len(terms)):
#    mat.append([])
#    ps.append([])
#    for j in range(len(G)):
#        mat[-1].append(len(set(G[j])&set(T[terms[i]])))
#        ps[-1].append(1-sum([choose(ccC[i],v)*(ccG[j]*1./lug)**v*(1-ccG[j]*1./lug)**(ccC[i]-v) for v in xrange(mat[i][j])]))
#
#mat_rat=[[max(mat[i][j]*1./ccC[i],0.000001) for j in xrange(len(mat[i]))] for i in xrange(len(mat))]
#expected=[val*1./lug for val in ccG]
#mean_levels=m.mean_gene_level(terms,[],T,levels)
#p_values=m.p_values_of_list(terms, T, G)
#
#fold=[[round(1./expected[j]*mat_rat[i][j],2) for j in xrange(len(mat[i]))] for i in xrange(len(mat))]
#
#filename='Heatplots/heatplot_test_foldchange_'+m.out_str+'_s'+str(m.s)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)
#
#f_out = open(filename+'.txt', 'w')
#string='ID in code\t GO term \t MCMC density (mean)\t MCMC density (std)\t Mean Gene Level\t Hypergeometric p-value of term\t Explained genes\t Total annotations per term & Total genes per category\t'
#for val in ccG:
#    string+=str(val)+'\t'
#f_out.write(string[:-1]+'\n')
#for i in range(len(fold)):
#    string=str(terms[i])+'\t'+term_names[terms[i]]+'\t'+str(round(means[i],4)) +'\t'+str(round(stds[i],4)) +'\t'+str(round(mean_levels[i],4))+'\t'+str(p_values[i])+'\t'+str(cc[terms[i]]-len(set(G[-1])&set(T[terms[i]])))+'\t'+str(cc[terms[i]])+'\t'
#    for val in fold[i]:
#        string+=str(val)+'\t'
#    f_out.write(string[:-1]+'\n')
#f_out.flush()
#f_out.close()
#
#filename='Heatplots/heatplot_foldchange_pvalue_'+m.out_str+'_s'+str(m.s)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)
#
#f_out = open(filename+'.txt', 'w')
#string='ID in code\t GO term \t MCMC density (mean)\t MCMC density (std)\t Mean Gene Level\t Hypergeometric p-value of term\t Explained genes\t Total annotations per term & Total genes per category\t'
#for val in ccG:
#    string+=str(val)+'\t'
#f_out.write(string[:-1]+'\n')
#for i in range(len(fold)):
#    string=str(terms[i])+'\t'+term_names[terms[i]]+'\t'+str(round(means[i],4)) +'\t'+str(round(stds[i],4)) +'\t'+str(round(mean_levels[i],4))+'\t'+str(p_values[i])+'\t'+str(cc[terms[i]]-len(set(G[-1])&set(T[terms[i]])))+'\t'+str(cc[terms[i]])+'\t'
#    for val in ps[i]:
#        string+=str(val)+'\t'
#    f_out.write(string[:-1]+'\n')
#f_out.flush()
#f_out.close()

##Create gene list with information about which gene is explained by what process
#upto=sum([mean>0.5 for mean in means])
#explained_by=[[] for i in range(lug)]
#for count,t in enumerate(terms[:upto]):
#    for gene in T[t]:
#        explained_by[gene].append((t,count))
#
#(symbol2num,num2symbol) = dicter.dict_transform(m.gene_file)
#output=["" for i in range(lug)]
#ranks=[levels_ordinal[i] for i in range(lug)]
#sorted_ranks=sorted(range(lug), reverse=False, key=lambda k: ranks[k])
#for count,i in enumerate(sorted_ranks):
#    output[count]=unique_genes[i]+'\t'+num2symbol[unique_genes[i]]
#    if explained_by[i]!=[]:
#        output[count]+="\t"+str(round(means[explained_by[i][0][1]],3))
#        for j in range(len(explained_by[i])):
#            output[count]+="\t"+term_names[explained_by[i][0][0]]+" ("+str(round(means[explained_by[i][j][1]],3))+")"
#            
#g=file('which_genes_explained.txt','w')
#for out in output:
#    g.write(out+'\n')
#g.close()