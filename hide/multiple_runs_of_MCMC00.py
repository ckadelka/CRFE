import GeneSetRankedEnrichment31q_web as tools
import numpy as np
import cPickle as pickle

def jaccard(liste):
    res=[[0 for i in xrange(len(liste))] for j in xrange(len(liste))]
    for i in xrange(len(liste)):
        for j in xrange(len(liste)):#i+1,len(liste)):
            if i==j:
                continue
            if T[liste[i]]!=[] or T[liste[j]]!=[]:
                res[i][j]=len(set(T[liste[i]])&set(T[liste[j]]))*1./len(set(T[liste[i]])|set(T[liste[j]]))
            else:
                res[i][j]=0
    return res

nsim=5

cutoff=5
top_cutoff=200
param_learning=1
MCMC=1

m = tools.GeneSetRankedEnrichment(2, 0, cutoff, top_cutoff, 0.1, 
                            0.25, 1e-50, 5, 0, 0.08016865, 1000, 1000, 100000, int(2*1e5), 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12.txt',0,
                                   'msigdb.txt',
                                   'output/', 'LIVER')
(genes_list, level_list) = m.get_sets_and_genes()
(T, unique_genes, large_goids, term_names) = m.load_all_new(genes_list)
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
    (C[i],MCMC_distr[i])=m.runMe(0,0,param_learning,0,MCMC,0,0)

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

res=[[0 for i in range(lug)] for i in range(nsim)]
counter=[0]*nsim
for i in range(nsim):
    for j in range(len(C[i])):
        if MCMC_distr[i][j]<0.9:
            break
        for gene in T[C[i][j]]:
            res[i][levels_ordinal[gene]]+=1
        counter[i]+=1
        
n_intervals=10
width=lug*1./n_intervals
explained=[[0 for i in range(n_intervals)] for i in range(nsim)]
for i in range(nsim):
    for j in range(n_intervals):
        explained[i][j]=sum([el>0 for el in res[i][int(width*j):int(width*(j+1))]])

mat_expl=np.matrix(explained)
mean_expl=np.mean(mat_expl,0)
std_expl=np.std(mat_expl,0)

#save data
#f=open('saved_data/multiple_runs_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_nsim'+str(nsim)+'.txt','w+')
#pickle.dump([C,MCMC_distr], f)
#f.close()

##load data
#f=open('saved_data/multiple_runs_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_nsim'+str(nsim)+'.txt','rb')
#[C,MCMC_distr]=pickle.load(f)
#f.close()