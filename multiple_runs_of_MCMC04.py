import GeneSetRankedEnrichment31q_web as tools
import numpy as np
import cPickle as pickle

def jaccard_total_genes(C,T,MCMC_distr,cutoff=0):
    import numpy as np
    res=[[0 for i in range(len(C))] for j in range(len(C))]
    cutoffs=[sum([distr>=cutoff for distr in MCMC_distr[i]]) for i in range(len(MCMC_distr))]
    T_C=[set.union(*[set([])] + [T[c] for c in C[i][:cutoffs[i]]]) for i in range(len(C))] 
    for i in xrange(len(C)):
        if cutoffs[i]==0:
            continue
        for j in xrange(i+1,len(C)):
            if i==j:
                continue
            if cutoffs[j]==0:
                continue
            res[i][j]=len(T_C[i]&T_C[j])*1./len(T_C[i]|T_C[j])
            res[j][i]=res[i][j]
    #to remove the zeros on the diagonale
    for i in xrange(len(C)):
        res[i].pop(i)
    return (np.matrix(res),cutoffs)

def jaccard_term_by_term(C,T,MCMC_distr,cutoff=0):
    #useless, right now
    import numpy as np
    cutoffs=[sum([distr>=cutoff for distr in MCMC_distr[i]]) for i in range(len(MCMC_distr))]
    res=[]
    for i in range(len(C)):
        res.append(np.mean(m.jaccard(C[i][:cutoffs[i]],T,[],False))*2*cutoffs[i]/(cutoffs[i]-1))
    return (np.array(res),cutoffs)

def nr_explained(C,T,MCMC_distr,levels_ordinal,cutoff=0):
    import numpy as np
    res=[[0 for i in range(lug)] for i in range(len(C))]
    counter=[0]*len(C)
    for i in range(len(C)):
        for j in range(len(C[i])):
            if MCMC_distr[i][j]<cutoff:
                break
            for gene in T[C[i][j]]:
                res[i][levels_ordinal[gene]]+=1
            counter[i]+=1
    explained=[[0 for i in range(2)] for i in range(len(C))]
    for i in range(nsim):
        explained[i][0]=sum([el>0 for el in res[i][:lug-len(G[-1])]])
        explained[i][1]=sum([el>0 for el in res[i][lug-len(G[-1]):]])
    return (np.matrix(explained),counter)

def categories_explained(C,T,MCMC_distr,levels_ordinal,n_intervals=10,cutoff=0):
    import numpy as np
    res=[[0 for i in range(lug)] for i in range(len(C))]
    counter=[0]*len(C)
    for i in range(len(C)):
        for j in range(len(C[i])):
            if MCMC_distr[i][j]<cutoff:
                break
            for gene in T[C[i][j]]:
                res[i][levels_ordinal[gene]]+=1
            counter[i]+=1
        
    width=(lug-len(G[-1]))*1./n_intervals
    explained=[[0 for i in range(n_intervals)] for i in range(len(C))]
    for i in range(nsim):
        for j in range(n_intervals):
            explained[i][j]=sum([el>0 for el in res[i][int(width*j):int(width*(j+1))]])
    means=np.mean(np.matrix(explained),0)
    means=np.array([means[0,i] for i in range(n_intervals)])
    stds=np.std(np.matrix(explained),0)
    stds=np.array([stds[0,i] for i in range(n_intervals)])
    return (means,stds,width)

def plotte_categories_explained(data_vec,legend="",title=""):
    import numpy as np
    import matplotlib.pyplot as plt
    
    colors=['b','r']#,'g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    
    N = len(data_vec)
    n_intervals=len(data_vec[0][0])
    
    ind = np.arange(n_intervals)  # the x locations for the groups
    width=0.85/N       # the width of the bars
    
    rects=[]
    fig, ax = plt.subplots()
    for i in range(N):
        rects.append(ax.bar(ind+i*width, data_vec[i][0]/data_vec[i][2], width, yerr=data_vec[i][1]/data_vec[i][2], color=colors[i]))
    
    # add some
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,0,1))
    ax.set_ylabel('Proportion of perturbed genes that are explained')
    ax.set_xlabel('All perturbed genes grouped by rank')
    ax.set_title(title)
    ax.set_xticks(ind+0.425)
    ax.set_xticklabels( [str(i+1)+'/'+str(n_intervals) for i in range(n_intervals)] )
    #ax.get_children()[7+2*n_intervals].set_color('k')
    
    for i in range(N):
        ax.plot(np.arange(n_intervals+1),[np.mean(data_vec[i][0])/data_vec[i][2]]*(n_intervals+1), color=colors[i], linestyle='--',linewidth=2)
    
    if legend!="":
        ax.legend( rects, legend )
    filename='aaasd.eps'
    plt.savefig(filename, bbox_inches=0)
    #plt.show()
    return ax

def freq_explained(C,T,MCMC_distr,levels_ordinal,cutoff=0):
    import numpy as np
    res=[[0 for i in range(lug)] for i in range(len(C))]
    counter=[0]*len(C)
    for i in range(len(C)):
        for j in range(len(C[i])):
            if MCMC_distr[i][j]<cutoff:
                break
            for gene in T[C[i][j]]:
                res[i][levels_ordinal[gene]]+=1
            counter[i]+=1
    explained=[[[0 for ii in range(6)] for i in range(2)] for j in range(len(C))]
    for i in range(nsim):
        for j in range(5):
            explained[i][0][j]=sum([el==j for el in res[i][:lug-len(G[-1])]])
            explained[i][1][j]=sum([el==j for el in res[i][lug-len(G[-1]):]])
        explained[i][0][-1]=sum([el>4 for el in res[i][:lug-len(G[-1])]])
        explained[i][1][-1]=sum([el>4 for el in res[i][lug-len(G[-1]):]])
    return (np.array(explained),counter)

def different_cutoffs(C,T,MCMC_distr,levels_ordinal,N=20):
    cutoffs=np.linspace(0,1,N)
    mean_JI=np.zeros(N)
    std_JI=np.zeros(N)
    mean_pert=np.zeros(N)
    std_pert=np.zeros(N)
    mean_unpert=np.zeros(N)
    std_unpert=np.zeros(N)
    mean_nr_terms=np.zeros(N)
    std_nr_terms=np.zeros(N)
    mean_freq=[np.zeros((2,6)) for i in range(N)]
    std_freq=[np.zeros((2,6)) for i in range(N)]
    i=0
    for cutoff in cutoffs:
        ((mean_JI[i],std_JI[i]),(mean_pert[i],std_pert[i]),(mean_unpert[i],std_unpert[i]),(mean_nr_terms[i],std_nr_terms[i]),(mean_freq[i],std_freq[i]))=stats(C,T,MCMC_distr,levels_ordinal,cutoff)
        i+=1
    return ((mean_JI,std_JI),(mean_pert,std_pert),(mean_unpert,std_unpert),(mean_nr_terms,std_nr_terms),(mean_freq,std_freq))
    
def stats(C,T,MCMC_distr,levels_ordinal,cutoff=0):
    a=jaccard_total_genes(C,T,MCMC_distr,cutoff)
    b=nr_explained(C,T,MCMC_distr,levels_ordinal,cutoff)
    c=freq_explained(C,T,MCMC_distr,levels_ordinal,cutoff)
    return ((np.mean(a[0]),np.std(a[0])),(np.mean(b[0],0)[0,0],np.std(b[0],0)[0,0]),(np.mean(b[0],0)[0,1],np.std(b[0],0)[0,1]),(np.mean(a[1]),np.std(a[1])),(np.mean(c[0],0),np.std(c[0],0)))

def plotte(data,N,what=[0],data1=[]):
    import matplotlib.pyplot as plt
    x=np.linspace(0,1,N)
    if data1==[]:
        colors=['b','r','g','m']
    else:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    for ind in what:
        plt.plot(x,data[ind][0],linestyle='-',color=colors[ind])
        plt.plot(x,data[ind][0]+data[ind][1],linestyle='--',color=colors[ind])
        plt.plot(x,data[ind][0]-data[ind][1],linestyle='--',color=colors[ind])
        if data1!=[]:
            plt.plot(x,data1[ind][0],linestyle='-',color=colors[ind+4])
            plt.plot(x,data1[ind][0]+data1[ind][1],linestyle='--',color=colors[ind+4])
            plt.plot(x,data1[ind][0]-data1[ind][1],linestyle='--',color=colors[ind+4])

def plotte_freq(data,N,what=[0]):
    import matplotlib.pyplot as plt
    x=np.linspace(0,1,N)
    colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    perturbed=1
    means=[]
    stds=[]
    for ind in range(len(data[4][0][0][0])):
        means.append(np.array([data[4][0][i][1-perturbed][ind] for i in range(N)]))
        stds.append(np.array([data[4][1][i][1-perturbed][ind] for i in range(N)]))
    for ind in what:
        plt.plot(x,means[ind],linestyle='-',color=colors[ind])
        plt.plot(x,means[ind]+stds[ind],linestyle='--',color=colors[ind])
        plt.plot(x,means[ind]-stds[ind],linestyle='--',color=colors[ind])

nsim=20

cutoff=5
top_cutoff=200
param_learning=1
MCMC=1
belief=5
burnin=int(1*1e5)
steps=int(1*1e6)
data2=[[[0 for iii in range(2)] for ii in range(2)] for i in range(3)]
for hdhsfs in [0,1,2]:
    for s in [0,1]:
        for weight_kind in [0,1]:
            #if s==0 and weight_kind==0 or weight_kind==0 and s==1 and hdhsfs==0:
            #    continue 
            if hdhsfs==0:
                m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                            0.25, 1e-50, belief, weight_kind, 0.559, 1000, 1000, burnin, steps, 0.5,0.2,20,20,'bkz_gene_list.csv',0,
                                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                            'output1/', 'HDF')
            elif hdhsfs==1:
                m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                        0.25, 1e-50, belief, weight_kind, 0.08016865, 1000, 1000, burnin, steps, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt',0,
                                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                            'output/', 'LIVER')
            elif hdhsfs==2:
                m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                        0.25, 1e-50, belief, weight_kind, 0.08016865, 1000, 1000, burnin, steps, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12.txt',0,
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
            
            #load data
            try:
                nsim=20
                if m.go_file=='msigdb.txt':
                    f=open('saved_data/multiple_runs_msigdb_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_nsim'+str(nsim)+'.txt','rb')
                else:
                    f=open('saved_data/multiple_runs_go_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_nsim'+str(nsim)+'.txt','rb')
            except IOError:
                nsim=50
                if m.go_file=='msigdb.txt':
                    f=open('saved_data/multiple_runs_msigdb_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_nsim'+str(nsim)+'.txt','rb')
                else:
                    f=open('saved_data/multiple_runs_go_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_nsim'+str(nsim)+'.txt','rb')
            [C,MCMC_distr]=pickle.load(f)
            f.close()
            
            print hdhsfs,s,weight_kind
            data2[hdhsfs][s][weight_kind]=categories_explained(C,T,MCMC_distr,levels_ordinal,10,0.1)