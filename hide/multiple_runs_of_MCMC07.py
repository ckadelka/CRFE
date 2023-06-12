import GeneSetRankedEnrichment31q_web as tools
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

def jaccard_term_by_term(C,T,MCMC_distr,cutoff=0,n_intervals=10,kind=2):
    #category: 0 = perturbed genes, 1 = unperturbed genes
    import numpy as np
    if kind==1:
        cutoffs=[cutoff for i in range(len(MCMC_distr))]
    else:
        cutoffs=[sum([distr>=cutoff for distr in MCMC_distr[i]]) for i in range(len(MCMC_distr))]
    res=[]
    for i in range(len(C)):
        mat=m.jaccard(C[i][:cutoffs[i]],T,[],False)
        maxmat=[max(row) for row in mat]
        res.append(tools.hist(maxmat,n_intervals,0,1)[0])
    res_mat=np.matrix(res)
    res_mat=res_mat*1./np.sum(res_mat,1)
    return (res_mat,cutoffs)

def jaccard_term_by_term_cat(C,T,MCMC_distr,G,category=[0],cutoff=0,n_intervals=10,kind=2):
    #category: 0 = perturbed genes, 1 = unperturbed genes
    import numpy as np
    if kind==1:
        cutoffs=[cutoff for i in range(len(C))]
    else:
        cutoffs=[sum([distr>=cutoff for distr in MCMC_distr[i]]) for i in range(len(MCMC_distr))]
    cat=[]
    if 0 in category:
        cat.extend(range(len(G)-1))
    if 1 in category:
        cat.append(len(G)-1)
    res=[]
    for i in range(len(C)):
        mat=m.jaccard_comp(C[i][:cutoffs[i]],cat,G,T,3)
        if np.size(mat)>1:
            maxmat=[max(row) for row in mat]
        else:
            maxmat=[mat]
        res.append(tools.hist(maxmat,n_intervals,0,1)[0])
    res_mat=np.matrix(res)
    res_mat=res_mat*1./np.sum(res_mat,1)
    return (res_mat,cutoffs)
    
def plot_jaccard_term_by_term(C,T,MCMC_distr,G,category=[0,1],cutoff=0,n_intervals=10,kind=2,color='r',title=""):
    import matplotlib.pyplot as plt
    #kind: 1 = always use top nr_terms terms in each run, 2 = use all terms above nr_terms cutoff
    if category==[0,1]:
        res=jaccard_term_by_term(C,T,MCMC_distr,cutoff=cutoff,n_intervals=n_intervals,kind=kind)[0]
    else:
        res=jaccard_term_by_term_cat(C,T,MCMC_distr,G,category=category,cutoff=cutoff,n_intervals=n_intervals,kind=kind)[0]
    means=np.mean(res,0)
    stds=np.std(res,0)
    means=[means[0,i] for i in range(np.size(means,1))]
    stds=[stds[0,i] for i in range(np.size(stds,1))]
    stepwidth=1./n_intervals
    x=np.array([i*stepwidth for i in range(n_intervals)])
    fig, ax = plt.subplots()
    ax.bar(x, means, width=stepwidth, yerr=stds, error_kw={'ecolor':color, 'capsize':4}, color=color)
    plt.axis((0,1,0,1))
    ax.set_ylabel('Percent of Terms')
    ax.set_xlabel('Jaccard Index')
    if title=="":
        if kind==2:
            title=r'Jaccard Similarity of all terms with $\tau \geq$'+str(cutoff)
        else:
            title=r'Jaccard Similarity of the top '+str(cutoff)+' terms'
        if category==[0,1]:
            title+=' (all genes)'
        elif category==[0]:
            title+=' (perturbed genes)'
        elif category==[1]:
            title+=' (unperturbed genes)'    
    ax.set_title(title)
    ax.set_yticklabels([str(int(el*100))+'%' for el in ax.get_yticks()]) 

def plot_jaccard_term_by_term2(C,T,MCMC_distr,G,cutoff=0,n_intervals=10,kind=2,colors=[[1,0.5,0.5],[0.5,1,0.5]],title="",identifier="",SAVE=False):
    import matplotlib.pyplot as plt
    #kind: 1 = always use top nr_terms terms in each run, 2 = use all terms above nr_terms cutoff
    res,means,stds=[],[],[]
    res.append(jaccard_term_by_term_cat(C,T,MCMC_distr,G,category=[0],cutoff=cutoff,n_intervals=n_intervals,kind=kind))
    res.append(jaccard_term_by_term_cat(C,T,MCMC_distr,G,category=[1],cutoff=cutoff,n_intervals=n_intervals,kind=kind))
    mean_terms=np.mean(res[0][1])
    for j in [0,1]:
        means.append(np.mean(res[j][0],0))
        stds.append(np.std(res[j][0],0))
        means[j]=[means[j][0,i] for i in range(np.size(means[j],1))]
        stds[j]=[stds[j][0,i] for i in range(np.size(stds[j],1))]
    stepwidth=1./n_intervals
    x=np.array([i*stepwidth for i in range(n_intervals)])
    fig, ax = plt.subplots(2, sharex=True)
    for j in [0,1]:
        ax[j].bar(x, means[j], width=stepwidth, yerr=stds[j], error_kw={'ecolor':colors[j], 'capsize':4}, color=colors[j])
        ax[j].set_ylabel('Perturbed genes' if j==0 else 'Unperturbed genes')
        ax[j].set_ylim((0,1))   
        ax[j].set_yticklabels([str(int(el*100))+'%' for el in ax[j].get_yticks()])
        ax[j].set_xlim((0,1)) 
    ax[1].set_xlabel('Jaccard Index')
    if kind==2:
        title=r'Jaccard similarity of all terms with $\tau \geq$'+str(cutoff)+ ' (~'+str(mean_terms)+' terms)'+(' for\n'+title if title!="" else "")
    else:
        title=r'Jaccard similarity of the top '+str(cutoff)+' terms' + (' for\n'+title if title!="" else "")
    ax[0].set_title(title)
    if SAVE:
        filename="jaccard_sim_"+identifier
        filename+="_cut"+str(int(cutoff*100) if kind==2 else cutoff)+"_k"+str(kind)+"_N"+str(n_intervals)
        filename+='.eps'
        plt.savefig(filename, bbox_inches=0)
    return (means,stds)
    
    
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

def nr_annotations(C,T,MCMC_distr,levels_ordinal,n_intervals=10,cutoff=0):
    import numpy as np
    res=[0 for i in range(lug)]
    for t in T:
        for gene in t:
            res[levels_ordinal[gene]]+=1
        
    width=(lug-len(G[-1]))*1./n_intervals
    explained=[0 for i in range(n_intervals)]
    for j in range(n_intervals):
        explained[j]=sum([el for el in res[int(width*j):int(width*(j+1))]])
    return (explained,width)

def plotte_categories_explained(data_vec,legend="",title=""):
    import numpy as np
    import matplotlib.pyplot as plt
    
    colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    
    N = len(data_vec)
    n_intervals=len(data_vec[0][0])
    
    ind = np.arange(n_intervals)  # the x locations for the groups
    width=0.85/N       # the width of the bars
    
    rects=[]
    fig, ax = plt.subplots()
    for i in range(N):
        rects.append(ax.bar(ind+i*width, data_vec[i][0]/data_vec[i][2], width, yerr=data_vec[i][1]/data_vec[i][2], error_kw={'ecolor':str(colors[i]), 'capsize':4}, color=colors[i]))
    
    # add some
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,0,1))
    ax.set_ylabel('Proportion of perturbed genes that are explained')
    ax.set_xlabel('All perturbed genes grouped by rank')
    if title!="":
        ax.set_title(title)
    ax.set_xticks(ind+0.425)
    ax.set_xticklabels( [str(i+1) for i in range(n_intervals)] )
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
weight_kind=0
data1=[[[[0 for iiii in range(3)] for iii in range(3)] for ii in range(2)] for i in range(4)]
Cs=[[[[0 for iiii in range(3)] for iii in range(3)] for ii in range(2)] for i in range(4)]
MCMC_distrs=[[[[0 for iiii in range(3)] for iii in range(3)] for ii in range(2)] for i in range(4)]
beliefs=[2,5,10]
min_levels=[10,20,30]
for hdhsfs in [0,1,2]:
    for s in [0]:
        count=0
        for ind in [1]:
            for j in [0,1,2]:
                #if s==0 and weight_kind==0 or weight_kind==0 and s==1 and hdhsfs==0:
                #    continue 
                belief=beliefs[j]
                min_level=min_levels[ind]
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
                elif hdhsfs==3:
                    m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                0.25, 1e-50, belief, weight_kind, 0.559, 1000, 1000, burnin, steps, 0.5,0.2,20,20,'datasets/disease_alzheimers.txt',0,
                                'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                'output1/', 'DISEASE')                                           
                                                                                                                                        
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
                    m.out_str=m.gene_file.split('_')[1][:-4]             
                                        
                #load data
                try:
                    nsim=20
                    if m.go_file=='msigdb.txt':
                        f=open('saved_data/multiple_runs_msigdb_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
                    else:
                        f=open('saved_data/multiple_runs_go_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
                except IOError:
                    nsim=50
                    if m.go_file=='msigdb.txt':
                        f=open('saved_data/multiple_runs_msigdb_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
                    else:
                        f=open('saved_data/multiple_runs_go_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
                [Cs[hdhsfs][s][ind][j],MCMC_distrs[hdhsfs][s][ind][j]]=pickle.load(f)
                f.close()
                
                if hdhsfs==0:
                    title='HDF data'
                    identifier='hdf'
                elif hdhsfs==1:
                    title='3-cell vs. 2-cell - GO'
                    identifier='liver_go'
                elif hdhsfs==2:
                    title='3-cell vs. 2-cell - MSigDB'
                    identifier='liver_msigdb'
                elif hdhsfs==3:
                    title='Alzheimers - GO'
                    identifier='alzheimers_go'
                if s==0:
                    title+=", MCMC"
                elif s==1:
                    title+=", GOing Bayesian"
                title+=" ("+str(min_level)+"% perturbed genes)"
                identifier+='_s'+str(s)+'_pert'+str(min_level)
                a=plot_jaccard_term_by_term2(Cs[hdhsfs][s][ind][j],T,MCMC_distrs[hdhsfs][s][ind][j],G,cutoff=20,n_intervals=10,kind=1,colors=[[1,0.7,0.7],[0.7,1,0.7]],title=title,identifier=identifier,SAVE=1)
                
                print hdhsfs,s,min_level,belief
                data1[hdhsfs][s][ind][j]=categories_explained(Cs[hdhsfs][s][ind][j],T,MCMC_distrs[hdhsfs][s][ind][j],levels_ordinal,10,0.01)
                count+=1