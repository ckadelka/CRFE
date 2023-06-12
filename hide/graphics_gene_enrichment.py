# -*- coding: utf-8 -*-
#Collection of tools to create plots for the Gene Enrichment paper
import numpy as np
import GeneSetRankedEnrichment32b_web as tools
import math
from matplotlib.pyplot import *
matplotlib.rcParams['pdf.fonttype'] = 3
matplotlib.rcParams['ps.fonttype'] = 3


def rgb(vector,g=-1,b=-1):
    #Color for perturbed genes, rgb(39,93,144), color for explained genes: rgb(226,168,165)

    if type(vector==str) and vector in ['E','e','expl','explained']:
        vector=[226,168,165]
    elif type(vector==str) and vector in ['P','p','pert','perturbed']:
        vector=[39,93,144]      
    
    if type(vector) in [list,np.ndarray]:
        return list(np.array(vector)*1./255)
    elif type(vector)==int and 0<=g<=255 and 0<=b<=255:
        return list(np.array([vector,g,b])*1./255)


def rg(x):
    if x>=100:
        return int(round(x))
    elif x>=10:
        return round(x,1)
    elif x>=0.005:
        return round(x,2)
    else:
        return x

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

def jaccard_total_genes(C,T,MCMC_distr,cutoff=0,kind=2):
    res=[[0 for i in range(len(C))] for j in range(len(C))]
    if kind==1:
        cutoffs=[cutoff for i in range(len(C))]
    else:
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

def jaccard_term_by_term(C,T,MCMC_distr,m,cutoff=0,n_intervals=10,kind=2):
    #category: 0 = perturbed genes, 1 = unperturbed genes
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

def jaccard_term_by_term_cat(C,T,MCMC_distr,G,m,category=[0],cutoff=0,n_intervals=10,kind=2):
    #category: 0 = perturbed genes, 1 = unperturbed genes
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
        mat=m.jaccard_comp(C[i][:cutoffs[i]],cat,G,T,3)
        if np.size(mat)>1:
            maxmat=[max(row) for row in mat]
        else:
            maxmat=[mat]
        res.append(tools.hist(maxmat,n_intervals,0,1)[0])
    res_mat=np.matrix(res)
    res_mat=res_mat*1./np.sum(res_mat,1)
    return (res_mat,cutoffs)
    
def plot_jaccard_term_by_term(C,T,MCMC_distr,G,m,category=[0,1],cutoff=0,n_intervals=10,kind=2,color='r',title=""):
    import matplotlib.pyplot as plt
    #kind: 1 = always use top nr_terms terms in each run, 2 = use all terms above nr_terms cutoff
    if category==[0,1]:
        res=jaccard_term_by_term(C,T,MCMC_distr,m,cutoff=cutoff,n_intervals=n_intervals,kind=kind)[0]
    else:
        res=jaccard_term_by_term_cat(C,T,MCMC_distr,G,m,category=category,cutoff=cutoff,n_intervals=n_intervals,kind=kind)[0]
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

def plot_jaccard_term_by_term2(C,T,MCMC_distr,G,m,min_level=20,cutoff=0,n_intervals=10,kind=2,colors=[[1,0.7,0.7],[0.7,1,0.7]],title="",identifier="",SAVE=False,fontsize=12,STD=True):
    import matplotlib.pyplot as plt
    import matplotlib
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    #kind: 1 = always use top nr_terms terms in each run, 2 = use all terms above nr_terms cutoff
    res,means,stds=[],[],[]
    res.append(jaccard_term_by_term_cat(C,T,MCMC_distr,G,m,category=[0],cutoff=cutoff,n_intervals=n_intervals,kind=kind))
    res.append(jaccard_term_by_term_cat(C,T,MCMC_distr,G,m,category=[1],cutoff=cutoff,n_intervals=n_intervals,kind=kind))
    mean_terms=np.mean(res[0][1])
    print res
    for j in [0,1]:
        means.append(np.mean(res[j][0],0))
        stds.append(np.std(res[j][0],0))
        means[j]=[means[j][0,i] for i in range(np.size(means[j],1))]
        stds[j]=[stds[j][0,i] for i in range(np.size(stds[j],1))]
    stepwidth=1./n_intervals
    x=np.array([i*stepwidth for i in range(n_intervals)])
    fig, ax = plt.subplots(2, sharex=True)
    for j in [0,1]:
        if STD:
            ax[j].bar(x, means[j], width=stepwidth, yerr=stds[j], error_kw={'lw':2, 'ecolor':colors[j], 'capthick':2, 'capsize':8}, color=colors[j])
        else:
            ax[j].bar(x, means[j], width=stepwidth, color=colors[j])          
        ax[j].set_ylabel('Perturbed genes' if j==0 else 'Unperturbed genes')
        ax[j].set_ylim((0,1))   
        ax[j].set_yticklabels([str(int(el*100))+'%' for el in ax[j].get_yticks()])
        ax[j].set_xlim((0,1)) 
    ax[1].set_xlabel('Jaccard Index')
    #if kind==2:
    #    title=r'Jaccard similarity of all terms with $\tau \geq$'+str(cutoff)+ ' (~'+str(mean_terms)+' terms)'+(' for\n'+title if title!="" else "")
    #else:
    #    title=r'Jaccard similarity of top '+str(cutoff)+' terms' + (' for '+title if title!="" else "")
    ax[0].set_title(title)
    if SAVE:
        filename="jaccard_sim_"+identifier
        filename+="_cutoff"+str(m.cutoff)+"_topcutoff"+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+"_bel"+str(m.belief)+"_pert"+str(min_level)+"_cut"+str(int(cutoff*100) if kind==2 else cutoff)+"_k"+str(kind)+"_N"+str(n_intervals)
        filename+='.eps'
        plt.savefig(filename, bbox_inches=0)
    return (means,stds)

def plot_jaccard_term_by_term2a(C,T,MCMC_distr,G,m,min_level=30,cutoff=0,n_intervals=10,kind=2,colors=[],title="",identifier="",SAVE=False,fontsize=12,STD=True,bbox="",pad=0):
    import matplotlib.pyplot as plt
    
    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : fontsize}

    matplotlib.rc('font', **font)  
    
    #kind: 1 = always use top nr_terms terms in each run, 2 = use all terms above nr_terms cutoff
    res,means,stds=[],[],[]
    res.append(jaccard_term_by_term_cat(C,T,MCMC_distr,G,m,category=[0,1],cutoff=cutoff,n_intervals=n_intervals,kind=kind))
    for j in [0]:
        means.append(np.mean(res[j][0],0))
        stds.append(np.std(res[j][0],0))
        means[j]=[means[j][0,i] for i in range(np.size(means[j],1))]
        stds[j]=[stds[j][0,i] for i in range(np.size(stds[j],1))]
    stepwidth=1./n_intervals
    x=np.array([i*stepwidth for i in range(n_intervals)])
    fig, ax = plt.subplots()
    #fig.set_size_inches(15,8)
    for j in [0]:
        if STD:
            ax.bar(x, means[j], width=stepwidth, yerr=stds[j], error_kw={'lw':2, 'ecolor':colors[j], 'capthick':2, 'capsize':8}, color=colors[j])
        else:
            ax.bar(x, means[j], width=stepwidth, color=colors[j])            
        ax.set_ylim((0,1)) 
        ax.set_yticks([0.25,0.5,0.75,1])     
        ax.set_yticklabels([str(int(el*100))+'%' for el in ax.get_yticks()])
        ax.set_xlim((0,1)) 
    ax.set_xlabel('Redundancy')
    if fontsize>18:
        plt.subplots_adjust(left=0.15,top=0.85,bottom=0.2)
    if title!="":
        ax.set_title(title,y=1.03)
    if SAVE:
        filename="jaccard_sim_"+identifier
        filename+="_cutoff"+str(m.cutoff)+"_topcutoff"+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+"_bel"+str(m.belief)+"_pert"+str(min_level)+"_cut"+str(int(cutoff*100) if kind==2 else cutoff)+"_k"+str(kind)+"_N"+str(n_intervals)
        filename+='.eps'
        if bbox!="":
            plt.savefig(filename, pad_inches=pad, bbox_inches=bbox)
        else:
            plt.savefig(filename, pad_inches=pad)            
    return (means,stds)
    
def nr_explained(C,T,MCMC_distr,levels_ordinal,nr_perturbed,nr_genes,cutoff=0,kind=2):
    if kind==1:
        C=C[:cutoff]
    res=[[0 for i in range(nr_genes)] for i in range(len(C))]
    counter=[0]*len(C)
    for i in range(len(C)):
        for j in range(len(C[i])):
            if kind==2 and MCMC_distr[i][j]<cutoff:
                break
            for gene in T[C[i][j]]:
                res[i][levels_ordinal[gene]]+=1
            counter[i]+=1
    explained=[[0 for i in range(2)] for i in range(len(C))]
    for i in range(len(C)):
        explained[i][0]=sum([el>0 for el in res[i][:nr_perturbed]])
        explained[i][1]=sum([el>0 for el in res[i][nr_perturbed:]])
    return (np.matrix(explained),counter)

def categories_explained(C,T,MCMC_distr,levels_ordinal,nr_perturbed,nr_genes,n_intervals=10,cutoff=0):
    res=[[0 for i in range(nr_genes)] for i in range(len(C))]
    counter=[0]*len(C)
    for i in range(len(C)):
        for j in range(len(C[i])):
            if MCMC_distr[i][j]<cutoff:
                break
            for gene in T[C[i][j]]:
                res[i][levels_ordinal[gene]]+=1
            counter[i]+=1
        
    width=nr_perturbed*1./n_intervals
    explained=[[0 for i in range(n_intervals)] for i in range(len(C))]
    for i in range(len(C)):
        for j in range(n_intervals):
            explained[i][j]=sum([el>0 for el in res[i][int(width*j):int(width*(j+1))]])
    means=np.mean(np.matrix(explained),0)
    means=np.array([means[0,i] for i in range(n_intervals)])
    stds=np.std(np.matrix(explained),0)
    stds=np.array([stds[0,i] for i in range(n_intervals)])
    return (means,stds,width)

def categories(C,T,MCMC_distr,levels_ordinal,nr_perturbed,nr_genes,n_intervals=10,cutoff=0,kind=2):
    if kind==1:
        print "a",cutoff,len(C)
    res=[[0 for i in range(nr_genes)] for i in range(len(C))]
    counter=[0]*len(C)
    nr_unperturbed_explained=[0]*len(C)
    nr_unperturbed=nr_genes-nr_perturbed
    for i in range(len(C)):
        for j in range(len(C[i])):
            if kind==2 and MCMC_distr[i][j]<cutoff:
                break
            elif kind==1 and j>=cutoff:
                break
            elif kind==3 and nr_unperturbed_explained[i]*1./nr_unperturbed>cutoff:
                print MCMC_distr[i][j-1]
                break
            for gene in T[C[i][j]]:
                index=levels_ordinal[gene]
                res[i][index]+=1
                if kind==3 and res[i][index]==1 and index>=nr_perturbed:
                    nr_unperturbed_explained[i]+=1
            counter[i]+=1
    #print "aaaaaaaaaaaaa",counter
    width=nr_perturbed*1./n_intervals
    explained=[[0 for i in range(n_intervals+1)] for i in range(len(C))]
    for i in range(len(C)-1,-1,-1):
        if nr_unperturbed_explained[i]*1./nr_unperturbed<=cutoff: #drop this run from mean value
            explained.pop(i)
            continue
        for j in range(n_intervals):
            explained[i][j]=sum([el>0 for el in res[i][int(width*j):int(width*(j+1))]])
        explained[i][n_intervals]=sum([el>0 for el in res[i][nr_perturbed:]])
    #print "bbbbbbbbbbbbbbb",explained

    means=np.mean(np.matrix(explained),0)
    means=np.array([means[0,i] for i in range(n_intervals+1)])
    stds=np.std(np.matrix(explained),0)
    stds=np.array([stds[0,i] for i in range(n_intervals+1)])
    return (means,stds,width)

def nr_annotations(C,T,MCMC_distr,levels_ordinal,nr_perturbed,nr_genes,n_intervals=10,cutoff=0):
    res=[0 for i in range(nr_genes)]
    for t in T:
        for gene in t:
            res[levels_ordinal[gene]]+=1
        
    width=nr_perturbed*1./n_intervals
    explained=[0 for i in range(n_intervals)]
    for j in range(n_intervals):
        explained[j]=sum([el for el in res[int(width*j):int(width*(j+1))]])
    return (explained,width)

def plotte_categories_explained(data_vec,m,min_level=20,legend="",title="",identifier="",MCMC_threshold=0,colors=[]):
    import numpy as np
    import matplotlib.pyplot as plt
    
    if colors==[]:
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
    ax.set_ylabel('Proportion of explained, perturbed genes')
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
    filename=identifier+"_cutoff"+str(m.cutoff)+"_topcutoff"+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+"_pert"+str(min_level)+'_tau'+str(int(round(100*MCMC_threshold)))+'.eps'
    plt.savefig(filename, bbox_inches=0)
    #plt.show()
    return ax

def plotte_categories_explained2(data_vec,m,nr_unperturbed,min_level=20,legend="",title="",identifier="",MCMC_threshold=0, fontsize=12, ymax=1,colors=[]):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    if colors==[]:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    
    N = len(data_vec)
    n_intervals=len(data_vec[0][0])-1
    
    ind = np.arange(n_intervals)  # the x locations for the groups
    width=0.85/N       # the width of the bars
    
    
    rects=[]
    fig, ax = plt.subplots()
    for i in range(N):
        print data_vec[i]
        rects.append(ax.bar(ind+i*width, data_vec[i][0][:-1]/data_vec[i][2], width, yerr=data_vec[i][1][:-1]/data_vec[i][2], error_kw={'ecolor':colors[i], 'capsize':4}, color=colors[i]))
        ax.bar([n_intervals+i*width+0.5], [data_vec[i][0][-1]*1./nr_unperturbed], width, yerr=data_vec[i][1][-1]*1./nr_unperturbed, error_kw={'ecolor':colors[i], 'capsize':4}, color=colors[i])    
    ax.bar([n_intervals-0.075], [1], 0.01, color='k')    

    # add some
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,0,ymax))
    ax.set_ylabel('Proportion of explained genes')
    #ax.set_xlabel('Perturbed genes grouped by rank')
    if title!="":
        ax.set_title(title)
    xticks=np.arange(n_intervals+1)+0.425
    xticks[-1]+=0.5
    ax.set_xticks(xticks)
    xticklabels=[str(i+1) for i in range(n_intervals)]
    xticklabels.append('Unperturbed\nGenes')
    ax.set_xticklabels(xticklabels)
    #ax.get_children()[7+2*n_intervals].set_color('k')
    
    for i in range(N):
        ax.plot(np.arange(n_intervals+1),[np.mean(data_vec[i][0][:-1])/data_vec[i][2]]*(n_intervals+1), color=colors[i], linestyle='--',linewidth=2)
    
    if legend!="":
        ax.legend( rects, legend )
    filename=identifier+"_cutoff"+str(m.cutoff)+"_topcutoff"+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+"_bel"+str(m.belief)+"_pert"+str(min_level)+'_tau'+str(int(round(100*MCMC_threshold)))+'.eps'
    plt.savefig(filename, bbox_inches=0)
    #plt.show()
    return ax

def freq_explained(C,T,MCMC_distr,levels_ordinal,nr_perturbed,nr_genes,cutoff=0,kind=2):
    import numpy as np
    if kind==1:
        C=C[:cutoff]
    res=[[0 for i in range(nr_genes)] for i in range(len(C))]
    counter=[0]*len(C)
    for i in range(len(C)):
        for j in range(len(C[i])):
            if kind==2 and MCMC_distr[i][j]<cutoff:
                break
            for gene in T[C[i][j]]:
                res[i][levels_ordinal[gene]]+=1
            counter[i]+=1
    explained=[[[0 for ii in range(6)] for i in range(2)] for j in range(len(C))]
    for i in range(len(C)):
        for j in range(5):
            explained[i][0][j]=sum([el==j for el in res[i][:nr_perturbed]])
            explained[i][1][j]=sum([el==j for el in res[i][nr_perturbed:]])
        explained[i][0][-1]=sum([el>4 for el in res[i][:nr_perturbed]])
        explained[i][1][-1]=sum([el>4 for el in res[i][nr_perturbed:]])
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
    
def stats(C,T,MCMC_distr,levels_ordinal,nr_perturbed,nr_genes,cutoff=0,kind=2):
    a=jaccard_total_genes(C,T,MCMC_distr,cutoff,kind=kind)
    b=nr_explained(C,T,MCMC_distr,levels_ordinal,nr_perturbed,nr_genes,cutoff,kind)
    c=freq_explained(C,T,MCMC_distr,levels_ordinal,nr_perturbed=nr_perturbed,nr_genes=nr_genes,cutoff=cutoff,kind=kind)
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

 
def plotte4(data_full,data_G_full=[],what=[0],which=1,s=0,perc_active=[10,20,30],N=100,STDS=True,PROPORTIONS=False,pos_leg=0):
    import matplotlib.pyplot as plt
    if type(which)==int:
        which=[which]
    if type(s)==int:
        s=[s]
    if type(perc_active)==int:
        perc_active=[perc_active]
    for val in [10,20,30]:
        if val in perc_active:
            perc_active[perc_active.index(val)]=val/10-1
    data=[]
    data_G=[]
    labels=[]
    title=[]
    if len(which)>1:
        for w in which:
            data.append(data_full[w][s[0]][perc_active[0]])
            data_G.append(data_G_full[w][s[0]][perc_active[0]])            
            if w==0:
                labels.append('HDF Data')
            elif w==1:
                labels.append('3-cell vs 2-cell - GO')
            elif w==2:
                labels.append('3-cell vs 2-cell - MSigDB')
        title_end="different datasets"
    else:
        title.append(['HDF Data','3-cell vs 2-cell - GO','3-cell vs 2-cell - MSigDB'][which[0]])
    if len(which)==1 and len(s)>1:
        for w in s:
            data.append(data_full[which[0]][w][perc_active[0]])
            data_G.append(data_G_full[which[0]][w][perc_active[0]])            
            if w==0:
                labels.append('MCMC')
            elif w==1:
                labels.append('GOing Bayesian')
        title_end="different algorithms"
    else:
        title.append(['MCMC','GOing Bayesian'][s[0]])
    if len(which)==1 and len(s)==1 and len(perc_active)>1:
        for w in perc_active:
            data.append(data_full[which[0]][s[0]][w])    
            data_G.append(data_G_full[which[0]][s[0]][w])            
            labels.append(str((w+1)*10)+'%')
        title_end="different cutoffs"
    else:
        title.append(['10% perturbed genes','20% perturbed genes','30% perturbed genes'][perc_active[0]])
    title=title[0]+', '+title[1]+' and '+title_end
    if STDS==True:
        linestyles=['-','-','-']
    else:
        linestyles=['-','--','-.']
    if PROPORTIONS:
        for i in range(N):
            for j in range(len(data)):
                data[j][1][0][i]*=1./(data_G[j][0]-data_G[j][1])
                data[j][1][1][i]*=1./(data_G[j][0]-data_G[j][1])
                data[j][2][0][i]*=1./data_G[j][1]
                data[j][2][1][i]*=1./data_G[j][1]
        label=['Jaccard Similarity of explained genes for different runs','Proportion explained perturbed genes','Proportion explained unperturbed genes','Number of explaining processes']
    else:
        label=['Jaccard Similarity of explained genes for different runs','Explained perturbed genes','Explained unperturbed genes','Number of explaining processes']
    x=np.linspace(0,1,N)
    #fig, ax1 = plt.subplots()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    if len(data)==2:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    else:
        colors=['b','r','g','m',[0.4,0.4,1],[1,0.4,0.4],[0.4,1,0.4],[1,0.4,1],[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]        
    ind=what[0]
    lines=[]
    for j in range(len(data)):
        lines.append(ax1.plot(x,data[j][ind][0],linestyle=linestyles[j],color=colors[ind+4*j],label=labels[j]))
        if STDS:
            ax1.plot(x,data[j][ind][0]+data[j][ind][1],linestyle='--',color=colors[ind+4*j])
            ax1.plot(x,data[j][ind][0]-data[j][ind][1],linestyle='--',color=colors[ind+4*j])
    ax1.set_xlabel(r"Posterior Probability Cutoff $\tau$")
    ax1.set_ylabel(label[ind], color=colors[ind])
    for tl in ax1.get_yticklabels():
        tl.set_color(colors[ind])
    if ind==0:
        x1,x2,y1,y2 = ax1.axis()
        ax1.axis((x1,x2,0,1))
    else:
        x1,x2,y1,y2 = ax1.axis()
        ax1.axis((x1,x2,0,y2))
    if len(what)>1:
        ind=what[1]
        ax2 = ax1.twinx()
        for j in range(len(data)):
            lines.append(ax2.plot(x,data[j][ind][0],linestyle=linestyles[j],color=colors[ind+4*j],label=labels[j]))
            if STDS:
                ax2.plot(x,data[j][ind][0]+data[j][ind][1],linestyle='--',color=colors[ind+4*j])
                ax2.plot(x,data[j][ind][0]-data[j][ind][1],linestyle='--',color=colors[ind+4*j])
        ax2.set_ylabel(label[ind], color=colors[ind])
        for tl in ax2.get_yticklabels():
            tl.set_color(colors[ind])
        if ind==0:
            x1,x2,y1,y2 = ax2.axis()
            ax2.axis((x1,x2,0,1))
        else:
            x1,x2,y1,y2 = ax2.axis()
            ax2.axis((x1,x2,0,y2))
    if PROPORTIONS:
        for i in range(N):
            for j in range(len(data)):
                data[j][1][0][i]*=(data_G[j][0]-data_G[j][1])
                data[j][1][1][i]*=(data_G[j][0]-data_G[j][1])
                data[j][2][0][i]*=data_G[j][1]
                data[j][2][1][i]*=data_G[j][1]
    filename='aaabd.eps'
    ls=lines[0]
    for i in range(1,len(lines)):
        ls+=lines[i]
    labs = [l[0].get_label() for l in lines]
    ax1.legend(ls,labs,loc=pos_leg)
    plt.title(title)
    plt.savefig(filename, bbox_inches=0)

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
        
def create_table_row(C,T,MCMC_distr,levels_ordinal,nr_perturbed,nr_genes,method_name,cutoff=0,kind=2,alles=0):
    if type(C[0])==int:
        ((mean_JI,std_JI),(mean_pert,std_pert),(mean_unpert,std_unpert),(mean_nr_terms,std_nr_terms),(mean_freq,std_freq))=stats([C],T,MCMC_distr,levels_ordinal,nr_perturbed,nr_genes,cutoff,kind)
        mean_nr_terms=len(C)
    else:
        ((mean_JI,std_JI),(mean_pert,std_pert),(mean_unpert,std_unpert),(mean_nr_terms,std_nr_terms),(mean_freq,std_freq))=stats(C,T,MCMC_distr,levels_ordinal,nr_perturbed,nr_genes,cutoff,kind)
    if alles==True:
        if method_name.lower() in ["gorilla","greedy","fa2.0"]:
            stri=(method_name+" & "+str(int(mean_nr_terms))+" & "+str(rg(mean_pert))+r" \newline "+str(rg(100*mean_pert*1./nr_perturbed))+r"\% & "+str(rg(mean_unpert))+r" \newline "+str(rg(100*mean_unpert*1./(nr_genes-nr_perturbed)))+r"\% \\ \hline") 
        else:
            stri=(method_name+r" \newline ($\tau="+str(cutoff)+"$) & "+str(rg(mean_nr_terms))+r" $\pm$ "+str(rg(std_nr_terms))+" & "+str(rg(mean_pert))+r" $\pm$ "+str(rg(std_pert))+r" \newline "+str(rg(100*mean_pert*1./nr_perturbed))+r"\% $\pm$ "+str(rg(100*std_pert*1./nr_perturbed))+r"\% & "+str(rg(mean_unpert))+r"\% \pm "+str(rg(std_unpert))+r" \newline "+str(rg(100*mean_unpert*1./(nr_genes-nr_perturbed)))+r"\% $\pm$ "+str(rg(100*std_unpert*1./(nr_genes-nr_perturbed)))+r"\% \\ \hline") 
    else:
        if method_name.lower() in ["gorilla","greedy","fa2.0"]:
            stri=(method_name+" & "+str(int(mean_nr_terms))+" & "+str(rg(100*mean_pert*1./nr_perturbed))+r"\% & "+str(rg(100*mean_unpert*1./(nr_genes-nr_perturbed)))+r"\% \\ \hline") 
        else:
            stri=(method_name+" & "+str(rg(mean_nr_terms))+r" $\pm$ "+str(rg(std_nr_terms))+" & "+str(rg(100*mean_pert*1./nr_perturbed))+r"\% $\pm$ "+str(rg(100*std_pert*1./nr_perturbed))+r"\% & "+str(rg(100*mean_unpert*1./(nr_genes-nr_perturbed)))+r"\% $\pm$ "+str(rg(100*std_unpert*1./(nr_genes-nr_perturbed)))+r"\% \\ \hline") 

    return stri
    
def pert_vs_unpert(termss,G,T,legend=""):
    import matplotlib.pyplot as plt
    import numpy as np
    pep=[] #proportion explained perturbed genes
    peu=[] #proportion explained unperturbed genes
    P=set.union(*[set([])] + [G[i] for i in range(len(G)-1)])
    U=set(G[-1])
    sizeP=len(P)
    sizeU=len(U)
    for j,terms in enumerate(termss):
        E=set([])
        pep.append([])
        peu.append([])
        print j, terms[:10]
        for i in range(len(terms)):
            E=set.union(E,set(T[terms[i]]))
            pep[j].append(len(E&P))
            peu[j].append(len(E)-pep[j][-1])
            
    fig, ax = plt.subplots()
    rects=[]
    for j in range(len(termss)):
        rects.append(ax.plot(np.array(peu[j])*1./sizeU,np.array(pep[j])*1./sizeP))
    ax.set_xlabel("Proportion Explained Unperturbed Genes")
    ax.set_ylabel("Proportion Explained Perturbed Genes")
    if legend!="":
        ax.legend(legend,loc=4)
    return (peu,pep)
    
def pert_vs_unpert2(termss,G,T,legend="",means=[],coarsity=0.1,colors=[],opacity=1,show_annotations=True):
    import matplotlib.pyplot as plt
    import numpy as np
    
    if colors==[]:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    
    pep=[] #proportion explained perturbed genes
    peu=[] #proportion explained unperturbed genes
    P=set.union(*[set([])] + [G[i] for i in range(len(G)-1)])
    U=set(G[-1])
    sizeP=len(P)
    sizeU=len(U)
    for j,terms in enumerate(termss):
        E=set([])
        pep.append([])
        peu.append([])
        print j, terms[:10]
        for i in range(len(terms)):
            E=set.union(E,set(T[terms[i]]))
            pep[j].append(len(E&P))
            peu[j].append(len(E)-pep[j][-1])
            
    fig, ax = plt.subplots()
    rects=[]
    for j in range(len(termss)):
        x=np.array(peu[j])*1./sizeU
        y=np.array(pep[j])*1./sizeP
        if type(legend)==list and j%(len(termss)/len(legend))==1: #len(legend)==len(termss)
            rects.append(ax.plot(x,y,color=colors[j],alpha=opacity,label=legend[j/(len(termss)/len(legend))]))
        else:
            rects.append(ax.plot(x,y,color=colors[j],alpha=opacity))
    counter=0
    for j in range(len(termss)):
        x=np.array(peu[j])*1./sizeU
        y=np.array(pep[j])*1./sizeP
        if means[j]!=[]:
            x=list(x)
            print max(x)*1./coarsity,1./coarsity,int(math.floor(max(x)*1./coarsity))
            for jj in range(int(math.floor(max(x)*1./coarsity))):
                dummy=[abs(el-(jj+1)*coarsity) for el in x]
                ind=dummy.index(min(dummy))
                #ax.plot(x,[y[ind]]*len(x),color=colors[j])
                ax.plot([(jj+1)*coarsity,(jj+1)*coarsity,0],[0,y[ind],y[ind]],color=colors[j],linestyle='--')
            for jj in [0.75,0.5,0.1,0.01]:
                ind=sum(np.array(means[j])>jj)
                if ind==0:
                    continue
                print j,jj
                #print ind,x[ind],y[ind],jj
                ax.plot([x[ind]],[y[ind]],color=colors[j],marker='o', ls='')
                if show_annotations:
                    ax.annotate(r"$\tau=$"+str(jj), (x[ind]+0.005,y[ind]-counter*0.03))
            counter+=1
    ax.set_xlabel("Proportion Explained Unperturbed Genes")
    ax.set_ylabel("Proportion Explained Perturbed Genes")
    if legend!="":
        ax.legend(legend,loc=4)
        #handles, labels = ax.get_legend_handles_labels()
        #ax.legend(handles,labels,loc=4)
        #ax.legend([rects[ind] for ind in uniq_colors_indices],legend,loc=4)
    return (peu,pep)
    
def pert_vs_unpert2_ratio(termss,G,T,legend="",means=[],coarsity=0.1,colors=[],opacity=1,show_annotations=True):
    import matplotlib.pyplot as plt
    import numpy as np
    
    if colors==[]:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    
    pep=[] #proportion explained perturbed genes
    peu=[] #proportion explained unperturbed genes
    P=set.union(*[set([])] + [G[i] for i in range(len(G)-1)])
    U=set(G[-1])
    sizeP=len(P)
    sizeU=len(U)
    for j,terms in enumerate(termss):
        E=set([])
        pep.append([])
        peu.append([])
        print j, terms[:10]
        for i in range(len(terms)):
            E=set.union(E,set(T[terms[i]]))
            pep[j].append(len(E&P))
            peu[j].append(len(E)-pep[j][-1])
            
    fig, ax = plt.subplots()
    rects=[]
    for j in range(len(termss)):
        x=np.array(peu[j])*1./sizeU
        y=np.array(pep[j])*1./sizeP
        if type(legend)==list and j%(len(termss)/len(legend))==1: #len(legend)==len(termss)
            rects.append(ax.plot(x,y*1./x,color=colors[j],alpha=opacity,label=legend[j/(len(termss)/len(legend))]))
        else:
            rects.append(ax.plot(x,y*1./x,color=colors[j],alpha=opacity))
    counter=0
    for j in range(len(termss)):
        x=np.array(peu[j])*1./sizeU
        y=np.array(pep[j])*1./sizeP
        if means[j]!=[]:
            x=list(x)
            print max(x)*1./coarsity,1./coarsity,int(math.floor(max(x)*1./coarsity))
            for jj in range(int(math.floor(max(x)*1./coarsity))):
                dummy=[abs(el-(jj+1)*coarsity) for el in x]
                ind=dummy.index(min(dummy))
                #ax.plot(x,[y[ind]]*len(x),color=colors[j])
                ax.plot([(jj+1)*coarsity,(jj+1)*coarsity,0],[0,y[ind]*1/x[ind],y[ind]*1/x[ind]],color=colors[j],linestyle='--')
            for jj in [0.75,0.5,0.1,0.01]:
                try:
                    ind=sum(np.array(means[j])>jj)
                    if ind==0:
                        continue
                    print ind,x[ind],y[ind],jj
                    ax.plot([x[ind]],[y[ind]*1/x[ind]],color=colors[j],marker='o', ls='')
                    if show_annotations:
                        ax.annotate(r"$\tau=$"+str(jj), (x[ind]+0.005,y[ind]*1/x[ind]-counter*0.03))
                except:
                    pass
            counter+=1
    ax.set_xlabel("Proportion Explained Unperturbed Genes")
    ax.set_ylabel("Ratio (Explained Perturbed)/(Explained Unperturbed) Genes")
    #ax.set_ylabel(r"Ratio $\frac{\text{Explained Perturbed Genes}}{\text{Explained Unperturbed Genes}}$")
    if legend!="":
        ax.legend(legend,loc=0)
    return (peu,pep)

def pert_vs_unpert2a_ratio(termss,G,T,legend="",means=[],coarsity=0.1,colors=[],opacity=1,highlights=[0.75,0.5,0.1,0.01],show_annotations=True,max_x=1,fontsize=16):
    import matplotlib.pyplot as plt
    import numpy as np

    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : fontsize}

    matplotlib.rc('font', **font)    
            
    if colors==[]:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    
    pep=[] #proportion explained perturbed genes
    peu=[] #proportion explained unperturbed genes
    P=set.union(*[set([])] + [G[i] for i in range(len(G)-1)])
    U=set(G[-1])
    sizeP=len(P)
    sizeU=len(U)
    for j,terms in enumerate(termss):
        E=set([])
        pep.append([])
        peu.append([])
        print j, terms[:10]
        for i in range(len(terms)):
            E=set.union(E,set(T[terms[i]]))
            pep[j].append(len(E&P))
            peu[j].append(len(E)-pep[j][-1])
            
    fig, ax = plt.subplots()
    rects=[]
    for j in range(len(termss)):
        x=np.array(peu[j])*1./sizeU
        y=np.array(pep[j])*1./sizeP
        if type(legend)==list and j%(len(termss)/len(legend))==1: #len(legend)==len(termss)
            rects.append(ax.plot(x,y*1./x,linewidth=2,color=colors[j],alpha=opacity,label=legend[j/(len(termss)/len(legend))])[0])
        else:
            rects.append(ax.plot(x,y*1./x,linewidth=2,color=colors[j],alpha=opacity)[0])
    counter=0
    symbols=[]
    markers=['o','d','s','^']
    for j in range(len(termss)):
        x=np.array(peu[j])*1./sizeU
        y=np.array(pep[j])*1./sizeP
        if means[j]!=[]:
            x=list(x)
            print max(x)*1./coarsity,1./coarsity,int(math.floor(max(x)*1./coarsity))
            for jj in range(int(math.floor(max(x)*1./coarsity))):
                dummy=[abs(el-(jj+1)*coarsity) for el in x]
                ind=dummy.index(min(dummy))
                #ax.plot(x,[y[ind]]*len(x),color=colors[j])
                ax.plot([(jj+1)*coarsity,(jj+1)*coarsity,0],[0,y[ind]*1/x[ind],y[ind]*1/x[ind]],color=colors[j],linestyle='--')
            count=-1
            for jj in highlights:
                count+=1
                try:
                    ind=sum(np.array(means[j])>jj)
                    if ind==0:
                        continue
                    print ind,x[ind],y[ind],jj
                    if j==0:
                        symbols.append(ax.plot([x[ind]],[y[ind]*1/x[ind]],color='k',marker=markers[count], ms=8,ls='')[0])
                    ax.plot([x[ind]],[y[ind]*1/x[ind]],color=colors[j],marker=markers[count], ms=8,ls='')
                    if show_annotations:
                        ax.annotate(r"$\tau=$"+str(jj), (x[ind]+0.005,y[ind]*1/x[ind]-counter*0.03))
                except:
                    pass
            counter+=1
    ax.set_xlabel("Proportion Explained Unperturbed Genes")
    ax.set_ylabel("Quality")# (Explained Perturbed)/(Explained Unperturbed) Genes")
    #ax.set_ylabel(r"Ratio $\frac{\text{Explained Perturbed Genes}}{\text{Explained Unperturbed Genes}}$")
    if max_x<1:
        x1,x2,y1,y2 = plt.axis()
        plt.axis((0,max_x,1,y2))
    if legend!="":
        l1 = ax.legend(rects,legend,loc=0)
        print symbols,[r"$\tau=$"+str(el) for el in highlights]
        l2 = ax.legend(symbols,[r"$\tau=$"+str(el) for el in highlights],loc=5,numpoints = 1)
        gca().add_artist(l1)
    return (peu,pep)
    
def pert_vs_unpert2b_ratio(termss,G,T,top_proportion=0,bottom_proportion=1,legend="",means=[],coarsity=0.1,colors=[],opacity=1,highlights=[0.75,0.5,0.1,0.01],show_annotations=True,max_x=1,title="",identifier="",fontsize=16,PLOT=True,SAVE=False,m=[],min_level=30,pert_or_unpert_on_xaxis='unpert',second_legend=False,show_legend=True):
    import matplotlib.pyplot as plt
    import numpy as np

    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : fontsize}

    matplotlib.rc('font', **font)    
            
    if colors==[]:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    
    pep=[] #proportion explained perturbed genes
    peu=[] #proportion explained unperturbed genes
    P=set.union(*[set([])] + [G[i] for i in range(int((len(G)-1)*top_proportion),int((len(G)-1)*bottom_proportion))])
    U=set(G[-1])
    sizeP=len(P)
    sizeU=len(U)
    for j,terms in enumerate(termss):
        E=set([])
        pep.append([])
        peu.append([])
        print j, terms[:10]
        for i in range(len(terms)):
            E=set.union(E,set(T[terms[i]]))
            pep[j].append(len(E&P))
            peu[j].append(len(E&U))
    
    if PLOT==True:
        fig, ax = plt.subplots()
        rects=[]
        for j in range(len(termss)):
            x=np.array(peu[j])*1./sizeU
            y=np.array(pep[j])*1./sizeP
            if pert_or_unpert_on_xaxis=='unpert':
                x_axis=x.copy()
            else:
                x_axis=y.copy()
            if type(legend)==list and j%(len(termss)/len(legend))==1: #len(legend)==len(termss)
                rects.append(ax.plot(x_axis,y*1./x,linewidth=2,color=colors[j],alpha=opacity,label=legend[j/(len(termss)/len(legend))])[0])
            else:
                rects.append(ax.plot(x_axis,y*1./x,linewidth=2,color=colors[j],alpha=opacity)[0])
        counter=0
        symbols=[]
        markers=['o','d','s','^']
        for j in range(len(termss)):
            x=np.array(peu[j])*1./sizeU
            y=np.array(pep[j])*1./sizeP
            if pert_or_unpert_on_xaxis=='unpert':
                x_axis=x.copy()
            else:
                x_axis=y.copy()
            if means[j]!=[]:
                x_axis=list(x_axis)
                print max(x),max(x)*1./coarsity,1./coarsity,int(math.floor(max(x)*1./coarsity))
                for jj in range(int(math.floor(max(x)*1./coarsity))):
                    dummy=[abs(el-(jj+1)*coarsity) for el in x_axis]
                    ind=dummy.index(min(dummy))
                    #ax.plot(x,[y[ind]]*len(x),color=colors[j])
                    ax.plot([(jj+1)*coarsity,(jj+1)*coarsity,0],[0,y[ind]*1/x[ind],y[ind]*1/x[ind]],color=colors[j],linestyle='--')
                count=-1
                for jj in highlights:
                    count+=1
                    try:
                        ind=sum(np.array(means[j])>jj)
                        if ind==0:
                            continue
                        print ind,x[ind],y[ind],jj
                        if j==0:
                            symbols.append(ax.plot([x_axis[ind]],[y[ind]*1/x[ind]],color='k',marker=markers[count], ms=8,ls='')[0])
                        ax.plot([x_axis[ind]],[y[ind]*1/x[ind]],color=colors[j],marker=markers[count], ms=8,ls='')
                        if show_annotations:
                            ax.annotate(r"$\tau=$"+str(jj), (x_axis[ind]+0.005,y[ind]*1/x[ind]-counter*0.03))
                    except:
                        pass
                counter+=1
        if pert_or_unpert_on_xaxis=='unpert':
            ax.set_xlabel("Proportion of Unperturbed Genes that are Explained")
        else:
            ax.set_xlabel("Proportion of Perturbed Genes that are Explained")            
        ylabel="Quality"
        filename_add=""
        if 0<bottom_proportion<1 and top_proportion==0:
            #ylabel+=" of Top "+str(int(bottom_proportion*100))+"% Perturbed Genes"
            filename_add="_top"+str(int(bottom_proportion*100))
            if pert_or_unpert_on_xaxis!='unpert':
                ax.set_xlabel("Proportion of Top "+str(int(bottom_proportion*100))+"% Perturbed Genes that are Explained")    
        elif bottom_proportion==1 and top_proportion!=0:
            #ylabel+=" of Lowest "+str(int(top_proportion*100))+"% Perturbed Genes"
            filename_add="_lowest"+str(int(top_proportion*100))
            if pert_or_unpert_on_xaxis!='unpert':
                ax.set_xlabel("Proportion of Lowest "+str(int(top_proportion*100))+"% Perturbed Genes that are Explained")
        elif bottom_proportion!=0 and top_proportion!=0:
            #ylabel+=" of "+str(int(top_proportion*100))+" - "+str(int(bottom_proportion*100))+"% Perturbed Genes"
            filename_add="_"+str(int(bottom_proportion*100))+"to"+str(int(top_proportion*100))
            if pert_or_unpert_on_xaxis!='unpert':
                ax.set_xlabel("Proportion of "+str(int(top_proportion*100))+" - "+str(int(bottom_proportion*100))+"% Perturbed Genes that are Explained")
            
        ax.set_ylabel(ylabel)# (Explained Perturbed)/(Explained Unperturbed) Genes")
        #ax.set_ylabel(r"Ratio $\frac{\text{Explained Perturbed Genes}}{\text{Explained Unperturbed Genes}}$")
        x1,x2,y1,y2 = plt.axis()
        plt.axis((0,min(max_x,x2),1,y2))
        if legend!="":
            l1 = ax.legend(rects,legend,loc=0)
            print symbols,[r"$\tau=$"+str(el) for el in highlights]
            if second_legend==True:
                l2 = ax.legend(symbols,[r"$\tau=$"+str(el) for el in highlights],loc=5,numpoints = 1)
                gca().add_artist(l1)
        plt.title(title)
        if SAVE:
            filename="quality_"+identifier
            filename+="_cutoff"+str(m.cutoff)+"_topcutoff"+str(m.top_cutoff)+"_pert"+str(min_level)
            filename+=filename_add+'.eps'
            plt.savefig(filename, bbox_inches=0)
    return (peu,pep)

def pert_vs_unpert2bb_ratio(termss,G,T,top_proportions=[0,0,0],bottom_proportions=[0.25,0.5,1],legend="",means=[],coarsity=0.1,colors=[],opacity=1,highlights=[0.75,0.5,0.1,0.01],show_annotations=True,max_x=1,title="",identifier="",fontsize=16,PLOT=True,SAVE=False,m=[],min_level=30,pert_or_unpert_on_xaxis='unpert',second_legend=False):
    #Legend only perfectlyworks for len(termss)==1, sufficient for paper figures
    import matplotlib.pyplot as plt
    import numpy as np

    print len(termss),legend

    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : fontsize}

    matplotlib.rc('font', **font)    
            
    if colors==[]:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    linestyles=['--','-',':','-.']
    pep=[] #proportion explained perturbed genes
    peu=[] #proportion explained unperturbed genes
    sizePs=[]
    for top_proportion,bottom_proportion in zip(top_proportions,bottom_proportions):
        P=set.union(*[set([])] + [G[i] for i in range(int((len(G)-1)*top_proportion),int((len(G)-1)*bottom_proportion))])
        U=set(G[-1])
        sizePs.append(len(P))
        sizeU=len(U)
        pep.append([])
        peu.append([])
        for j,terms in enumerate(termss):
            E=set([])
            pep[-1].append([])
            peu[-1].append([])
            print j, terms[:10]
            for i in range(len(terms)):
                E=set.union(E,set(T[terms[i]]))
                pep[-1][j].append(len(E&P))
                peu[-1][j].append(len(E&U))
    
    if PLOT==True:
        fig, ax = plt.subplots()
        rects=[]
        for i in range(len(top_proportions)):
            for j in range(len(termss)):
                x=np.array(peu[i][j])*1./sizeU
                y=np.array(pep[i][j])*1./sizePs[i]
                if pert_or_unpert_on_xaxis=='unpert':
                    x_axis=x.copy()
                else:
                    x_axis=y.copy()
                rects.append(ax.plot(x_axis,y*1./x,linewidth=2,color=colors[j],linestyle=linestyles[i],alpha=opacity,label=r"$p=$"+str(int(round(bottom_proportions[i]*100))))[0])
            counter=0
            symbols=[]
            markers=['o','d','s','^']
            for j in range(len(termss)):
                x=np.array(peu[i][j])*1./sizeU
                y=np.array(pep[i][j])*1./sizePs[i]
                if pert_or_unpert_on_xaxis=='unpert':
                    x_axis=x.copy()
                else:
                    x_axis=y.copy()
                if means[j]!=[]:
                    x_axis=list(x_axis)
                    print max(x),max(x)*1./coarsity,1./coarsity,int(math.floor(max(x)*1./coarsity))
                    for jj in range(int(math.floor(max(x)*1./coarsity))):
                        dummy=[abs(el-(jj+1)*coarsity) for el in x_axis]
                        ind=dummy.index(min(dummy))
                        #ax.plot(x,[y[ind]]*len(x),color=colors[j])
                        ax.plot([(jj+1)*coarsity,(jj+1)*coarsity,0],[0,y[ind]*1/x[ind],y[ind]*1/x[ind]],color=colors[j],linestyle='--')
                    count=-1
                    for jj in highlights:
                        count+=1
                        try:
                            ind=sum(np.array(means[j])>jj)
                            if ind==0:
                                continue
                            print ind,x[ind],y[ind],jj
                            if j==0:
                                symbols.append(ax.plot([x_axis[ind]],[y[ind]*1/x[ind]],color='k',marker=markers[count], ms=8,ls='')[0])
                            ax.plot([x_axis[ind]],[y[ind]*1/x[ind]],color=colors[j],marker=markers[count], ms=8,ls='')
                            if show_annotations:
                                ax.annotate(r"$\tau=$"+str(jj), (x_axis[ind]+0.005,y[ind]*1/x[ind]-counter*0.03))
                        except:
                            pass
                    counter+=1
        ax.set_xlabel(r"Proportion of Top $p\%$ Perturbed Genes that are Explained")            
            
        ax.set_ylabel("Quality")# (Explained Perturbed)/(Explained Unperturbed) Genes")
        #ax.set_ylabel(r"Ratio $\frac{\text{Explained Perturbed Genes}}{\text{Explained Unperturbed Genes}}$")
        x1,x2,y1,y2 = plt.axis()
        plt.axis((0,min(max_x,1),1,y2))
        if legend!="":
            l1 = ax.legend(title=legend[0],loc=0)
            print symbols,[r"$\tau=$"+str(el) for el in highlights]
            if second_legend==True:
                l2 = ax.legend(symbols,[r"$\tau=$"+str(el) for el in highlights],loc=5,numpoints = 1)
                gca().add_artist(l1)
        plt.title(title)
        if SAVE:
            filename="quality_"+identifier
            filename+="_cutoff"+str(m.cutoff)+"_topcutoff"+str(m.top_cutoff)+"_pert"+str(min_level)+'.eps'
            plt.savefig(filename, bbox_inches=0)
    return (peu,pep)

def pert_vs_unpert2b_ratio_separate(Cs,G,T,top_proportion=0,bottom_proportion=1,legend="",MCMC_distrs=[],coarsity=0.1,colors=[],opacity=1,highlights=[0.75,0.5,0.1,0.01],show_annotations=True,max_x=1,fontsize=16,PLOT=True):
    peus,peps,x_domain,xs,ys,ns=[],[],[],[],[]
    for i in range(len(Cs)):
        nsim=len(Cs[i])
        (peu,pep)=gr.pert_vs_unpert2b_ratio(Cs[i],G,T,means=MCMC_distrs[i],PLOT=False)
        peus.append(peu)
        peps.append(pep)
        min_min_xval=min([min(el) for el in peu])
        max_min_xval=max([min(el) for el in peu])
        min_max_xval=min([max(el) for el in peu])
        max_max_xval=max([max(el) for el in peu])
        x_domain.append(range(max_max_xval+1))
        ns.append(len(peu)*np.ones(max_max_xval+1,dtype=int))
        ys.append([])
        for j in range(nsim):
            for ii in range(peu[j][0]):
                ns[i][ii]-=1
            for ii in range(peu[j][-1],max_max_xval+1):
                ns[i][ii]-=1
            ys[i].append([-1 for indind in range(max_max_xval+1)])
            for ii in range(len(peu[j])):
                ys[i][j][peu[j][ii]]=pep[j][ii]
        to_be_removed=[]
        for ii in range(max_max_xval+1):
            if sum([ys[i][j][ii]==-1 for j in range(nsim)])==nsim:
                to_be_removed.append(ii)
        xs.append(range(max_max_xval+1))
        for ind in to_be_removed:
            pass
            #else:
            #    for j in range(nsim):
            #        try:
            #            ys[i][j][ii]=
                    
    return (peps,peus,x_domain,ys,ns)
    
def combine_pep_lines(Cs,G,T,top_proportion=0,bottom_proportion=1,MCMC_distrs=[]):
    peus,peps,xs,ys,ns,dict_xs=[],[],[],[],[],[]
    for i in range(len(Cs)):
        count=0
        dict_xs.append({})
        nsim=len(Cs[i])
        (peu,pep)=pert_vs_unpert2b_ratio(Cs[i],G,T,means=MCMC_distrs[i],top_proportion=top_proportion,bottom_proportion=bottom_proportion,PLOT=False)
        peus.append(peu)
        peps.append(pep)
        min_min_xval=min([min(el) for el in peu])
        max_min_xval=max([min(el) for el in peu])
        min_max_xval=min([max(el) for el in peu])
        max_max_xval=max([max(el) for el in peu])
        ys_pre=[[-1 for ii in range(max_max_xval)] for j in range(nsim)]
        for j in range(nsim):
            for ii in range(len(peu[j])):
                try:
                    ind=dict_xs[i][peu[j][ii]]
                except KeyError:
                    dict_xs[i].update({peu[j][ii]:count})
                    ind=count
                    count+=1
                ys_pre[j][ind]=pep[j][ii]
        xs.append(dict_xs[i].keys())
        ys.append([])
        ns.append([0 for ii in range(count)])
        to_estimate=[]
        for j in range(nsim):
            ys[i].append([])
            for ii,x in enumerate(xs[i]):
                if x<peu[j][0]:
                    ys[i][j].append(-1)
                elif x>peu[j][-1]:
                    ys[i][j].append(-1)
                    break
                else:
                    ind=dict_xs[i][x]
                    ns[i][ii]+=1
                    if ys_pre[j][ind]!=-1:
                        next_val_x=x
                        next_val_y=ys_pre[j][ind]
                        #estimate missing values
                        for (iiii,xx) in to_estimate:
                            ys[i][j][iiii]=last_val_y+(next_val_y-last_val_y)*(xx-last_val_x)*1./(next_val_x-last_val_x)
                        to_estimate=[]
                        ys[i][j].append(ys_pre[j][ind])
                        last_val_x=x
                        last_val_y=ys_pre[j][ind]
                    else: #need to estimate
                        to_estimate.append((ii,x))
                        ys[i][j].append(-1)
    return (peps,peus,xs,ys,ns)

def pert_vs_unpert2b_ratio_separate2(xs,ys,peus,G,T,top_proportion=0,bottom_proportion=1,legend="",MCMC_distrs=[],coarsity=0.1,colors=[],opacity=1,highlights=[0.75,0.5,0.1,0.01],max_x=1,title="",identifier="",fontsize=16,PLOT=True,SAVE=False,m=[],min_level=30):    
    import matplotlib.pyplot as plt
    import numpy as np

    #Analyze
    x,y,y_mean,y_std=[],[],[],[]
    fig, ax = plt.subplots()
    rects=[]
    P=set.union(*[set([])] + [G[i] for i in range(int((len(G)-1)*top_proportion),int((len(G)-1)*bottom_proportion))])
    U=set(G[-1])
    sizeP=len(P)
    sizeU=len(U)
    for i in range(len(xs)):
        x_min=xs[i].index(max([min(el) for el in peus[i]]))
        x_max=xs[i].index(min([max(el) for el in peus[i]]))
        x.append(np.array(xs[i])[x_min:x_max+1])
        y.append(np.matrix(ys[i])[:,x_min:x_max+1])
        y_mean.append(np.asarray(np.mean(y[i],0)).reshape(-1))
        y_std.append(np.asarray(np.std(y[i],0)).reshape(-1))
        if type(legend)==list: #len(legend)==len(termss)
            rects.append(ax.plot(x[i]*1./sizeU,y_mean[i]*1./sizeP*sizeU/x[i],linewidth=2,color=colors[i],alpha=opacity,label=legend[i])[0])
        else:
            rects.append(ax.plot(x[i]*1./sizeU,y_mean[i]*1./sizeP*sizeU/x[i],linewidth=2,color=colors[i],alpha=opacity)[0])        
        symbols=[]
        markers=['o','d','s','^']
    for i in range(len(xs)):
        if MCMC_distrs!=[] and MCMC_distrs[i][0]!=[]:
            count=-1
            for jj in highlights:
                count+=1
                av_x=0
                for j in range(len(ys[i])):
                    try:
                        ind=sum(np.array(MCMC_distrs[i][j])>jj)
                        #print i,jj,ind
                        if ind==0:
                            continue
                        else:
                            av_x+=xs[i][ind]
                    except:
                        pass
                av_x=av_x*1./len(ys[i])
                ind_av_x=sum(np.array(xs[i])<av_x)
                print i,jj,av_x,ind_av_x
                if ind_av_x==0:
                    continue
                av_y=y_mean[i][ind_av_x]+(y_mean[i][ind_av_x+1]-y_mean[i][ind_av_x])*(av_x-xs[i][ind_av_x])*1./(xs[i][ind_av_x+1]-xs[i][ind_av_x])
                
                if i==0:
                    symbols.append(ax.plot([av_x*1./sizeU],[av_y*1./sizeP*sizeU/av_x],color='k',marker=markers[count], ms=8,ls='')[0])
                ax.plot([av_x*1./sizeU],[av_y*1./sizeP*sizeU/av_x],color=colors[i],marker=markers[count], ms=8,ls='')
    ax.set_xlabel("Proportion Explained Unperturbed Genes")
    ylabel="Quality"
    filename_add=""
    if 0<bottom_proportion<1 and top_proportion==0:
        ylabel+=" of Top "+str(int(bottom_proportion*100))+"% Perturbed Genes"
        filename_add="_top"+str(int(bottom_proportion*100))
    elif bottom_proportion==1 and top_proportion!=0:
        ylabel+=" of Lowest "+str(int(top_proportion*100))+"% Perturbed Genes"
        filename_add="_lowest"+str(int(top_proportion*100))
    elif bottom_proportion!=0 and top_proportion!=0:
        ylabel+=" of "+str(int(top_proportion*100))+" - "+str(int(bottom_proportion*100))+"% Perturbed Genes"
        filename_add="_"+str(int(bottom_proportion*100))+"to"+str(int(top_proportion*100))
            
    ax.set_ylabel(ylabel)# (Explained Perturbed)/(Explained Unperturbed) Genes")
    #ax.set_ylabel(r"Ratio $\frac{\text{Explained Perturbed Genes}}{\text{Explained Unperturbed Genes}}$")
    x1,x2,y1,y2 = plt.axis()
    plt.axis((0,min(max_x,x2),1,y2))
    if legend!="":
        l1 = ax.legend(rects,legend,loc=0)
        if MCMC_distrs!=[]:
            print symbols,[r"$\tau=$"+str(el) for el in highlights]
            l2 = ax.legend(symbols,[r"$\tau=$"+str(el) for el in highlights],loc=5,numpoints = 1)
            gca().add_artist(l1)  
    plt.title(title)     
    if SAVE:
        filename="quality_separate_"+identifier
        filename+="_cutoff"+str(m.cutoff)+"_topcutoff"+str(m.top_cutoff)+"_pert"+str(min_level)
        filename+=filename_add+'.eps'
        plt.savefig(filename, bbox_inches=0)                             
    return (y_mean,y_std)

def pert_vs_unpert2b_ratio_separate3(xs,ys,peus,G,T,top_proportion=0,bottom_proportion=1,legend="",MCMC_distrs=[],coarsity=0.1,colors=[],opacity=1,highlights=[0.75,0.5,0.1,0.01],show_annotations=True,max_x=1,title="",identifier="",fontsize=16,PLOT=True,SAVE=False,m=[]):    
    import matplotlib.pyplot as plt
    import numpy as np

    #Analyze
    x,y,y_mean,y_std=[],[],[],[]
    fig, ax = plt.subplots()
    rects=[]
    P=set.union(*[set([])] + [G[i] for i in range(int((len(G)-1)*top_proportion),int((len(G)-1)*bottom_proportion))])
    U=set(G[-1])
    sizeP=len(P)
    sizeU=len(U)
    for i in range(len(xs)):
        x_min=xs[i].index(max([min(el) for el in peus[i]]))
        x_max=xs[i].index(min([max(el) for el in peus[i]]))
        x.append(np.array(xs[i])[x_min:x_max+1])
        y.append(np.matrix(ys[i])[:,x_min:x_max+1])
        y_mean.append(np.asarray(np.mean(y[i],0)).reshape(-1))
        y_std.append(np.asarray(np.std(y[i],0)).reshape(-1))
        if type(legend)==list: #len(legend)==len(termss)
            rects.append(ax.plot(y_mean[i]*1./sizeP,y_mean[i]*1./sizeP*sizeU/x[i],linewidth=2,color=colors[i],alpha=opacity,label=legend[i/(len(termss)/len(legend))])[0])
        else:
            rects.append(ax.plot(y_mean[i]*1./sizeP,y_mean[i]*1./sizeP*sizeU/x[i],linewidth=2,color=colors[i],alpha=opacity)[0])        
        symbols=[]
        markers=['o','d','s','^']
    for i in range(len(xs)):
        if MCMC_distrs!=[] and MCMC_distrs[i][0]!=[]:
            count=-1
            for jj in highlights:
                count+=1
                av_x=0
                for j in range(len(ys[i])):
                    try:
                        ind=sum(np.array(MCMC_distrs[i][j])>jj)
                        if ind==0:
                            continue
                        else:
                            av_x+=xs[i][ind]
                    except:
                        pass
                av_x=av_x*1./len(ys[i])
                ind_av_x=sum(np.array(xs[i])<av_x)
                if ind_av_x==0:
                    continue
                av_y=y_mean[i][ind_av_x]+(y_mean[i][ind_av_x+1]-y_mean[i][ind_av_x])*(av_x-xs[i][ind_av_x])*1./(xs[i][ind_av_x+1]-xs[i][ind_av_x])
                print i,jj,av_x,av_y
                if i==0:
                    symbols.append(ax.plot([av_y*1./sizeP],[av_y*1./sizeP*sizeU/av_x],color='k',marker=markers[count], ms=8,ls='')[0])
                ax.plot([av_y*1./sizeP],[av_y*1./sizeP*sizeU/av_x],color=colors[i],marker=markers[count], ms=8,ls='')
    ax.set_xlabel("Proportion Explained Perturbed Genes")
    ylabel="Quality"
    if 0<bottom_proportion<1 and top_proportion==0:
        ylabel+=" of Top "+str(int(bottom_proportion*100))+"% Perturbed Genes"
    elif bottom_proportion==1 and top_proportion!=0:
        ylabel+=" of Lowest "+str(int(top_proportion*100))+"% Perturbed Genes"
    elif bottom_proportion!=0 and top_proportion!=0:
        ylabel+=" of "+str(int(top_proportion*100))+" - "+str(int(bottom_proportion*100))+"% Perturbed Genes"
        
    ax.set_ylabel(ylabel)# (Explained Perturbed)/(Explained Unperturbed) Genes")
    #ax.set_ylabel(r"Ratio $\frac{\text{Explained Perturbed Genes}}{\text{Explained Unperturbed Genes}}$")
    x1,x2,y1,y2 = plt.axis()
    plt.axis((0,min(max_x,x2),1,y2))
    if legend!="":
        l1 = ax.legend(rects,legend,loc=0)
        print symbols,[r"$\tau=$"+str(el) for el in highlights]
        l2 = ax.legend(symbols,[r"$\tau=$"+str(el) for el in highlights],loc=5,numpoints = 1)
        gca().add_artist(l1)  
    plt.title(title)     
                             
    return (y_mean,y_std)

def welchs_t_test(xs,ys,ns,ind1,ind2):
    import math
    import numpy as np

    #performs a pointwise between datasets ind1 and ind2 (usually in 0,1,2,3)
    t=[]
    v=[]
    y_mean,y_std=[],[]
    for i in [ind1,ind2]:
        ymat=np.matrix(ys[i])
        x=xs[i]
        y_mean.append(np.asarray(np.mean(ymat,0)).reshape(-1))
        y_std.append(np.asarray(np.std(ymat,0)).reshape(-1))
    for i in range(len(xs[0])):
        t.append((y_mean[0][i]-y_mean[1][i])*1./math.sqrt(y_std[0][i]**2*1./ns[ind1][i]+y_std[1][i]**2*1./ns[ind2][i]))
        v.append((y_std[0][i]**2*1./ns[ind1][i]+y_std[1][i]**2*1./ns[ind2][i])**2/(y_std[0][i]**4*1./ns[ind1][i]**2/(ns[ind1][i]-1)+y_std[1][i]**4*1./ns[ind2][i]**2/(ns[ind2][i]-1)))
    return (t,v)
    
def welchs_t_test2(xs,ys,ns,ind1,ind2):
    import math
    import numpy as np
    from scipy.stats import t
    
    #performs a pointwise between datasets ind1 and ind2 (usually in 0,1,2,3)
    t_val=[]
    df=[]
    p=[]
    y_mean,y_std,dict_xs=[],[],[]
    for i in [ind1,ind2]:
        dict_xs.append({})
        for j,x in enumerate(xs[i]):
            dict_xs[-1].update({x:j})
        ymat=np.matrix(ys[i])
        y_mean.append(np.asarray(np.mean(ymat,0)).reshape(-1))
        y_std.append(np.asarray(np.std(ymat,0)).reshape(-1))
    x=list(set(xs[ind1])&set(xs[ind2]))
    for el in x:
        pos1=dict_xs[0][el]
        pos2=dict_xs[1][el]
        t_val.append((y_mean[0][pos1]-y_mean[1][pos2])*1./math.sqrt(y_std[0][pos1]**2*1./ns[ind1][pos1]+y_std[1][pos2]**2*1./ns[ind2][pos2]))
        df.append((y_std[0][pos1]**2*1./ns[ind1][pos1]+y_std[1][pos2]**2*1./ns[ind2][pos2])**2/(y_std[0][pos1]**4*1./ns[ind1][pos1]**2/(ns[ind1][pos1]-1)+y_std[1][pos2]**4*1./ns[ind2][pos2]**2/(ns[ind2][pos2]-1)))
        p.append(1-t.cdf(t_val[-1],df[-1]))
    return (x,t_val,df,p)

def pert_vs_unpert3(Cs,G,T,legend="",means=[],coarsity=0.1):
    import matplotlib.pyplot as plt
    import numpy as np
    
    def find_nearest_index(array,value):
        if type(array)==set:
            array=np.array(list(array))
        elif type(array)==list:
            array=np.array(array)
        ind = (np.abs(array-value)).argmin()
        return ind
    
    colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    
    peps=[] #proportion explained perturbed genes
    peus=[] #proportion explained unperturbed genes
    P=set.union(*[set([])] + [G[i] for i in range(len(G)-1)])
    U=set(G[-1])
    sizeP=len(P)
    sizeU=len(U)
    for jj,termss in enumerate(Cs):
        peps.append([])
        peus.append([])
        for j,terms in enumerate(termss):
            E=set([])
            peps[jj].append([0])
            peus[jj].append([0])
            for i in range(len(terms)):
                E=set.union(E,set(T[terms[i]]))
                peps[jj][j].append(len(E&P))
                peus[jj][j].append(len(E)-peps[jj][j][-1])
    pep=[]
    peu=[]
    xs=[]
    ys=[]
    count=[]
    for jj,termss in enumerate(Cs):
        xs.append(set.union(*[set([])] + [set(el) for el in peus[jj]]))
        ys.append([])
        pep.append([])
        peu.append([])
        count.append([])
        for n_x,x in enumerate(xs[jj]):
            y=0
            count[jj].append(0)
            for j in range(len(peps[jj])):
                ind=find_nearest_index(peus[jj][j],x)
                count[jj][n_x]+=1
                if peus[jj][j][ind]==x:
                    y+=peps[jj][j][ind]
                elif peus[jj][j][ind]>x:
                    if peus[jj][j][ind]-peus[jj][j][ind-1]>0:
                        y+=((x-peus[jj][j][ind-1])*peps[jj][j][ind]+(peus[jj][j][ind]-x)*peps[jj][j][ind-1])/(peus[jj][j][ind]-peus[jj][j][ind-1])
                    else:
                        y+=(peps[jj][j][ind]+peps[jj][j][ind-1])/2.
                elif peus[jj][j][ind]<x and ind<len(peus[jj][j])-1:
                    if peus[jj][j][ind+1]-peus[jj][j][ind]>0:
                        y+=((x-peus[jj][j][ind])*peps[jj][j][ind+1]+(peus[jj][j][ind+1]-x)*peps[jj][j][ind])/(peus[jj][j][ind+1]-peus[jj][j][ind])
                    else:
                        y+=(peps[jj][j][ind+1]+peps[jj][j][ind])/2.
                else: #if peus[jj][j][ind]<x and ind=len(peus[jj][j])-1
                    count[jj][n_x]-=1
            ys[jj].append(y*1./count[jj][n_x])
    fig, ax = plt.subplots()
    rects=[]
    for j in range(len(Cs)):
        x=np.array(list(xs[jj]))*1./sizeP
        y=np.array(list(ys[jj]))*1./sizeU
        rects.append(ax.plot(x,y,color=colors[j]))
    counter=0
    for j in range(len(Cs)):
        x=np.array(peu[j])*1./sizeU
        y=np.array(pep[j])*1./sizeP
        if means[j]!=[]:
            x=list(x)
            print max(x)*1./coarsity,1./coarsity,int(math.floor(max(x)*1./coarsity))
            for jj in range(int(math.floor(max(x)*1./coarsity))):
                dummy=[abs(el-(jj+1)*coarsity) for el in x]
                ind=dummy.index(min(dummy))
                #ax.plot(x,[y[ind]]*len(x),color=colors[j])
                ax.plot([(jj+1)*coarsity,(jj+1)*coarsity,0],[0,y[ind],y[ind]],color=colors[j],linestyle='--')
            for jj in [0.75,0.5,0.1,0.01]:
                ind=sum(np.array(means[j])>jj)
                if ind==0:
                    continue
                print ind,x[ind],y[ind],jj
                ax.plot([x[ind]],[y[ind]],color=colors[j],marker='o', ls='')
                ax.annotate(r"$\tau=$"+str(jj), (x[ind]+0.005,y[ind]-counter*0.03))
            counter+=1
    ax.set_xlabel("Proportion Explained Unperturbed Genes")
    ax.set_ylabel("Proportion Explained Perturbed Genes")
    if legend!="":
        ax.legend(legend,loc=4)
    return (peu,pep)
    
def jaccard_similarity_tau(termss,G,T,legend="",meanss=[],max_nr_terms=-1,colors=[],opacity=1,show_annotations=True):
    import matplotlib.pyplot as plt
    import numpy as np
    
    if colors==[]:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]

    if max_nr_terms>=0:
        for i in range(len(termss)):
            termss[i]=termss[i][:max_nr_terms]
            if meanss!=[]:
                meanss[i]=meanss[i][:max_nr_terms]

    if meanss==[]:
        invert_xaxis=False
        for terms in termss:
            meanss.append(range(1,1+len(terms)))
        xlabel="Number of Terms in Explanatory Set"
    else:
        invert_xaxis=True
        xlabel="MCMC Posterior Threshold"

    sim=[]
    taus=[]
    for j,terms in enumerate(termss):
        sim.append([0])
        taus.append([meanss[j][0]])
        max_jaccard=[0]
        for i in range(1,len(terms)):
            max_jaccard.append(0)
            for ii,term in enumerate(terms[:i]):
                val=len(set(T[terms[i]])&set(T[term]))*1./len(set(T[terms[i]])|set(T[term]))
                if val>max_jaccard[ii]:
                    max_jaccard[ii]=val
                if val>max_jaccard[i]:
                    max_jaccard[i]=val
            sim[j].append(sum(max_jaccard)*1./(i+1))
            
    fig, ax = plt.subplots()
    rects=[]
    for j in range(len(termss)):
        x=np.array(meanss[j])
        y=np.array(sim[j])
        if type(legend)==list and j%(len(termss)/len(legend))==1: #len(legend)==len(termss)
            rects.append(ax.plot(x,y,color=colors[j],alpha=opacity,label=legend[j/(len(termss)/len(legend))]))
        else:
            rects.append(ax.plot(x,y,color=colors[j],alpha=opacity))
    if invert_xaxis==True:
        ax.invert_xaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Average Maximal Jaccard Similarity Among Terms")
    #ax.set_ylabel(r"Ratio $\frac{\text{Explained Perturbed Genes}}{\text{Explained Unperturbed Genes}}$")
    if legend!="":
        ax.legend(legend,loc=0)
    return sim

def plot_optimal_beta(E,weights,N=100,normalize=True):    
    betas=np.linspace(0,1,N+1)[1:-1]
    y=[]
    k=len(weights)
    if normalize==True:
        sumE=sum(E)
        for i in range(len(E)):
            E[i]=E[i]*1./sumE
    for i in range(len(betas)):
        while betas[i]*weights[-1]>=1 and abs(betas[i]-betas[i-1])>1e-9:
            betas[i]=(betas[i]+betas[i-1])/2
        dummy=np.log(betas[i])*E[-1]
        for j in range(k):
            dummy+=np.log(1-betas[i]*weights[j])*E[j]
        y.append(dummy)
    return (betas,y)
    
def create_optimal_beta_plot_old(proportion_perturbed,weights,ratios,N=100,fontsize=12,colors=[],y_min=[]):
    import matplotlib.pyplot as plt

    if colors==[]:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]

    k=len(weights)
    res=[]
    for i in range(len(ratios)):
        E=[proportion_perturbed*ratios[i]*1./k for j in range(k)]
        E.append(1-proportion_perturbed)
        res.append(plot_optimal_beta(E,weights,N))
    for i in range(len(res)):
        hold
        plt.plot(res[i][0],res[i][1],color=colors[i],lw=2)
    for i in range(len(res)):
        max_index=res[i][1].index(max(res[i][1]))
        plt.plot(res[i][0][max_index],res[i][1][max_index],color=colors[i],marker='o', ls='')
    if y_min!=[]:    
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,y_min,y2))
        
def create_optimal_beta_plot(proportion_perturbed,weights,ratios,N=100,fontsize=12,colors=[],linestyles=[],y_min=[],y_max=10,x_max=0,legende=True):
    import matplotlib.pyplot as plt

    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)      
            
    #for weights use get_weights(np.ones(100,dtype=int),0,5)
    if type(proportion_perturbed) in [float,int]:
        proportion_perturbed=[proportion_perturbed]
    if colors==[]:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    if linestyles==[]:
        linestyles=['-','--',':']
    k=len(weights)
    res=[]
    for i in range(len(ratios)):
        res.append([])
        for j in range(len(proportion_perturbed)):
            E=[proportion_perturbed[j]*ratios[i]*1./k for iiii in range(k)]
            E.append(1-proportion_perturbed[j])
            res[i].append(plot_optimal_beta(E,weights,N))
    el_legend1,el_legend2=[],[]
    fig, ax = plt.subplots()
    for i in range(len(res)):
        for j in range(len(res[0])):
            hold
            if i==0:
                el_legend2.append(plt.plot(res[i][j][0],res[i][j][1],color='k',ls=linestyles[j],lw=2)[0])
            if j==0:
                el_legend1.append(plt.plot(res[i][j][0],res[i][j][1],color=colors[i],ls=linestyles[j],lw=2)[0])
            if j>0:
                plt.plot(res[i][j][0],res[i][j][1],color=colors[i],ls=linestyles[j],lw=2)
    for i in range(len(res)):
        for j in range(len(res[0])):
            max_index=res[i][j][1].index(max(res[i][j][1]))
            plt.plot(res[i][j][0][max_index],res[i][j][1][max_index],color=colors[i],marker='o', ls='')
    if y_min!=[]:    
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,max(x2,x_max),y_min,min(y_max,y2)))
    plt.xlabel(r'False Negative Rate $\beta$')
    plt.ylabel(r'Linear Transform of Likelihood Function')
    if legende==True:
        l1 = ax.legend(el_legend1,[r"Ratio$=$"+str(el) for el in ratios],loc=3)
        l2 = ax.legend(el_legend2,[r"$|P|/|V| =$"+str(el) for el in proportion_perturbed],loc=4)
        gca().add_artist(l1)

def get_weights(G_lengths,nr_categories,belief):
    s=nr_categories
    if nr_categories==0:
        s=sum(G_lengths)
    rates=[]
    if nr_categories>0: #FIXED error in formula, 32d_web and later
        nr_perturbed=sum(G_lengths[:nr_categories])
        summe=0
        for i in xrange(nr_categories):
            summe += (nr_categories-1-i+belief*i)*G_lengths[i]
        rates.append((nr_categories-1)*nr_perturbed*1./summe)
        #denom_rate_1=1 #OLD WRONG CODE, 32c_web and prior
        #for i in xrange(1,s):
        #    denom_rate_1 *= ((self.belief-1)*1.*i/(s-1) + 1)**(G_lengths[i]*1./nr_perturbed)
        #rates.append(falserate*1./denom_rate_1)
    else: #For singleton case, a simpler formula suffices since G_lengths[i]==1 for all i
        rates.append(2./(belief+1))
    for i in xrange(1,s): #Once rates[0] is found, multiply by linearly increasing constant to get rates[1],...,rates[s-1]
        rates.append(rates[0]*((belief-1)*1.*i/(s-1)+ 1))
    return rates

def find_optimal_beta(E,belief,N=100,accuracy=1e-4):    
    betas=np.array([0.001,0.25,0.5,0.5+0.999*0.5/belief])
    y=[]
    weights=get_weights(np.ones(N,dtype=int),0,belief)
    k=len(weights)
    while betas[-1]-betas[0]>accuracy:
        Ls=[sum([np.log(1-beta*weights[j])*E[j] for j in range(k)])+np.log(beta)*E[-1] for beta in betas]
        max_ind=Ls.index(max(Ls))
        if Ls[max_ind]<sum([np.log(1-(betas[max_ind]+1e-12)*weights[j])*E[j] for j in range(k)])+np.log(betas[max_ind]+1e-12)*E[-1]:
            betas=np.linspace(betas[max_ind],betas[max_ind+1],4)
        else:
            betas=np.linspace(betas[max_ind-1],betas[max_ind],4)
    return betas[max_ind]
        
            
def optimal_beta_against_quality(proportion_perturbed,beliefs,ratio_min=1,ratio_max=10,N=90,fontsize=12,colors=[],linestyles=[]):
    import matplotlib.pyplot as plt

    font = {'size'   : fontsize}

    matplotlib.rc('font', **font) 
    if colors==[]:
        colors=[rgb(212,175,255),rgb(176,232,249),rgb(255,224,191),rgb(255,176,176),rgb(169,235,196),rgb(153,76,0)]
    if linestyles==[]:
        linestyles=['-','--',':']
        
    ratios=np.linspace(ratio_min,ratio_max,N+1)
    res=[]
    for i in range(N+1):
        E=[proportion_perturbed*ratios[i]*1./N for iiii in range(N)]
        E.append(1-proportion_perturbed)
        res.append([])
        for belief in beliefs:        
            res[i].append(find_optimal_beta(E,belief,N))
    
    res_mat=np.matrix(res)
    legend_text=[]
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for i in range(len(beliefs)):
        legend_text.append(r'CRFE (belief =%i)' % beliefs[i] if beliefs[i]>1 else 'MGSA')
        y=np.asarray(res_mat[:,i]).reshape(-1)
        plt.plot(ratios,y,color=colors[i])
    plt.legend(legend_text)
    ax1.set_xlabel(r'|EP| / |EU|')
    ax1.set_ylabel(r'Optimal false negative rate $\beta^\star$')
    ax2 = ax1.twiny()
    ax2.axis(ax1.axis())

    new_tick_locations = np.array([str(round(val*(1-proportion_perturbed)/proportion_perturbed,2)) for val in ax1.get_xticks()])

    ax2.set_xticklabels(new_tick_locations)
    ax2.set_xlabel("Quality")
    
def optimal_beta_against_quality2(proportions_perturbed,beliefs,ratio_min=1,ratio_max=10,N=90,fontsize=12,colors=[],linestyles=[]):
    import matplotlib.pyplot as plt

    font = {'size'   : fontsize}

    matplotlib.rc('font', **font) 
    if colors==[]:
        colors=[rgb(212,175,255),rgb(176,232,249),rgb(255,224,191),rgb(255,176,176),rgb(169,235,196),rgb(153,76,0)]
    if linestyles==[]:
        linestyles=['-','--',':']
        
    ratios=np.linspace(ratio_min,ratio_max,N+1)
    res=[]
    for proportion_perturbed in proportions_perturbed:
        res.append([])
        for i in range(N+1):
            E=[proportion_perturbed*ratios[i]*1./N for iiii in range(N)]
            E.append(1-proportion_perturbed)
            res[-1].append([])
            for belief in beliefs:        
                res[-1][i].append(find_optimal_beta(E,belief,N))

    res_mat=[]
    for j in range(len(proportions_perturbed)):
        res_mat.append(np.matrix(res[j]))

    legend_text,symbols,rects=[],[],[]
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for j in range(len(proportions_perturbed)):
        for i in range(len(beliefs)):
            if j==0:
                legend_text.append(r'CRFE (belief =%i)' % beliefs[i] if beliefs[i]>1 else 'MGSA')
            y=np.asarray(res_mat[j][:,i]).reshape(-1)
            if i==0:
                symbols.append(ax1.plot(ratios,y,color='k',ls=linestyles[j])[0])
            if j==0:
                rects.append(ax1.plot(ratios,y,color=colors[i],ls=linestyles[j],lw=2)[0])
            else:
                ax1.plot(ratios,y,color=colors[i],ls=linestyles[j],lw=2)
    l1=ax1.legend(rects,legend_text)
    ax1.set_xlabel(r'|EP| / |EU|')
    ax1.set_ylabel(r'Optimal false negative rate $\beta^\star$')
    l2 = ax1.legend(symbols,[r"$|P|/|V|=$"+str(el) for el in proportions_perturbed],loc=5)#,numpoints = 1)
    gca().add_artist(l1)
    


    
    