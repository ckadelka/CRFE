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

def plotte_explained(C,T,MCMC_distr,levels_ordinal,n_intervals=10,cutoff=0):
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
    return (means,counter,width)
    
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

def plotte(data,N,what=[0],data1=[],STDS=True,data_G=[],data_G1=[],PROPORTIONS=False):
    import matplotlib.pyplot as plt
    if PROPORTIONS:
        for i in range(N):
            data[1][0][i]*=1./(data_G[0]-data_G[1])
            data[1][1][i]*=1./(data_G[0]-data_G[1])
            data[2][0][i]*=1./data_G[1]
            data[2][1][i]*=1./data_G[1]
            data1[1][0][i]*=1./(data_G1[0]-data_G1[1])
            data1[1][1][i]*=1./(data_G1[0]-data_G1[1])
            data1[2][0][i]*=1./data_G1[1]
            data1[2][1][i]*=1./data_G1[1]
    x=np.linspace(0,1,N)
    if data1==[]:
        colors=['b','r','g','m']
    else:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    for ind in what:
        plt.plot(x,data[ind][0],linestyle='-',color=colors[ind])
        if STDS:
            plt.plot(x,data[ind][0]+data[ind][1],linestyle='--',color=colors[ind])
            plt.plot(x,data[ind][0]-data[ind][1],linestyle='--',color=colors[ind])
        if data1!=[]:
            plt.plot(x,data1[ind][0],linestyle='-',color=colors[ind+4])
            if STDS:
                plt.plot(x,data1[ind][0]+data1[ind][1],linestyle='--',color=colors[ind+4])
                plt.plot(x,data1[ind][0]-data1[ind][1],linestyle='--',color=colors[ind+4])
    if PROPORTIONS:
        for i in range(N):
            data[1][0][i]*=(data_G[0]-data_G[1])
            data[1][1][i]*=(data_G[0]-data_G[1])
            data[2][0][i]*=data_G[1]
            data[2][1][i]*=data_G[1]
            data1[1][0][i]*=(data_G1[0]-data_G1[1])
            data1[1][1][i]*=(data_G1[0]-data_G1[1])
            data1[2][0][i]*=data_G1[1]
            data1[2][1][i]*=data_G1[1]

def plotte2(data,N,what=[0],data1=[],STDS=True,data_G=[],data_G1=[],PROPORTIONS=False):
    import matplotlib.pyplot as plt
    if PROPORTIONS:
        for i in range(N):
            data[1][0][i]*=1./(data_G[0]-data_G[1])
            data[1][1][i]*=1./(data_G[0]-data_G[1])
            data[2][0][i]*=1./data_G[1]
            data[2][1][i]*=1./data_G[1]
            data1[1][0][i]*=1./(data_G1[0]-data_G1[1])
            data1[1][1][i]*=1./(data_G1[0]-data_G1[1])
            data1[2][0][i]*=1./data_G1[1]
            data1[2][1][i]*=1./data_G1[1]
        label=['Jaccard Index of explained genes of different runs','Proportion Explained perturbed genes','Proportion Explained unperturbed genes','Number of explaining processes']
    else:
        label=['Jaccard Index of explained genes of different runs','Explained perturbed genes','Explained unperturbed genes','Number of explaining processes']
    x=np.linspace(0,1,N)
    fig, ax1 = plt.subplots()
    if data1==[]:
        colors=['b','r','g','m']
    else:
        colors=['b','r','g','m',[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    ind=what[0]
    ax1.plot(x,data[ind][0],linestyle='-',color=colors[ind])
    if STDS:
        ax1.plot(x,data[ind][0]+data[ind][1],linestyle='--',color=colors[ind])
        ax1.plot(x,data[ind][0]-data[ind][1],linestyle='--',color=colors[ind])
    if data1!=[]:
        ax1.plot(x,data1[ind][0],linestyle='-',color=colors[ind+4])
        if STDS:
            ax1.plot(x,data1[ind][0]+data1[ind][1],linestyle='--',color=colors[ind+4])
            ax1.plot(x,data1[ind][0]-data1[ind][1],linestyle='--',color=colors[ind+4])
    ax1.set_xlabel(r"$\tau$")
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
        ax2.plot(x,data[ind][0],linestyle='-',color=colors[ind])
        if STDS:
            ax2.plot(x,data[ind][0]+data[ind][1],linestyle='--',color=colors[ind])
            ax2.plot(x,data[ind][0]-data[ind][1],linestyle='--',color=colors[ind])
        if data1!=[]:
            ax2.plot(x,data1[ind][0],linestyle='-',color=colors[ind+4])
            if STDS:
                ax2.plot(x,data1[ind][0]+data1[ind][1],linestyle='--',color=colors[ind+4])
                ax2.plot(x,data1[ind][0]-data1[ind][1],linestyle='--',color=colors[ind+4])
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
            data[1][0][i]*=(data_G[0]-data_G[1])
            data[1][1][i]*=(data_G[0]-data_G[1])
            data[2][0][i]*=data_G[1]
            data[2][1][i]*=data_G[1]
            data1[1][0][i]*=(data_G1[0]-data_G1[1])
            data1[1][1][i]*=(data_G1[0]-data_G1[1])
            data1[2][0][i]*=data_G1[1]
            data1[2][1][i]*=data_G1[1]
    filename='aaabd.eps'
    plt.savefig(filename, bbox_inches=0)
    
def plotte3(data,N,what=[0],STDS=True,data_G=[],PROPORTIONS=False):
    import matplotlib.pyplot as plt
    if PROPORTIONS:
        for i in range(N):
            for j in range(len(data)):
                data[j][1][0][i]*=1./(data_G[j][0]-data_G[j][1])
                data[j][1][1][i]*=1./(data_G[j][0]-data_G[j][1])
                data[j][2][0][i]*=1./data_G[j][1]
                data[j][2][1][i]*=1./data_G[j][1]
        label=['Jaccard Index of explained genes of different runs','Proportion Explained perturbed genes','Proportion Explained unperturbed genes','Number of explaining processes']
    else:
        label=['Jaccard Index of explained genes of different runs','Explained perturbed genes','Explained unperturbed genes','Number of explaining processes']
    x=np.linspace(0,1,N)
    #fig, ax1 = plt.subplots()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    colors=['b','r','g','m',[0.4,0.4,1],[1,0.4,0.4],[0.4,1,0.4],[1,0.4,1],[0.7,0.7,1],[1,0.7,0.7],[0.7,1,0.7],[1,0.7,1]]
    ind=what[0]
    lines=[]
    for j in range(len(data)):
        lines.append(ax1.plot(x,data[j][ind][0],linestyle='-',color=colors[ind+4*j],label=str((j+1)*10)+'%'))
        if STDS:
            ax1.plot(x,data[j][ind][0]+data[j][ind][1],linestyle='--',color=colors[ind+4*j])
            ax1.plot(x,data[j][ind][0]-data[j][ind][1],linestyle='--',color=colors[ind+4*j])
    ax1.set_xlabel(r"$\tau$")
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
            lines.append(ax2.plot(x,data[j][ind][0],linestyle='-',color=colors[ind+4*j],label=str((j+1)*10)+'%'))
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
    ax1.legend(ls,labs,loc=0)
    #ax2.legend(['10%','20%','30%'],loc=2)
    plt.savefig(filename, bbox_inches=0)
    return lines
    
def plotte4(data_full,m,data_G_full=[],what=[0],which=1,s=0,perc_active=[10,20,30],N=100,STDS=True,PROPORTIONS=False,pos_leg=0):
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
    if what==[0,3] or what==[3,0]:
        identifier="jaccard"
    elif what==[1,2] or what==[2,1]:
        identifier="perturbed"
    if len(perc_active)>1:
        identifier+="_levels"
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
            elif w==2:
                labels.append(m.out_str.capitalize()+' - GO')
        title_end="different datasets"
    else:
        title.append(['HDF Data','3-cell vs 2-cell - GO','3-cell vs 2-cell - MSigDB','Alzheimers - GO'][which[0]])
        identifier+='_'+['hdf','liver_go','liver_msigdb',m.out_str+'_go'][which[0]]
    if len(which)==1 and len(s)>1:
        for w in s:
            data.append(data_full[which[0]][w][perc_active[0]])
            data_G.append(data_G_full[which[0]][w][perc_active[0]])            
            if w==0:
                labels.append('many categories')
            elif w==1:
                labels.append('1 category')
        title_end="different algorithms"
    else:
        title.append(['many categories','1 category'][s[0]])
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
    filename=identifier+"_cutoff"+str(m.cutoff)+"_topcutoff"+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+"_pert"+str(min_level)+'.eps'
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

nsim=20

cutoff=20
top_cutoff=200
param_learning=1
MCMC=1
belief=5
burnin=int(1*1e5)
steps=int(1*1e6)
weight_kind=0
category_kind=0
data1=[[[0 for iii in range(3)] for ii in range(2)] for i in range(4)]
data_full=[[[0 for iii in range(3)] for ii in range(2)] for i in range(4)]
data_G=[[[0 for iii in range(3)] for ii in range(2)] for i in range(4)]
for hdhsfs in [1,2,3]:
    for s in [0,1]:
        count=0
        for min_level in [10,20,30]:
            #if s==0 and weight_kind==0 or weight_kind==0 and s==1 and hdhsfs==0:
            #    continue 
            if hdhsfs==0:
                m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                            0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'bkz_gene_list.csv',0,
                                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                            'output1/', 'HDF')
            elif hdhsfs==1:
                m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                        0.25, 1e-50, belief, weight_kind, category_kind, 0.08016865, burnin, steps, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt',0,
                                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                            'output/', 'LIVER')
            elif hdhsfs==2:
                m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                        0.25, 1e-50, belief, weight_kind, category_kind, 0.08016865, burnin, steps, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12.txt',0,
                                            'msigdb.txt',
                                            'output/', 'LIVER')
            elif hdhsfs==3:
                m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                            0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'datasets/disease_obesity.txt',0,
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
                    f=open('saved_data/multiple_runs_msigdb_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
                else:
                    f=open('saved_data/multiple_runs_go_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
            except IOError:
                nsim=50
                if m.go_file=='msigdb.txt':
                    f=open('saved_data/multiple_runs_msigdb_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
                else:
                    f=open('saved_data/multiple_runs_go_'+m.out_str+'_s'+str(m.s)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+'_nsim'+str(nsim)+'.txt','rb')
            [C,MCMC_distr]=pickle.load(f)
            f.close()
            
            help=[max(MCMC_distr[i]) for i in range(len(MCMC_distr))]
            print hdhsfs,s,count,nsim,max(help),min(help)
            data_full[hdhsfs][s][count]=different_cutoffs(C,T,MCMC_distr,levels_ordinal,N=100)
            data1[hdhsfs][s][count]=stats(C,T,MCMC_distr,levels_ordinal,cutoff=0.5)#different_cutoffs(C,T,MCMC_distr,levels_ordinal,N=100)
            data_G[hdhsfs][s][count]=(lug,len(G[-1]))
            count+=1
            
def show(data,data_G,which,s,count):
    d=data[which][s][count]
    dG=data_G[which][s][count]

    print 'Jaccard Similarity:',d[0][0],'+-',d[0][1]
    print
    print 'Explained perturbed:',d[1][0],'+-',d[1][1]
    print 'Explained unperturbed:',d[2][0],'+-',d[2][1]
    print
    print 'Explained perturbed (in %):',d[1][0]*1./(dG[0]-dG[1]),'+-',d[1][1]*1./(dG[0]-dG[1])
    print 'Explained unperturbed (in %):',d[2][0]*1./dG[1],'+-',d[2][1]  *1./dG[1]  
    print
    print 'Processes:',d[3][0],'+-',d[3][1]
    
    


#def pr50(data):
#    for i in [49]:
#        print (data[1][0][i]+data[1][0][i+1])/2,(data[1][1][i]+data[1][1][i+1])/2
#        print (data[2][0][i]+data[2][0][i+1])/2,(data[2][1][i]+data[2][1][i+1])/2
#        print (data[3][0][i]+data[3][0][i+1])/2,(data[3][1][i]+data[3][1][i+1])/2
#        print (data[1][0][i]+data[1][0][i+1])/2/(lug-len(G[-1])),(data[1][1][i]+data[1][1][i+1])/2/(lug-len(G[-1]))
#        print (data[2][0][i]+data[2][0][i+1])/2/len(G[-1]),(data[2][1][i]+data[2][1][i+1])/2/len(G[-1])
#        print lug
#
#def pr10(data):
#    for i in [10]:
#        print (data[1][0][i]),(data[1][1][i])
#        print (data[2][0][i]),(data[2][1][i])
#        print (data[3][0][i]),(data[3][1][i])
#        print (data[1][0][i])/(lug-len(G[-1])),(data[1][1][i])/(lug-len(G[-1]))
#        print (data[2][0][i])/len(G[-1]),(data[2][1][i])/len(G[-1])

#for i in [10]:
#    print (data[1][0][i]),(data[1][1][i])
#    print (data[2][0][i]),(data[2][1][i])
#    print (data[3][0][i]),(data[3][1][i])
#    print (data[1][0][i])/(len(G)-1),(data[1][1][i])/(len(G)-1)
#    print (data[2][0][i])/(lug-len(G)+1),(data[2][1][i])/(lug-len(G)+1)