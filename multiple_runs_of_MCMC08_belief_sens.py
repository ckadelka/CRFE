import GeneSetRankedEnrichment31s_web as tools
import numpy as np
import cPickle as pickle
import graphics_gene_enrichment as gr

nsim=50

cutoff=5
top_cutoff=200
param_learning=1
MCMC=1
weight_kind=1
category_kind=0
s=0
which=0
belief=5
burnin=int(1*1e5)
steps=int(1*1e6)
min_level=20

Cs,MCMC_distrs,mean_expl,std_expl,meanss,termss,stdss=[[0 for i in range(4)] for j in range(4)],[[0 for i in range(4)] for j in range(4)],[[0 for i in range(4)] for j in range(4)],[[0 for i in range(4)] for j in range(4)],[[0 for i in range(4)] for j in range(4)],[[0 for i in range(4)] for j in range(4)],[[0 for i in range(4)] for j in range(4)]

for which in [1,2,3]:
    for jj in [0,1,2]:
        belief=[2,5,10][jj]
        if which==0:
            m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                        0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'bkz_gene_list.csv',0,
                                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                        'output1/', 'HDF')
        elif which==1:
            m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                    0.25, 1e-50, belief, weight_kind, category_kind, 0.08016865, burnin, steps, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt',0,
                                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                        'output/', 'LIVER')
        elif which==2:
            m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                                    0.25, 1e-50, belief, weight_kind, category_kind, 0.08016865, burnin, steps, 0.5,0.2,20,20,'3-cell-hep-d12_versus_2-cell-hep-d12.txt',0,
                                        'msigdb.txt',
                                        'output/', 'LIVER')
        elif which==3:
            m = tools.GeneSetRankedEnrichment(2, s, cutoff, top_cutoff, 0.1, 
                        0.25, 1e-50, belief, weight_kind, category_kind, 0.559, burnin, steps, 0.5,0.2,20,20,'datasets/disease_obesity.txt',0,
                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                        'output1/', 'DISEASE')
        
        
        (genes_list, level_list) = m.get_sets_and_genes()
        (T, unique_genes, large_goids, term_names) = m.load_all_new(genes_list)
        lug=len(unique_genes)
        m.min_level_for_activation=gr.correct_cutoff(min_level/100.,unique_genes,level_list,genes_list,m.kind_of_list)
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
        [Cs[which][jj],MCMC_distrs[which][jj]]=pickle.load(f)
        C=Cs[which][jj]
        MCMC_distr=MCMC_distrs[which][jj]
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
        termss[which][jj]=terms[:]
        meanss[which][jj]=means[:]
        stdss[which][jj]=stds[:]
        mean_expl[which][jj]=(np.mean(mat_expl,0))
        std_expl[which][jj]=(np.std(mat_expl,0))
        
        print which,jj

#Sensitivity Analysis w.r.t belief or whatever second variable (jj) was
def sensitivity_plot(which,fixated,termss,meanss,m,stdss,N=10,loc=3,mini=0,maxi=1.2,title="",identifier="",SAVE=False):
    import matplotlib.pyplot as plt

    y=[]
    y_std=[]
    for i in range(len(termss[which])):
        y.append([])
        y_std.append([])
        for j in range(N):
            if i!=fixated:
                ind=termss[which][i].index(termss[which][fixated][j])
            else:
                ind=j
            y[i].append(meanss[which][i][j])
            y_std[i].append(stdss[which][i][j])
            
    colors=[[0,0,0.7],[0.4,0.4,0.9],[0.7,0.7,1],[0.9,0.9,1]]
    
    ind = np.arange(N)  # the x locations for the groups
    width=0.85/len(termss[which])       # the width of the bars
    
    rects=[]
    fig, ax = plt.subplots()
    for i in range(len(y)):
        rects.append(ax.bar(ind+i*width, y[i], width, yerr=y_std[i], error_kw={'ecolor':colors[i], 'capsize':4}, color=colors[i]))
    ax.set_xticks(ind+0.425)
    ax.set_xticklabels( [str(i+1) for i in range(N)] )
    leg_text=['belief '+str(b) for b in [2,5,10]]
    ax.legend(leg_text,loc=loc)
    ax.set_ylabel('Posterior Probability')
    title=r'Posterior Probability of top terms for different belief parameters'+('\n '+title if title!="" else "")
    ax.set_title(title)
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,mini,maxi))
    if SAVE:
        filename="sens_belief_"+identifier
        filename+="_cutoff"+str(m.cutoff)+"_topcutoff"+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+"_pert"+str(min_level)+"_N"+str(N)
        filename+='.eps'
        plt.savefig(filename, bbox_inches=0)    
    return y

def sensitivity_plot2(which,termss,meanss,stdss,m,N=10,loc=3,mini=0,maxi=1.2,title="",identifier="",SAVE=False):
    import matplotlib.pyplot as plt

    y=[]
    y_std=[]
    terms=[]
    for i in range(len(termss[which])):
        y.append([0 for j in range(len(termss[which])*N)])
        y_std.append([0 for j in range(len(termss[which])*N)])
        for j in range(N):
            try:
                ind=terms.index(termss[which][i][j])
            except:
                ind=len(terms)
                terms.append(termss[which][i][j])
            y[i][ind]=meanss[which][i][j]
            y_std[i][ind]=stdss[which][i][j]
    
    N=len(terms)
    for i in range(len(y)):
        y[i]=y[i][:N]
        y_std[i]=y_std[i][:N]
    
    colors=[[0,0,0.7],[0.4,0.4,0.9],[0.7,0.7,1],[0.9,0.9,1]]
    
    ind = np.arange(N)  # the x locations for the groups
    width=0.85/len(termss[which])       # the width of the bars
    
    rects=[]
    fig, ax = plt.subplots()
    for i in range(len(y)):
        rects.append(ax.bar(ind+i*width, y[i], width, yerr=y_std[i], error_kw={'ecolor':colors[i], 'capsize':4}, color=colors[i]))
    ax.set_xticks(ind+0.425)
    ax.set_xticklabels( [str(i+1) for i in range(N)] )
    leg_text=['belief '+str(b) for b in [2,5,10]]
    ax.legend(leg_text,loc=loc)
    ax.set_ylabel('Posterior Probability')
    title=r'Posterior Probability of top terms for different belief parameters'+('\n '+title if title!="" else "")
    ax.set_title(title)
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,mini,maxi))
    if SAVE:
        filename="sens_belief_"+identifier
        filename+="_cutoff"+str(m.cutoff)+"_topcutoff"+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+"_pert"+str(min_level)+"_N"+str(N)
        filename+='.eps'
        plt.savefig(filename, bbox_inches=0)    
    return y