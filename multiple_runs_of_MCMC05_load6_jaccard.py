import GeneSetRankedEnrichment32d_web as tools
import numpy as np
import cPickle as pickle
import graphics_gene_enrichment as gr
from matplotlib.pyplot import *

OUTPUT_FOR_PUBLICATION=False
if OUTPUT_FOR_PUBLICATION:
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
else:
    matplotlib.rcParams['pdf.fonttype'] = 3
    matplotlib.rcParams['ps.fonttype'] = 3

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

for which in range(6,7):
    print 'WHICH',which
    #leave these parameters fixed
    cutoff=20
    top_cutoff=200
    burnin=int(1*1e5)
    steps=int(1*1e6)
    
    #you can change these parameters
    weight_kind=1 #1 is what we used
    category_kind=0
    
    #really change these parameters
    #which=17 #1: Liver, 6: Lung, >9: NCBI new datasets
    min_level=0.3 #10: 10% cutoff, 20: 20% cutoff, 30: 30% cutoff
    alpha_beta_top=0.5
    
    nsim=20
    
    gene_files=['bkz_gene_list.csv','3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12.txt','datasets/disease_glaucoma.txt','datasets/disease_obesity.txt','3-cell-hep-d12_versus_CS-d12_ID.txt','datasets/disease_adenocarcinoma_of_lung.txt','datasets/disease_ischemic_cardiomyopathy.txt','3-cell-hep-d12_versus_2-cell-hep-d12_pos_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12_neg_ID.txt']
    gene_files.extend(['gene_file_GDS3004.txt','gene_file_GDS4419.txt','gene_file_GDS4610.txt','gene_file_GDS4974.txt'])
    gene_files.extend(['gene_file_GDS1686.txt','gene_file_GDS2969.txt','gene_file_GDS3216.txt','gene_file_GDS3246.txt','gene_file_GDS3866.txt','gene_file_GDS3928.txt','gene_file_GDS3933.txt','gene_file_GDS4423.txt','gene_file_GDS4929.txt','gene_file_GDS5092.txt'])
    out_strs=['HDF','LIVER','LIVER','DISEASE','DISEASE','LIVER','DISEASE','DISEASE','Lpos','Lneg','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI']
    
    #for the category plot:
    n_intervals=10 #number of different groups of active categories
    threshold=0.2
    
    termss,meanss,stdss,data,cats,legend,cats,Cs,MCMC_distrs=[],[],[],[],[],[],[],[],[]
    
    for belief,s in [(10,0),(5,0),(2,0),(100,1),(9,0)]:
        min_level=0.3
        #code, don't change
        m = tools.GeneSetRankedEnrichment(s, cutoff, top_cutoff, 0.1, 
                                        0.25, 1e-50, belief, weight_kind, category_kind, min_level, 'percentage',burnin, steps, alpha_beta_top,0.2,20,20,gene_files[which],0,
                                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                        'output1/', out_strs[which])
        if which==2:
            m.go_file='msigdb.txt'                   
        elif m.out_str=='NCBI':
            if int(m.gene_file[13:17]) in [3004,4419,4610,4974]:
                m.go_file='ONTOLOGY_PARSER/annotation_file_human_association_human_biological_process.txt'
            elif int(m.gene_file[13:17]) in [2969,3866]:
                m.go_file='ONTOLOGY_PARSER/annotation_file_yeast_association_yeast_biological_process.txt'
            elif int(m.gene_file[13:17]) in [1686]:
                m.go_file='ONTOLOGY_PARSER/annotation_file_fly_association_fly_biological_process.txt'
            elif int(m.gene_file[13:17]) in [3928]:
                m.go_file='ONTOLOGY_PARSER/annotation_file_rat_association_rat_biological_process.txt'
            elif int(m.gene_file[13:17]) in [3246,4423,4929,5092]:
                m.go_file='ONTOLOGY_PARSER/annotation_file_mouse_association_mouse_biological_process.txt'
            elif int(m.gene_file[13:17]) in [3216,3933]:
                m.go_file='ONTOLOGY_PARSER/annotation_file_arabidopsis_association_arabidopsis_biological_process.txt'
            m.gene_file='datasets/'+m.gene_file
            m.out_str+='_'+m.gene_file.split('_')[-1][:-4]
                                
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
        Cs.append(C)
        MCMC_distrs.append(MCMC_distr)
        
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
        #cats.append(gr.categories(C,m.T,MCMC_distr,m.levels_ordinal,len(m.unique_genes)-len(m.G[-1]),len(m.unique_genes),n_intervals,cutoff=threshold,kind=3))
        if which==1:
            identifier="liver"
            title="Liver models"
        elif which==6:
            identifier="lung_cancer"
            title="Adenocarcinoma of Lung"
        elif which>9:
            identifier="NCBI_"+m.gene_file.split('_')[-1][:-4]
            title='NCBI - '+m.gene_file.split('_')[-1][:-4]
        if s==1:
            identifier2="MGSA"
        elif s==0.5:
            identifier2="FA 2.0"
        elif s==1.5:
            identifier2="GOrilla"        
        else:
            identifier2="CRFE"
        legend.append('CRFE (belief='+str(belief)+')' if s==0 else identifier2)
        #gr.plot_jaccard_term_by_term2a(C,m.T,MCMC_distr,m.G,m,min_level=20,cutoff=10,n_intervals=10,kind=1,colors=[gr.rgb('e')],title=legend[-1],identifier=identifier,SAVE=True,fontsize=54)
        #gr.plot_jaccard_term_by_term2a(C,m.T,MCMC_distr,m.G,m,min_level=20,cutoff=20,n_intervals=10,kind=1,colors=[gr.rgb('e')],title=legend[-1],identifier=identifier,SAVE=True,fontsize=54)
        #gr.plot_jaccard_term_by_term2a(C,m.T,MCMC_distr,m.G,m,min_level=20,cutoff=40,n_intervals=10,kind=1,colors=[gr.rgb('e')],title=legend[-1],identifier=identifier,SAVE=True,fontsize=54)
        #gr.plot_jaccard_term_by_term2a(C,m.T,MCMC_distr,m.G,m,min_level=30,cutoff=10,n_intervals=10,kind=1,colors=[gr.rgb('e')],title=legend[-1],identifier=identifier,SAVE=True,STD=False,fontsize=54)
        #gr.plot_jaccard_term_by_term2a(C,m.T,MCMC_distr,m.G,m,min_level=30,cutoff=20,n_intervals=10,kind=1,colors=[gr.rgb('e')],title=legend[-1],identifier=identifier,SAVE=True,STD=False,fontsize=54)
        #gr.plot_jaccard_term_by_term2a(C,m.T,MCMC_distr,m.G,m,min_level=30,cutoff=40,n_intervals=10,kind=1,colors=[gr.rgb('e')],title=legend[-1],identifier=identifier,SAVE=True,STD=False,fontsize=54)
        #
        #a=gr.plot_jaccard_term_by_term2(C,m.T,MCMC_distr,m.G,m,min_level,cutoff=20,kind=1,title=legend[-1],identifier=identifier+"_"+identifier2+"suppl",SAVE=True,fontsize=17,colors=[gr.rgb('p'),gr.rgb(235,235,235)],STD=True)    
    
    #colorlist=[gr.rgb(212,175,255),gr.rgb(176,232,249),gr.rgb(255,224,191),gr.rgb(255,176,176)]
    ##
    #for top_proportion,bottom_proportion in [[0,1]]:#zip([0,0,0,0,0.5,0],[1,0.1,0.25,0.5,1,0.25]):
    #    a=gr.pert_vs_unpert2b_ratio(termss,m.G,m.T,top_proportion,bottom_proportion,legend,meanss,fontsize=14,colors=colorlist,coarsity=1,show_annotations=0,identifier=identifier,title=title,m=m,SAVE=True,max_x=1.2,min_level=30,pert_or_unpert_on_xaxis='pert')
    
        #(peps,peus,xs,ys,ns)=gr.combine_pep_lines(Cs,m.G,m.T,top_proportion,bottom_proportion,MCMC_distrs=MCMC_distrs)
        #gr.pert_vs_unpert2b_ratio_separate2(xs,ys,peus,m.G,m.T,top_proportion,bottom_proportion,legend=legend,fontsize=14,MCMC_distrs=MCMC_distrs,coarsity=1,colors=colorlist,opacity=1,highlights=[0.75,0.5,0.1,0.01],identifier=identifier,title=title,m=m,SAVE=True,max_x=0.5,min_level=30)
    
    #pert25=set.union(*[set([])] + [m.G[i] for i in range(int((len(m.unique_genes)-len(m.G[-1]))*0.25))])
    #pert=set.union(*[set([])] + [m.G[i] for i in range(int((len(m.unique_genes)-len(m.G[-1]))))])
    #unpert=set(m.G[-1])
    
    #ind=0
    #EP,EU,EP25=set(),set(),set()
    #for i,t in enumerate(termss[ind][:50]):
    #    new_EP25=set.difference(set(m.T[t])&pert25,EP25)
    #    new_EP=set.difference(set(m.T[t])&pert,EP)
    #    new_EU=set.difference(set(m.T[t])&unpert,EU)
    #    EP25=set.union(EP25,new_EP25)
    #    EP=set.union(EP,new_EP)
    #    EU=set.union(EU,new_EU)
    #    dummy=[m.levels_ordinal[el] for el in m.T[t]]
    #    dummy.sort()
    #    print len(set(m.T[t])&pert25),len(set(m.T[t])&pert),len(set(m.T[t])&unpert),'|||',len(new_EP25),len(new_EP),len(new_EU),'|||',len(EP25)*1./len(EU)*7./3.*4,len(EP)*1./len(EU)*7./3,'|||',t,meanss[ind][i],'|||',mean(dummy[:len(set(m.T[t])&pert25)])/len(pert25),(mean(dummy[len(set(m.T[t])&pert25):len(set(m.T[t])&pert)])-len(pert25))/(len(pert)-len(pert25)),mean(dummy[:len(set(m.T[t])&pert)])/len(pert),'|||',len(set(m.T[termss[ind][i]])&pert&set.union(*[set([])] + [set(m.T[termss[ind][j]]) for j in range(100) if i!=j]))*1./len(set(m.T[t])&pert)
    
    
    #To produce a quality plot with CRFE, MGSA AND GOrilla and FA2.0, first run
    #%run Revigo_xls_output_transformer00.py
    #C_REVIGO=C[:]
    #%run FA_xls_output_transformer01.py
    #C_FA=C[:]
    #%run Gorilla_xls_output_transformer02.py
    #C_GOrilla=C[:]
    #%run GSEA_txt_output_transformer00.py
    #C_GSEA=C[:]
    #%run this file without any plotting
    
    colorlist=[gr.rgb(212,175,255),gr.rgb(176,232,249),gr.rgb(255,224,191),gr.rgb(255,176,176),gr.rgb(255,105,180),gr.rgb(153,76,0),gr.rgb(55,240,40),gr.rgb(100,100,255)]
    termss.extend([C_FA,C_GOrilla,C_GSEA,C_REVIGO])
    meanss.extend([[0 for i in range(len(C_FA))],[0 for i in range(len(C_GOrilla))],[0 for i in range(len(C_GSEA))],[0 for i in range(len(C_REVIGO))]])
    legend.extend(['FuncAssociate2.0','GOrilla','GSEA','GSEA+REVIGO'])
    for top_proportion,bottom_proportion in zip([0,0,0],[1,0.25,0.5]):
        a=gr.pert_vs_unpert2b_ratio(termss,m.G,m.T,top_proportion,bottom_proportion,"",meanss,fontsize=14,colors=colorlist,coarsity=1,show_annotations=0,identifier=identifier+'_all',title=title,m=m,SAVE=True,max_x=1.2,min_level=30,pert_or_unpert_on_xaxis='pert')

identifier2=['crfe10','crfe5','crfe2','mgsa','fa','gorilla','gsea','revigo']
for ind in range(len(termss)):
    a=gr.pert_vs_unpert2bb_ratio([termss[ind]],m.G,m.T,legend=[legend[ind]],top_proportions=[0,0,0],bottom_proportions=[0.25,0.5,1],means=[meanss[ind]],fontsize=14,colors=[colorlist[ind]],coarsity=1,show_annotations=0,identifier=identifier+'_'+identifier2[ind],title=title,m=m,SAVE=True,max_x=1,min_level=30,pert_or_unpert_on_xaxis='pert')



##Count annotations per gene
#counter=[0 for i in range(len(m.unique_genes))]
#for t in m.T:
#    for gene in t:
#        counter[gene]+=1
#
#sizeP=len(m.G)-1
#sizeP25=sizeP/4
#sizeP50=sizeP/2
#
#a=[]
#a25=[]
#a50=[]
#for i in range(sizeP):
#    if i<sizeP25:
#        a25.append(counter[list(m.G[i])[0]])
#    elif i<sizeP50:
#        a50.append(counter[list(m.G[i])[0]])
#    else:
#        a.append(counter[list(m.G[i])[0]])
#        
#print np.mean(a25),np.mean(a25)*0.5+np.mean(a50)*0.5,np.mean(a25)*0.25+np.mean(a50)*0.25+np.mean(a)*0.5
#    