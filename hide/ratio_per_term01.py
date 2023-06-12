import GeneSetRankedEnrichment32d_web as tools
import numpy as np
import graphics_gene_enrichment as gr
import matplotlib.pyplot as plt

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
min_level=0.3 #10: 10% cutoff, 20: 20% cutoff, 30: 30% cutoff
alpha_beta_top=0.5
belief=5
s=1

gene_files=['bkz_gene_list.csv','3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12.txt','datasets/disease_glaucoma.txt','datasets/disease_obesity.txt','3-cell-hep-d12_versus_CS-d12_ID.txt','datasets/disease_adenocarcinoma_of_lung.txt','datasets/disease_ischemic_cardiomyopathy.txt','3-cell-hep-d12_versus_2-cell-hep-d12_pos_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12_neg_ID.txt']
out_strs=['HDF','LIVER','LIVER','DISEASE','DISEASE','LIVER','DISEASE','DISEASE','Lpos','Lneg']

data_vec,qual_vec,ratios_vec=[],[],[]
nr_bins=40
    
#code, don't change
for which in [1,6]:
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
        
    if m.out_str=='DISEASE':
        m.out_str=m.gene_file.split('.')[0].split('datasets/disease_')[1]
        
    ratios=[]
    for i in range(len(m.T)):
        ratios.append(1-len(m.G[-1]&set(m.T[i]))*1./len(m.T[i]))
    qual_vec.append([el*1./(1-el)*(1./min_level-1) for el in ratios])
    ratios_vec.append(ratios)
    data_vec.append(tools.hist(ratios,nr_bins,0,1))
    
N = len(data_vec)
n_intervals=len(data_vec[0][0])

ind = np.arange(n_intervals)  # the x locations for the groups
width=0.85/N/nr_bins       # the width of the bars
colors=[gr.rgb(176,232,249),gr.rgb(255,176,176)]

rects=[]
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

for i in range(N):
    print data_vec[i]
    rects.append(ax1.bar((ind+i*width*nr_bins)/nr_bins, [el*1./sum(data_vec[i][0]) for el in data_vec[i][0]], width, color=colors[i]))
    #ax.bar([n_intervals+i*width+0.5], [data_vec[i][0][-1]*1./nr_unperturbed], width, yerr=data_vec[i][1][-1]*1./nr_unperturbed, error_kw={'ecolor':colors[i], 'capsize':4}, color=colors[i])    
ax1.set_xlabel(r'Proportion of Explained Genes that are Perturbed, $|EP|/|E|$')
ax1.set_ylabel('Proportion of Processes')
ax1.legend(rects,['Liver','Lung Cancer'])

#ax2.axis(ax1.axis())

new_tick_locations = np.array([val*min_level/(1+min_level*(val-1)) for val in [0.25,0.5,1,1.5,2,3,5,7]])

def tick_function(X):
    res = (X/min_level)/((1-X)/(1-min_level))
    #res = [round(el,2) for el in res]
    return [str(z) for z in res]

ax1.set_xticks(new_tick_locations)
ax1.set_xticklabels(tick_function(new_tick_locations))
ax1.set_xlabel(r"Quality of Explanatory Set, $(|EP||U|)/(|P||EU|)$")

x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,0.25))

ann=[0 for i in m.unique_genes]
for t in m.T:
    for gene in t:
        ann[gene]+=1
        
ind=sorted(range(len(ann)), reverse=False, key=lambda k: m.levels_ordinal[k])
ann_sorted=[ann[i] for i in ind]

mean_1=np.mean(ann_sorted[:int(len(m.unique_genes)*m.min_level*1./10*3)])
mean_2=np.mean(ann_sorted[int(len(m.unique_genes)*m.min_level*1./10*3):int(len(m.unique_genes)*m.min_level)])


