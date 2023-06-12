import GeneSetRankedEnrichment31m_web as tools
import random as r

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

MCMC=1

s=7
#m = tools.GeneSetRankedEnrichment(2, s, 5, 200, 0.05, 
#                                   0.25, 0.001, 2, 0, 0.15, 10, 25, 20000, 100000, 0.05157523,0.2,19,20, '3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt',
#                                   'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
#                                   'output/', 'LIVER')
#m = tools.GeneSetRankedEnrichment(2, s, 5, 200, 0.1, 
#                            0.25, 1e-5, 5, 0, 1.2, 1000, 1000, 100, int(1*1e5), 0.5,0.2,20,20,'datasets/disease_bipolar_disorder.txt',
#                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
#                            'output1/', 'DISEASE')
#m = tools.GeneSetRankedEnrichment(2, s, 5, 200, 0.05, 
#                                   0.25, 0.001, 2, 0, 0.15, 10, 25, 20000, 100000, 0.5,0.2,19,20, 'my_expression_data.txt',
#                                   'my_ontology.txt',
#                                   'output/', 'OWN')
m = tools.GeneSetRankedEnrichment(2, s, 5, 500, 0.1, 
                            0.25, 1e-10, 5, 0, 0.559, 1000, 1000, 20000, int(1e5), 0.5,0.2,20,20,'bkz_gene_list.csv',
                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                            'output1/', 'HDF')
#Usage:
#myGOAlgo = GeneSetRankedEnrichment(kind, s, bottom_cutoff, top_cutoff, alpha, beta, prob, belief, ?, min_level_for_activation, 
#kmax, Temp0, burnout, MCMC_steps, A,P, ?,?, gene_file, go_file, output_file, out_str)

#runMe(verbose, israndom, machine_learning, annealing_alg, MCMC, save_plots, show_plots)

filename='Heat_map_'+m.out_str+'_terms'

if MCMC==1:
    (C,MCMC_Distr)=m.runMe(0,0,1,0,MCMC,0,0)
else:
    (C,genes_in_C)=m.runMe(0,0,1,0,MCMC,0,0)
(genes_list, bounds, level_list) = m.get_sets_and_genes()

#bounds=[0,281]
#bounds.extend([481+200*i for i in xrange(s-1)])
#bounds.append(9330)

bounds=[]
bounds.extend([200*i for i in xrange(s+1)])
bounds.append(9330)

(T, unique_genes, large_goids, term_names) = m.load_all_new(genes_list)
#C=range(len(T))
lug=len(unique_genes)
Tset=set(range(len(T)))
cc=[len(T[i]) for i in xrange(len(T))]
if m.s==0:
    bounds[-1]=lug
(G,_, dict_G)=m.getG(unique_genes,genes_list,T,Tset,bounds)
#levels=[level_list[genes_list.index(unique_genes[i])] for i in xrange(len(unique_genes))]

glist=[genes_list.index(unique_genes[i]) for i in xrange(len(unique_genes))]
levels=[level_list[glist[i]] for i in xrange(len(unique_genes))]

#compare to heat plot of random ontology
for i in range(len(T)):
    T[i]=r.sample(range(lug),cc[i])
C=range(len(T))

mat=[]
for c in C:
    mat.append([])
    for g in G:
        mat[-1].append(len(set(g)&set(T[c])))
        
ccG=[len(g) for g in G]
ccC=[cc[c] for c in C]

mat_rat=[[max(mat[i][j]*1./ccC[i],0.000001) for j in xrange(len(mat[i]))] for i in xrange(len(mat))]
expected=[val*1./lug for val in ccG]

fold=[[round(1./expected[j]*mat_rat[i][j],2) for j in xrange(len(mat[i]))] for i in xrange(len(mat))]

f_out = open(filename+'.txt', 'w')
if MCMC:
    string='GO term \t MCMC density \t Total annotations per term & Total genes per category\t'
else:
    string='GO term \t Total annotations per term & Total genes per category\t'    
for val in ccG:
    string+=str(val)+'\t'
f_out.write(string[:-1]+'\n')
for i in range(len(fold)):
    if MCMC:
        string=term_names[C[i]]+'\t'+str(round(MCMC_Distr[i],4)) +'\t'+str(cc[C[i]])+'\t'
    else:
        string=term_names[C[i]]+'\t'+str(cc[C[i]])+'\t'
    for val in fold[i]:
        string+=str(val)+'\t'
    f_out.write(string[:-1]+'\n')

f_out.flush()
f_out.close()