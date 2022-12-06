import GeneSetRankedEnrichment30 as mytools
import optparse
import matplotlib.pyplot as plt

def runPrecRec(myGOAlgo,N=5,output=1,percentage_high_active=0.2,nr_active=10,MACHINE_LEARNING=False,MCMC=0):
    myGOAlgo.out_str='FULL'
    #(genes_list, _,_) = myGOAlgo.get_sets_and_genes()
    (T, unique_genes, large_goids, term_names) = myGOAlgo.load_all([])
    Tset=set(range(len(T)))
    (recall,precision,active_GOs,found_GOs,p_values,pos_prec,pos_prec_p)=myGOAlgo.prec_rec_plot3(unique_genes,T,Tset,term_names,N,output,percentage_high_active,nr_active,False,MACHINE_LEARNING,DATA_FROM_FILE=True,MCMC=MCMC)
    return (recall,precision,pos_prec,pos_prec_p)

def main():
    '''
        This kicks off the process and parses the options from the arguments.
    '''
    p = optparse.OptionParser()
    
    p.add_option('--cutoff', '-c', default ='200', help='only GO terms above this cutoff are considered (default 200)')
    p.add_option('--top_cutoff', '-z', default ='0', help='only GO terms below this cutoff are considered (default 0 - meaning all)')

    p.add_option('--kind', '-k', default = '2', help='Algorithm kind. Set to 1 or 2. 1: GenGO 2: GOing Bayesian (default)')
    #p.add_option('--activityLevels', '-l', default='2')
    
    p.add_option('--alpha', '-a', default = '0.01', help='alpha value for GOing Bayesian algo (default 0.01)')
    p.add_option('--beta', '-b', default = '0.1', help='beta value for GOing Bayesian algo (default 0.1)')
    p.add_option('--prob', '-p', default = '0.01', help='prob p value for GOing Bayesian algo (default 0.01)')
    
    p.add_option('--gene_file', '-e', default='db-brass-konig-zhou-mouse-ovn-ss-w-1-unweighted-net-predictions-annotated.csv', help='file string for a ranked list of genes')
   
    p.add_option('--go_file', '-g', default='human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv', help='csv file for the go terms')
    p.add_option('--out', '-o', default='results_GOEnrichment.txt', help='file string to save the output. Defaults to results_GOEnrichment.txt')
    p.add_option('--save', '-d' ,default = 'FULL', help='string to concatenate to the files created by the program (empty default)')
    p.add_option('--belief', '-l', default='2', help='Belief for how much more active the highest gene set is compared to the least active gene set (default 2)')
    p.add_option('-r', action='store_true', dest='israndom', help='When turned on, the algorithm starts with a random set of GO terms')
    p.add_option('-v', action='store_true', dest='verbose', help='Turns on verbose output. Will slow down the algorithm')
    options, arguments = p.parse_args()
    #cutoff = options.cutoff
    israndom = options.israndom
    verbose = options.verbose
    try:
        cutoff = int(options.cutoff)
    except ValueError:
        cutoff = 200
    try:
        top_cutoff = int(options.top_cutoff)
    except ValueError:
        top_cutoff = 0    
    try:
        wt = int(options.weight_function)
    except ValueError:
        wt = 0
    
    try:
        kind = int(options.kind)
    except ValueError:
        kind = 2
    if (verbose):
        print 'Running kind: ' + str(kind)
# try:
# s = int(options.activityLevels)
# except ValueError:
# s = 2

    try:
        alpha = int(options.alpha)
    except ValueError:
        alpha = 0.01
        
    try:
        beta = int(options.beta)
    except ValueError:
        beta = 0.1
        
    try:
        prob = int(options.prob)
    except ValueError:
        prob = 0.01
        
    try:
        belief = int(options.belief)
    except ValueError:
        belief = 2
    gene_file=options.gene_file
    go_file = options.go_file
    out_str = options.save
    output_file = options.out

    m = mytools.GeneSetRankedEnrichment(kind, 400, top_cutoff,wt, p, q, a, alpha, beta, prob, belief, gene_file, go_file, output_file, out_str)
    
    runPrecRec(m,100,0)

def display_plots(kind,N,a_prob,cutoff,top_cutoff=0,nr_active=1):
    alphas=[0.1,0.4]
    betas=[0.25,0.4]
    nr_actives=[10,20]
    f, A = plt.subplots(len(alphas), len(nr_actives), sharex='col', sharey='row')
    title=r'$N=$%s, $p=$%s, %s$\geq$size of GO terms$\geq$%s' % (N,a_prob*1./1000000,top_cutoff if top_cutoff>cutoff else r'$\infty$',cutoff)
    f.text(0.5,0.975,title,horizontalalignment='center',verticalalignment='top')
    plt.axis([0,1,0,1])
    for i in xrange(len(alphas)):
        #plt.ylabel("Precision"+str(i))
        for j in xrange(len(nr_actives)):
            [rec,prec]=mytools.GeneSetRankedEnrichment.load_rec_prec(kind,N,alphas[i],betas[i],a_prob,cutoff,top_cutoff,nr_actives[j])
            A[i][j].plot(rec,prec,'rx')
            A[i][j].set_title(r'$\alpha=$%s, $\beta=$%s, nr$=$%s' % (alphas[i],betas[i],nr_actives[j]), fontsize=12)
    plt.show()
    return

def gen_generative_data(N,nr_active,cutoff,top_cutoff,alpha,beta):
    import cPickle as pickle
    m = mytools.GeneSetRankedEnrichment(1, s, cutoff, top_cutoff, alpha, beta, prob, 2, 0.6, 1000, 1000, 20000, 1000000, 19,20, 'db-brass-konig-zhou-mouse-ovn-ss-w-1-unweighted-net-predictions-annotated.csv','human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv','output/', 'FULL')
    (T, unique_genes, large_goids, term_names) = m.load_all([])
    Tset=set(range(len(T)))
    Gs=[]
    active_GOs=[]
    for i in xrange(N):
        (G,Sn_max,active_GO)=m.getgenG2(unique_genes,T,Tset,percentage_high_active=0.2,nr_active=nr_active)
        Gs.append(G)
        active_GOs.append(active_GO)
    f1=open('saved_data/save_generative_data_G.txt','w')
    pickle.dump(Gs, f1)
    f1.close()
    f2=open('saved_data/save_generative_data_nr_active.txt','w')
    pickle.dump(active_GOs, f2)
    f2.close()
    print "Data generated"
    return
    
#if __name__ == '__main__':
#    main()

#FULL
N=50
cutoff=20
top_cutoff=200
s=10
prob=0.01

results=[]

for nr_active in [10]: #3
    for top_cutoff in [200]: #3
        for alpha,beta in zip([0.1],[0.25]):
            #generate generative data and use the same for all following precision-recall plots
            gen_generative_data(N,nr_active,cutoff,top_cutoff,alpha,beta)
            for kind in [2]:
                for MACHINE_LEARNING in [0]:
                    for MCMC in [0]:
                        for s in [nr_active,1]:
                            m = mytools.GeneSetRankedEnrichment(kind, s, cutoff, top_cutoff, alpha, beta, prob, 2, 0.6, 1000, 1000, 20000, 1000000, 19,20,'db-brass-konig-zhou-mouse-ovn-ss-w-1-unweighted-net-predictions-annotated.csv','human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv','output/', 'FULL')
                            results.append(runPrecRec(m,N,0,0.2,nr_active,MACHINE_LEARNING,MCMC))

##m.runMe(1,0,0,0,0,0,0)
unique_genes=m.loadvar('unique_genes')
T=m.loadvar('T')
Tset=set(range(len(T)))
large_goids=m.loadvar('large_goids')
term_names=m.loadvar('term_names')
cc=[len(T[i]) for i in xrange(len(T))]
#m.prec_rec_plot(unique_genes,T,Tset,term_names,N,m.p,m.q,1)

import cPickle as pickle
f1=open('saved_data/save_generative_data_G.txt','r')
Gs=pickle.load(f1)
f1.close()
f2=open('saved_data/save_generative_data_nr_active.txt','r')
active_GOs=pickle.load(f2)
f2.close()

f=open('saved_data/results_N'+str(N)+'_cutoff'+str(cutoff)+'_topcutoff'+str(top_cutoff)+'.txt','w')
pickle.dump(results, f)
f.close()

sizes=[[cc[c] for c in active_GOs[i]] for i in xrange(N)]

f=open('saved_data/results_N'+str(N)+'_cutoff'+str(cutoff)+'.txt','r')
res=pickle.load(f)
f.close()