import GeneSetRankedEnrichment23 as mytools

for kind in [1,2]:
    for topcutoff in [0,500,200]:
        for s in [0,1,2,4]:
            for belief in [2,5]:
                for prob in [0.01,0.001,0.0001]:
                    print kind,topcutoff,s,belief,prob
                    if prob!=0.01 and kind==1: #different probs only for Bayesian
                        continue
                    m = mytools.GeneSetRankedEnrichment(kind, s, 5, topcutoff, 0, 0.5, 0.15, 3, 0.15,
                                                        0.5, prob, belief, 0.6, 1000, 1000, 2,
                                                        'bkz_gene_list.csv',
                                                        'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                                        'set_bounds_file.txt',
                                                        'results_GOEnrichment.txt',
                                                        'HDF')
                    C=m.runMe(0,0,0,0,1,0)
                    if kind==1 and belief==2: #only do it once
                        unique_genes=m.loadvar('unique_genes')
                        lug=len(unique_genes)
                        T=m.loadvar('T')
                        large_goids=m.loadvar('large_goids')
                        Tset=set(range(len(T)))
                        term_names=m.loadvar('term_names')
                        cc=[len(T[i]) for i in xrange(len(T))]
                        (genes_list, bounds, level_list) = m.get_sets_and_genes()
                        if m.s==0:
                            bounds[-1]=lug
                        (G,_, dict_G)=m.getG(unique_genes,genes_list,T,Tset,bounds)
                        ##m.alg(T,G,len(unique_genes),sum([len(set(G[m.s])&set(T[cn])) for cn in Tset]),term_names,[],1)
                        if m.s==2:
                            category_mat=[[0],[1],[2],[0,1],[0,2],[1,2],[0,1,2]]
                        elif m.s==1:
                            category_mat=[[0],[1],[0,1]]
                        elif m.s>0:
                            category_mat=[[0],[m.s],range(m.s),range(m.s+1)]
                        else:
                            category_mat=[range(len(dict_G)),[len(dict_G)],range(len(dict_G)+1)]
                        C_Gorilla=mytools.format_gorilla_output('GO.csv',large_goids)
                        title=r'Gorilla results (p-value$<10^{-11}$)'
                        m.plot_jaccard_comp(C_Gorilla,category_mat,G,T,20,kind=1,title=title,dict_G=dict_G,save=True,algorithm='Gorilla')
