import csv
import cPickle as pickle

def uniq(input):
        temp = set(input)
        return list(temp)

def savevar( variable, v):
        f=open('saved_data/save_GO_BP_'+variable+'.txt','w+')
        pickle.dump(v, f)
        f.close()

def find_equal_terms(T,term_names):
    cc=[len(T[i]) for i in xrange(len(T))]
    ind=sorted(range(len(cc)), key=lambda k: cc[k])
    cc_ord=[cc[ind[i]] for i in xrange(len(cc))]
    res=[]
    for i in xrange(len(cc)):
        for j in xrange(i+1,len(cc)):
            if cc_ord[i]<cc_ord[j]:
                break
            elif jaccard([ind[i],ind[j]],T,term_names,0)==1:
                res.append((ind[i],ind[j],cc_ord[i]))
    return res

def remove_equal_terms(T,term_names,large_goids,output=False):
    #print "merge equal GO terms into one..."
    res=find_equal_terms(T,term_names)
    l0=[item[0] for item in res]
    l1=[item[1] for item in res]
    toberemoved=uniq(l1)
    toberemoved.sort(reverse=True)
    for i in xrange(len(l0)):
        term_names[l0[i]]+=' & '+term_names[l1[i]]
    for gene in toberemoved:
        T.pop(gene)
        large_goids.pop(gene)
        term_names.pop(gene)
    if output>0:
        "New combined nodes:\n"
        for name in term_names:
            if name.find('&')!=-1:
                print name
    return (T,term_names,large_goids)

def remove_genes_at_most_annotated_to_root(unique_genes,T,dict_genes,kind=1):
    #print 'remove genes at most annotated to the root: biological process...'
    lug=len(unique_genes)
    if kind==1:
        toberemoved=set(range(lug))-set.union(*[set([])] + [set(t) for t in T if len(t)<len(unique_genes)])
    else:
        toberemoved=set(range(lug))-set.union(*[set([])] + [set(t) for t in T])

    toberemoved=list(toberemoved)
    toberemoved.sort()

    if toberemoved==[]:
        return (unique_genes,T)
    if kind==1:
        cc=[len(T[i]) for i in xrange(len(T))]
        #remove genes that are only annotated by biological process 
        try:
            bp_index=cc.index(len(unique_genes))
            for gene in toberemoved:
                T[bp_index].remove(gene)
        except ValueError:
            pass
    toberemoved.append(lug)
    minus=[0 for j in xrange(toberemoved[0])]
    for i in xrange(1,len(toberemoved)):
        minus.extend([i for j in xrange(toberemoved[i]-toberemoved[i-1])])
    
    for j in xrange(len(T)):
        #T[j].sort()
        for i in xrange(len(T[j])):
            T[j][i]-=minus[T[j][i]]
    unique_genes=[unique_genes[i] for i in set(range(lug))-set(toberemoved)]
    return (unique_genes,T)
    
def jaccard(liste,T,term_names,output=1):
    """Returns a table, in which all terms in liste are pairwise compared w.r.t their Jaccard index"""
    res=[[0 for i in xrange(len(liste))] for j in xrange(len(liste))]
    for i in xrange(len(liste)):
        for j in xrange(len(liste)):#i+1,len(liste)):
            if i==j:
                continue
            if T[liste[i]]!=[] or T[liste[j]]!=[]:
                res[i][j]=len(set(T[liste[i]])&set(T[liste[j]]))*1./len(set(T[liste[i]])|set(T[liste[j]]))
            else:
                res[i][j]=0
            if output>0:
                print term_names[liste[i]],'&',term_names[liste[j]],len(set(T[liste[i]])&set(T[liste[j]]))*1./len(set(T[liste[i]])|set(T[liste[j]]))
    #print [len(T[c]) for c in liste]
    if len(liste)==2:
        return res[0][1]
    else:
        return res

f = open('human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv')
       
reader = csv.reader(f)

genes=[]
goids=[]
names=[]
for row in reader:
    a=row[0]
    a=a.split('\t')
    genes.append(a[0])
    goids.append(a[1])
    names.append(a[2])
f.close()
unique_goids=uniq(goids)
dict_goids=dict((unique_goids[i],i) for i in xrange(len(unique_goids)))
counter=[0 for i in xrange(len(unique_goids))]

print "Create unique_genes..."               
unique_genes=[]#slower: self.uniq(genes)
current=None
for i in xrange(len(genes)):
    if genes[i]!=current:
        current=None
    if current==None:
        unique_genes.append(genes[i])
        current=genes[i]
unique_genes.sort()    
dict_genes=dict((unique_genes[i],i) for i in xrange(len(unique_genes)))

print "Create large_goids..."
for i in xrange(len(genes)):
    #check whether one gene is not annotated by the same GO term multiple times (different evidence codes)
    try: #check whether gene actually is part of the considered universe of genes
        dict_genes[genes[i]]
    except KeyError:
        continue
    if (i>0 and goids[i]==goids[i-1] and genes[i]==genes[i-1]):
        continue
    else: #check whether one gene is not annotated by the same GO term multiple times
        index=dict_goids[goids[i]]#slower: index= unique_goids.index(goids[i])
        counter[index]+=1
large_goids=[unique_goids[i] for i in xrange(len(counter)) if counter[i]>=1]
dict_goids=dict((large_goids[i],i) for i in xrange(len(large_goids)))

print "Create term_names..."
term_names=[]
for goid in large_goids:
    term_names.append(names[goids.index(goid)])
    
print "Create T..."
T=[[] for i in xrange(len(large_goids))]

for i in xrange(len(genes)):
    try:
        index_goid=dict_goids[goids[i]]
        index_gene=dict_genes[genes[i]]

        if index_gene not in T[index_goid]:
            T[index_goid].append(index_gene)
    except ValueError:
        pass
    except KeyError:
        pass
    if i%100000==0:
        print round(100.*i/len(genes),0),'% finished'

#Check for genes only annotated to root, and for equal terms in GO
(unique_genes,T)=remove_genes_at_most_annotated_to_root(unique_genes,T,dict_genes)
(T,term_names,large_goids)=remove_equal_terms(T,term_names,large_goids)
savevar('T', T)
savevar('term_names', term_names)
savevar('large_goids', large_goids)
savevar('unique_genes', unique_genes)

