import GeneSetRankedEnrichment31q_web as tools

top_cutoff=0
cutoff=1

f=open('gene2go','r')
count=0
genes=[]
goids=[]
names=[]
genes2=[]
goids2=[]
names2=[]
for row in f:
    a=row.split('\t')
    if len(a)>1 and a[0]=='10116' and a[7]=='Process\n':
        genes.append(a[1])
        goids.append(a[2][3:])
        names.append(a[5])
    elif len(a)>1 and a[0]=='9606' and a[7]=='Process\n':
        genes2.append(a[1])
        goids2.append(a[2][3:])
        names2.append(a[5])
    count+=1
    
f.close()

unique_goids=tools.GeneSetRankedEnrichment.uniq(goids)
dict_goids=dict((unique_goids[i],i) for i in xrange(len(unique_goids)))
counter=[0 for i in xrange(len(unique_goids))]
#print "Create unique_genes..."               
#unique list of genes
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
        
#print "Create large_goids..."
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
if top_cutoff>cutoff:
    large_goids=[unique_goids[i] for i in xrange(len(counter)) if top_cutoff>=counter[i]>=cutoff]
else:
    large_goids=[unique_goids[i] for i in xrange(len(counter)) if counter[i]>=cutoff]
dict_goids=dict((large_goids[i],i) for i in xrange(len(large_goids)))
#print "Create term_names..."
term_names=[]
for goid in large_goids:
    term_names.append(names[goids.index(goid)])
#print "Create T..."
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
cc=[len(t) for t in T]
