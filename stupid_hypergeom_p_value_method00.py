import math

def hypergeometric_pmf(x, m, n, k):
    """
Given a population consisting of `m` items of class M and `n` items of class N,
this returns the probability of observing `x` items of class M when sampling
`k` times without replacement from the entire population (i.e., {M,N})
p(x) = (choose(m, x) * choose(n, k-x)) / choose(m+n, k)
"""
    a = math.log(binomial_coefficient(m, x))
    b = math.log(binomial_coefficient(n, k-x))
    c = math.log(binomial_coefficient(m+n, k))
    return math.exp(a+b-c)
      
def binomial_coefficient(population, sample):
    "Returns `population` choose `sample`."
    s = max(sample, population - sample)
    assert s <= population
    assert population > -1
    if s == population:
        return 1
    numerator = 1
    denominator = 1
    for i in range(s+1, population + 1):
        numerator *= i
        denominator *= (i - s)
    return numerator/denominator

ranked_list=[0 for i in range(len(m.unique_genes))]
count=0
for i in range(len(m.genes_list)):
    try:
        ranked_list[m.unique_genes.index(m.genes_list[i])]=count
        count+=1
    except ValueError:
        continue

Ns,Ms,Xs,Ps=[],[],[],[]
nr_genes=len(m.unique_genes)
for t in m.T:
    ranks=[ranked_list[gene] for gene in t]
    ranks.sort()
    len_term=len(t)
    p=1
    for el in [Ns,Ms,Xs,Ps]:
        el.append([0,0,0,0])
    for i,rank in enumerate(ranks):
        dummy=hypergeometric_pmf(i+1, len_term, nr_genes-len_term, rank+1)
        if dummy<p:
            p=dummy
            Ns[-1]=(i+1)
            Ms[-1]=(rank)
            Xs[-1]=(len_term)
            Ps[-1]=(dummy)
    print m.T.index(t)
            
            
    


m.glist=[m.genes_list.index(m.unique_genes[i]) for i in xrange(len(m.unique_genes))]
m.glist_inv=[m.unique_genes.index(m.genes_list[i]) for i in xrange(len(m.genes_list))]
