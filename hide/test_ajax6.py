import sys
import GeneSetRankedEnrichment31q_web as tools

##gene_file=sys.argv[1]
##min_level=sys.argv[2]
##cutoff=sys.argv[3]
##top_cutoff=sys.argv[4]
##kind=sys.argv[5]
##go_file=sys.argv[6]
##kind_of_file=sys.argv[7]
#
#gene_file="bkz_gene_list.csv"#3-cell-hep-d12_versus_HM-d12.txt"#"bkz_gene_list.csv"#"3-cell-hep-d12_versus_HM-d12.txt"
#min_level="0.2"
#cutoff="5"
#top_cutoff="200"
#kind="3"
#go_file='human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv'
#kind_of_file="0"
#
#outstr=("HDF" if gene_file.split('/')[-1]=="bkz_gene_list.csv" else "LIVER")
#add2=''
#if outstr=="LIVER" and go_file.split('/')[-1]=='human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv':
#    add2='_ID'
#
#gene_file=gene_file[:-4]+add2+gene_file[-4:]
#
#m = tools.GeneSetRankedEnrichment(2, 1, int(cutoff), int(top_cutoff), 0.1, 0.25, 0.01,
#                            5, 0, float(min_level), 1000, 1000, 20000, int(1e5), 0.5,0.2,19,20,
#                            gene_file,int(kind_of_file),go_file,'output/', outstr)

#(genes_list, level_list) = m.get_sets_and_genes()
#(T, unique_genes, large_goids, term_names) = m.load_all_new(genes_list)
#lug=len(unique_genes)
#Tset=set(range(len(T)))
#(G,_, dict_G)=m.getG(unique_genes,genes_list,T,Tset,level_list)

##print str(sum([l<0 for l in level_list]))+',0,-2'
#
#unique_genes=m.loadvar('unique_genes')
#if unique_genes==None:
#    (_,unique_genes,_,_)=m.load_all_new(genes_list)
#    m.savevar('unique_genes', unique_genes)

#print str(len(unique_genes))+","+str(len(level_list))+","+str(len(genes_list))

#print genes_list[0],level_list[0],genes_list[-1],level_list[-1]

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