import GeneSetRankedEnrichment31q_web as tools
import sys

gene_file=sys.argv[1]
cutoff=sys.argv[2]
top_cutoff=sys.argv[3]
go_file=sys.argv[4]
kind_of_file=sys.argv[5]

#gene_file="3-cell-hep-d12_versus_HM-d12.txt"#"bkz_gene_list.csv"#"3-cell-hep-d12_versus_HM-d12.txt"
#cutoff="5"
#top_cutoff="200"
#go_file='human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv'
#kind_of_file="0"

outstr=("HDF" if gene_file.split('/')[-1]=="bkz_gene_list.csv" else "LIVER")
add2=''
if outstr=="LIVER" and go_file.split('/')[-1]=='human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv':
    add2='_ID'

gene_file=gene_file[:-4]+add2+gene_file[-4:]

m = tools.GeneSetRankedEnrichment(2, 1, int(cutoff), int(top_cutoff), 0.1, 0.25, 0.01,
                            5, 0, 0.5, 1000, 1000, 20000, int(1e5), 0.5,0.2,19,20,
                            gene_file,int(kind_of_file),go_file,'output/', outstr)

(genes_list, level_list) = m.get_sets_and_genes()
(T, unique_genes, large_goids, term_names) = m.load_all_new(genes_list)

gg=[0 for g in unique_genes]
for t in T:
    for gene in t:
        gg[gene]+=1
hist_steps=max(gg)+1

(yy,xx)=tools.hist(gg,min(100,hist_steps),0)

jsondata="{\"data\": ["
for i in range(len(xx)):
   jsondata += "["+str(xx[i])+", "+str(yy[i])+"],"
jsondata=jsondata[:-1]
jsondata += "]}"

print jsondata