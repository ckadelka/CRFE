import sys
import csv
import GeneSetRankedEnrichment31o_web as tools

#gene_file=sys.argv[1]
#min_level=sys.argv[2]
#cutoff=sys.argv[3]
#top_cutoff=sys.argv[4]
#kind=sys.argv[5]
#go_file=sys.argv[6]

gene_file="2-cell-hep-d12_versus_HM-d12_flipped.txt"#"data/bkz_gene_list.csv"#"3-cell-hep-d12_versus_HM-d12.txt"
min_level="0.1"
cutoff="5"
top_cutoff="200"
kind="3"
go_file='human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv'


outstr=("HDF" if gene_file.split('/')[-1]=="bkz_gene_list.csv" else "LIVER")
add2=''
if outstr=="LIVER" and go_file.split('/')[-1]=='human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv':
    add2='_ID'
gene_file=gene_file[:-4]+add2+gene_file[-4:]
if gene_file[-14:-4] == "flipped_ID":
    gene_file = gene_file[:-14]+"ID_flipped.txt"

m = tools.GeneSetRankedEnrichment(2, 1, int(cutoff), int(top_cutoff), 0.1, 
                            0.25, 0.01, 5, 0, float(min_level), 1000, 1000, 20000, int(1e5), 0.5,0.2,19,20,gene_file,
                            go_file,
                            'output/', outstr)

genes_list=[]
level_list=[]
           
gf = open(m.gene_file)

if m.gene_file[-3:]=='csv':
    reader = csv.reader(gf, delimiter='\t')
else:
    reader = open(m.gene_file,'r')

for row in reader:
    if m.gene_file[-3:]=='txt':
        row=row.split('\t')
    try:
        level_list.append(float(row[1]))
        genes_list.append(row[0])
    except Exception:
        continue 

gf.close()

#print str(sum([l<0 for l in level_list]))+',0,-2'

unique_genes=m.loadvar('unique_genes')
if unique_genes==None:
    (_,unique_genes,_,_)=m.load_all_new(genes_list)
    m.savevar('unique_genes', unique_genes)

#print str(len(unique_genes))+","+str(len(level_list))+","+str(len(genes_list))

if kind=="1":
    counter=0
    for i in xrange(len(level_list)):
        if level_list[i]>=float(min_level):
            counter+=1
        else:
            break
    counter=len(set(unique_genes)&set(genes_list[:counter]))
    print min_level+","+str(counter)+","+str(round(counter*1./len(unique_genes),4))
elif kind=="2":
    counter=0
    BROKE=0
    for i in xrange(len(genes_list)):
        if genes_list[i] in unique_genes:
            counter+=1
            if counter>=int(min_level):
                broke_at=i
                BROKE=1
                break
    if BROKE==0:
        broke_at=len(genes_list)-1
        min_level=str(len(unique_genes))
    counter=len(set(unique_genes)&set(genes_list[:counter]))
    print str(level_list[broke_at])+","+min_level+","+str(round(float(min_level)/len(unique_genes),4))
elif kind=="3":
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
    print str(level_list[broke_at])+","+str(nr_active)+","+min_level