'''This program takes as input a gene file and an annotation file (f.e., created by ontology_parser.py).
It parses through the files and creates an annotation file that can be used by FuncAssociate.'''

gene_file='datasets/gene_file_GDS3004.txt'
annotation_file='ONTOLOGY_PARSER/annotation_file_human_association_human_biological_process.txt'
cutoff=20
top_cutoff=200

index_gene_in_gene_file=0
index_id_in_annotation_file=0
index_name_in_annotation_file=1
index_genes_in_annotation_file=2
gene_delimiter_in_annotation_file=' '

reader=open(gene_file,'r')
unique_genes=[]
dict_gene={}
i=0
for row in reader:
    row=row.split('\t')
    try:
        unique_genes.append(row[index_gene_in_gene_file])
        dict_gene.update({row[index_gene_in_gene_file]:i})
        i+=1
    except:
        pass
reader.close()

reader=open(annotation_file,'r')
T=[]
term_names=[]
term_ids=[]
term_annotations=[]
for row in reader:
    row=row.split('\t')
    genes=row[index_genes_in_annotation_file].split(gene_delimiter_in_annotation_file)
    if len(genes)<cutoff: #don't want terms with less than 'cutoff' gene annotations
        continue
    annotation_counter=0
    term_annotations.append([])
    for gene in genes:
        try:
            index_gene_in_unique_genes=dict_gene[gene]
            term_annotations[-1].append(gene)
        except KeyError:
            continue
        annotation_counter+=1
    if annotation_counter>=cutoff and (annotation_counter<=top_cutoff or top_cutoff<cutoff):
        term_names.append(row[index_name_in_annotation_file])
        term_ids.append(row[index_id_in_annotation_file])
    else:
        term_annotations.pop()
reader.close()

writer=open('annotation_file_funcass_'+gene_file.split('_')[-1].split('.')[0]+'.txt','w')
text=""
for i in range(len(term_names)):
    text+=term_ids[i]+'\t'+term_names[i]+'\t'
    for gene in term_annotations[i]:
        text+=gene+gene_delimiter_in_annotation_file
    text=text[:-1]+'\n'
writer.write(text)
writer.close()