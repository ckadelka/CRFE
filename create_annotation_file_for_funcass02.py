'''This program takes as input a gene file and an annotation file (f.e., created by ontology_parser.py).
It parses through the files and creates an annotation file that can be used by FuncAssociate.'''

gene_files=[]
gene_files.extend(['gene_file_GDS3004.txt','gene_file_GDS4419.txt','gene_file_GDS4610.txt','gene_file_GDS4974.txt'])
gene_files.extend(['gene_file_GDS1686.txt','gene_file_GDS2969.txt','gene_file_GDS3216.txt','gene_file_GDS3246.txt','gene_file_GDS3866.txt','gene_file_GDS3928.txt','gene_file_GDS3933.txt','gene_file_GDS4423.txt','gene_file_GDS4929.txt','gene_file_GDS5092.txt'])

for gene_file in gene_files:
    
    #gene_file='datasets/gene_file_GDS3004.txt'
    #annotation_file='ONTOLOGY_PARSER/annotation_file_human_association_human_biological_process.txt'

    if int(gene_file[13:17]) in [3004,4419,4610,4974]:
        annotation_file='ONTOLOGY_PARSER/annotation_file_human_combined_association_human_biological_process.txt'
    elif int(gene_file[13:17]) in [2969,3866]:
        annotation_file='ONTOLOGY_PARSER/annotation_file_yeast_combined_association_yeast_biological_process.txt'
    elif int(gene_file[13:17]) in [1686]:
        annotation_file='ONTOLOGY_PARSER/annotation_file_fly_combined_association_fly_biological_process.txt'
    elif int(gene_file[13:17]) in [3928]:
        annotation_file='ONTOLOGY_PARSER/annotation_file_rat_combined_association_rat_biological_process.txt'
    elif int(gene_file[13:17]) in [3246,4423,4929,5092]:
        annotation_file='ONTOLOGY_PARSER/annotation_file_mouse_combined_association_mouse_biological_process.txt'
    elif int(gene_file[13:17]) in [3216,3933]:
        annotation_file='ONTOLOGY_PARSER/annotation_file_arabidopsis_combined_association_arabidopsis_biological_process.txt'    
    gene_file='datasets/'+gene_file


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
            dict_gene[row[index_gene_in_gene_file].lower()] #if this works, the gene is a duplicate
            continue
        except KeyError:
            pass
        try:
            unique_genes.append(row[index_gene_in_gene_file].lower())
            dict_gene.update({row[index_gene_in_gene_file].lower():i})
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
        genes=row[index_genes_in_annotation_file].lower().split(gene_delimiter_in_annotation_file)
        if len(genes)<cutoff: #don't want terms with less than 'cutoff' gene annotations
            continue
        annotation_counter=0
        term_annotations.append([])
        already_added=[False for i in range(len(unique_genes))]
        for gene in genes:
            try:
                index_gene_in_unique_genes=dict_gene[gene]
                if already_added[index_gene_in_unique_genes]==False:
                    term_annotations[-1].append(gene)
                    already_added[index_gene_in_unique_genes]=True
            except KeyError:
                continue
            annotation_counter+=1
            if annotation_counter>top_cutoff and top_cutoff>=cutoff: #too many gene annotations
                break
        if annotation_counter>=cutoff and (annotation_counter<=top_cutoff or top_cutoff<cutoff):
            term_names.append(row[index_name_in_annotation_file])
            term_ids.append(row[index_id_in_annotation_file])
        else:
            term_annotations.pop()
    reader.close()

    tt=[len(t) for t in term_annotations]
    ind=sorted(range(len(tt)), key=lambda k: tt[k])
    tt_ord=[tt[ind[i]] for i in range(len(tt))]
    res=[]
    for i in range(len(tt_ord)):
        for j in range(i+1,len(tt_ord)):
            if tt_ord[i]<tt_ord[j]:
                break
            elif set(term_annotations[ind[i]])==set(term_annotations[ind[j]]): #set comparison
                res.append((ind[i],ind[j],tt_ord[i]))
    l0=[item[0] for item in res]
    l1=[item[1] for item in res]
    toberemoved=list(set(l1))
    toberemoved.sort(reverse=True)
    #for i in range(len(l0)):
    #    self.term_names[l0[i]]+=' & '+self.term_names[l1[i]]
    for term in toberemoved:
        term_annotations.pop(term)
        term_names.pop(term)
        term_ids.pop(term)
        
    print gene_file,len(toberemoved)
            
    writer=open('annotation_file_new_funcass_'+gene_file.split('_')[-1].split('.')[0]+'.txt','w')
    text=""
    for i in range(len(term_names)):
        text+=term_ids[i]+'\t'+term_names[i]+'\t'
        for gene in term_annotations[i]:
            text+=gene+gene_delimiter_in_annotation_file
        text=text[:-1]+'\n'
    text+="\t\t"
    for gene in set(unique_genes)-set.union(*[set([])] + [term_annotations[i] for i in range(len(term_annotations))]):
        text+=gene+gene_delimiter_in_annotation_file
    text=text[:-1]+'\n'
    writer.write(text)
    writer.close()