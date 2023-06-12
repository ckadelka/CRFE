#This file creates a dictionary num2symbol, which maps gene symbol to entrez ID
import csv

def dict_transform(gene_file='3-cell-hep-d12_versus_2-cell-hep-d12.txt'):
    g=open('list_of_entrez_genes','w')
    f=open(gene_file,'r')
    if gene_file[-3:]=='csv':
        reader = csv.reader(f, delimiter='\t')
    elif gene_file[-3:]=='txt':
        reader = f.read().splitlines()
    
    genes=[]
    values=[]
    for row in reader:
        if gene_file[-3:]=='txt':
            row=row.split('\t')
        try:
            if float(row[1])<-3:
                break
        except ValueError:
            continue
        #if int(row[0]) in unique_genes:
        genes.append(row[0])
        values.append(row[1])
        #if '///' in row[0]:
        #    g.write(row[0][:row[0].index('///')-1]+"\n")
        #else:
        #    g.write(row[0]+"\n")
        
    for gene,i in zip(genes,range(len(genes))):
        #if '///' in gene:
        #    g1,g2=gene.split(' /// ')
        #    if g1 not in genes:
        #        genes[i]=g1
        #    elif g2 not in genes:
        #        genes[i]=g2
        #    else:
        #        genes[i]='duplicate entry:'+g1+','+g2
        #g.write(genes[i]+"\n")
        #g.write(genes[i]+"\t"+values[i]+"\n")
        if '///' in gene:
            g1,g2=gene.split(' /// ')
            if g1 not in genes:
                genes[i]=g1
            elif g2 not in genes:
                genes[i]=g2
            else:
                genes[i]='duplicate entry:'+g1+','+g2
        g.write(genes[i]+"\n")
    g.close()
    
    
    ##Build gene ontology in format required by FuncAssociate2.0
    #(T, unique_genes, large_goids, term_names)=m.load_all_new('bkz_gene_list.csv')
    #text=""
    #f=open('for_gene_associate.txt','w')
    #for i in xrange(len(T)):
    #    text=str(large_goids[i])+"\t"+term_names[i]+"\t"
    #    for t in T[i]:
    #        text+=str(t)+" "
    #    text=text[:-1]
    #    f.write(text+'\n')
    #f.close()
    
    #f=open('EntrezGene_IDS-transfromed_by_MADgene_LIver.txt','r')
    #g=open('EntrezGene_IDS-transfromed_back_by_MADgene_LIver.txt','w')
    #reader=f.read().splitlines()
    #for row in reader:
    #    g.write(row+"\n")
    #g.close()
    #f.close()
    
    ##Open the transformed file (by Madgene & Synergizer) to optimize the dictionary between GeneSymbol and EntrezIDf=open('EntrezGene_IDS-transfromed_by_MADgene_LIver.txt','r')
    dict_gene={}
    def add_entry_if_int(dict_gene,key,value):
        try:
            dict_gene[key]=str(int(value))
        except:
            pass
        return dict_gene
    
    f=open('MADGENE_SYNERGIZER.csv','rU')
    reader = csv.reader(f, delimiter='\t')
    count=-1
        
    for row in reader:
        count+=1
        if count<3:
            continue
        row=row[0].split(',')
        if row[1]==row[3]:
            dict_gene=add_entry_if_int(dict_gene,row[0],row[1])
        elif row[1] in row[3]:
            dict_gene=add_entry_if_int(dict_gene,row[0],row[1])
        elif row[3] in row[1]:
            dict_gene=add_entry_if_int(dict_gene,row[0],row[3])
        elif row[1]=='Input not found' and row[3]!='I?':
            dict_gene=add_entry_if_int(dict_gene,row[0],row[3])
        elif row[1]!='Input not found':
            dict_gene=add_entry_if_int(dict_gene,row[0],row[1])
    
    #for row in reader:
    #    count+=1
    #    if count<3:
    #        continue
    #    row=row[0].split(',')
    #    if row[1]==row[3]:
    #        dict_gene[row[0]]=row[1]
    #    elif row[1] in row[3]:
    #        dict_gene[row[0]]=row[1]
    #    elif row[3] in row[1]:
    #        dict_gene[row[0]]=row[3]
    #    elif row[1]=='Input not found' and row[3]!='I?':
    #        dict_gene[row[0]]=row[3]
    #    elif row[1]!='Input not found':
    #        dict_gene[row[0]]=row[1]
    
    symbol2num=dict_gene
    num2symbol=dict((value, key) for key, value in dict_gene.items())
    return (symbol2num,num2symbol)