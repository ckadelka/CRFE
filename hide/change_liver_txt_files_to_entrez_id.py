import os
import creates_files_for_funcassociate as trans
for files in os.listdir("."):
    if files.endswith(".txt") and 'versus' in files:
        reader=open(files,'r')
        writer=open(files[:-4]+'_ID'+'.txt','w')
        c=0
        for row in reader:
            row=row.split('\t')
            try:
                writer.write(trans.dict_gene[row[0]]+'\t'+row[1])
            except Exception:
                c+=1
                continue