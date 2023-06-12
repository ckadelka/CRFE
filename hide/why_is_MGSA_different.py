f=open('MGSA_ontology_Liver_GO_20-200.txt','r')

lines=f.read().splitlines()
f.close()

TTT,ccc=[],[]

for line in lines:
    a=line.split('\t')
    TTT.append(a[1])
    ccc.append(len(a[2].split(" ")))