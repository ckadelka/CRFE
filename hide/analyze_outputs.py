import pylab as P
import numpy as np
import csv, math

"""Analyses the data generated from GeneSetRankedEnrichment.save_simulation_results"""

def str_to_vec(s):
    s=s.replace(" ","")
    if s=="[]":
        return []
    if s[0]=="[":
        s=s[1:]
    if s[-1]=="]":
        s=s[:-1]
    
    vec=[]
    for i in s.split(','):
        try:
            vec.append(int(i))
        except:
            vec.append(float(i))
    return vec

def str_to_mat(s):
    if s[0]=="[":
        s=s[1:]
    if s[-1]=="]":
        s=s[:-1]
    mat=[]
    s=s.replace(" ","")
    for i in s.split('],'):
        mat.append(str_to_vec(i))
    return mat

def analyze(filename,variable=[3],fixed=['HDF',0.6,'GenGO',2,5,200,5,0.9,0.01,3,0.01,0.15,0.01,1000,1000,0,0,0],CONSIDER_DOUBLE_ENTRIES=False):
    gf = open(filename)

    if filename[-3:]=='csv':
        reader = csv.reader(gf, delimiter='\t')
    else:
        reader = open(filename,'r')

    variable.sort()

    res=[]
    
    for row in reader:
        if filename[-3:]=='txt':
            row=row.split('\t')

        values_vars=[]
        ACCEPT=True
        #check whether this row should be analyzed
        for i in xrange(18):
            if row[2]=='Bayesian' and i in [7,8,9]:
                continue
            elif row[2]=='GenGO' and i in [10,11,12]:
                continue
            try:
                row[i]=int(row[i])
            except:
                try:
                    row[i]=float(row[i])
                except:
                    pass
            if i in variable:
                values_vars.append(row[i])
                continue
            elif row[i]!=fixed[i]:
                ACCEPT=False
                break

        if CONSIDER_DOUBLE_ENTRIES==False:
            for j in xrange(len(res)): #no double entries considered
                if sum([row[variable[i]]==res[j][12+i] for i in xrange(len(variable))])==len(variable):
                    ACCEPT=False
                    continue
        
        if ACCEPT==False:
            continue
        
        res.append([int(row[3]),int(row[18]),str_to_vec(row[19]),
                    int(row[20]),int(row[21]),int(row[22]),
                    str_to_vec(row[23]),float(row[24]),
                    str_to_vec(row[25]),str_to_vec(row[26]),
                    str_to_vec(row[27]),str_to_mat(row[28])])
        res[-1].extend(values_vars)

        #0: s, 1: lenC, 2: C, 3: sumC, 4: lenT, 5:total_annotated,
        #6: p-values, 7: L-value, 8: mean-gene level, 9: mean(expl)
        #10: mean(annotations), 11: max(jaccard comp)
    return res

#analyze jaccard comparison
def jaccard_comparison(res):
    mean_jaccard_comp=[]

    for cat in range(3):
        mean_jaccard_comp.append([])
        for item in res:
            if cat==0 and item[0]<2:
                ind=0
            elif cat==0 and item[0]==2 or cat==2 and item[0]>2:
                ind=3
            elif cat==0 and item[0]>2 or cat==1 and item[0]==2 or cat==2 and item[0]<2:
                ind=2
            elif cat==1 and item[0]!=2:
                ind=1
            elif cat==2 and item[0]==2:
                ind=6
            mean_jaccard_comp[-1].append(sum([item[11][ind][k]*k for k in range(20)]))

    return mean_jaccard_comp

a = analyze('simulation_results.txt',variable=[3,5])
