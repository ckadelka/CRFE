#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 17:21:23 2023

@author: ckadelka
"""

import pandas as pd
import numpy as np
import os

folder = 'data'

for root,subdirs,files in os.walk('data'):
    #print(root,subdirs,files)
    if root == folder:
        continue
    for filename in os.listdir(root):
        if filename.endswith('analytics.tsv'):
            print(filename)
            A = pd.read_csv(root+'/'+filename,delimiter=('\t'))
            columns = list(A.columns)
            indices_p_columns = [i for i,el in enumerate(columns) if 'p-value' in el]
            indices_foldchange_columns = [i for i,el in enumerate(columns) if 'foldchange' in el]
            names_contrasts = ['.'.join(el.split('.')[:-1]) for i,el in enumerate(columns) if 'p-value' in el]
            assert names_contrasts == ['.'.join(el.split('.')[:-1]) for i,el in enumerate(columns) if 'foldchange' in el], 'ensure the order of the contrasts is the same for p-values and foldchanges'
            GENE_HAS_SYMBOL = [1 if type(el)==str else 0 for el in A['Gene Name']]
            for index_p,index_fc,name_contrast in zip(indices_p_columns,indices_foldchange_columns,names_contrasts):
                P_VALUE_NOT_NAN = [1 if not np.isnan(el) else 0 for el in A[columns[index_p]]]
                which = np.bitwise_and(GENE_HAS_SYMBOL,P_VALUE_NOT_NAN)
                which = np.array(which,dtype=bool)
                genes = np.array(list(A['Gene Name']))[which]
                foldchanges = np.array(list(A[columns[index_fc]]))[which]
                p_values = np.array(list(A[columns[index_p]]))[which]
                values = np.array(foldchanges,dtype=np.float64) + 0.01* np.array([1-p if change >= 0 else p for change,p in zip(foldchanges,p_values)])
                
                #sort
                indices = np.argsort(values)[::-1]
                
                
                #write to file
                new_filename = root.split('/')[-1] + '_' + name_contrast
                f = open(folder + '/' + new_filename+'.txt','w')
                for index in indices:
                    f.write(genes[index] + '\t' + str(values[index]) + '\n')
                f.close()