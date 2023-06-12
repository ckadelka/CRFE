import crfe_v007 as tools

import numpy as np
import sys
if sys.version_info[0] < 3:
    import cPickle as pickle
else:
    import pickle
import optparse

def main():
    p = optparse.OptionParser('DESCRIPTION')
    p.add_option('--nsim', '-n', default='20', help='Number of simulations.')
    p.add_option('--which', '-w', default='10', help='Which dataset should be analyzed.')
    p.add_option('--belief', '-b', default='5', help='Belief parameter.')
    p.add_option('--nr_categories', '-s', default='0', help='Number of categories.')
    
    options, arguments = p.parse_args() 
       
    try:
        nsim = int(options.nsim)
    except ValueError:
        nsim = 20
    try:
        which = int(options.which)
    except ValueError:
        which = 12
    try:
        belief = int(options.belief)
    except ValueError:
        belief = 5
    try:
        s = int(options.nr_categories)
    except ValueError:
        s = 0                     
    
    nsim = 1
    repeats = 1
    cutoff=20
    top_cutoff=200
    burnin=int(1*1e5)
    steps=int(1*1e6)
    alpha_beta_max=0.5
    
    gene_files=['bkz_gene_list.csv','3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12.txt','datasets/disease_glaucoma.txt','datasets/disease_obesity.txt','3-cell-hep-d12_versus_CS-d12_ID.txt','datasets/disease_adenocarcinoma_of_lung.txt','datasets/disease_ischemic_cardiomyopathy.txt','3-cell-hep-d12_versus_2-cell-hep-d12_pos_ID.txt','3-cell-hep-d12_versus_2-cell-hep-d12_neg_ID.txt']
    gene_files.extend(['gene_file_GDS3004.txt','gene_file_GDS4419.txt','gene_file_GDS4610.txt','gene_file_GDS4974.txt'])
    gene_files.extend(['gene_file_GDS1686.txt','gene_file_GDS2969.txt','gene_file_GDS3216.txt','gene_file_GDS3246.txt','gene_file_GDS3866.txt','gene_file_GDS3928.txt','gene_file_GDS3933.txt','gene_file_GDS4423.txt','gene_file_GDS4929.txt','gene_file_GDS5092.txt'])
    out_strs=['HDF','LIVER','LIVER','DISEASE','DISEASE','LIVER','DISEASE','DISEASE','Lpos','Lneg','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI','NCBI']
    for hdhsfs in [which]:
        for alpha_beta_max in [0.5]:
            m = tools.CRFE(repeats, 2, cutoff, top_cutoff, belief, 0.3, 'proportion', burnin, steps, alpha_beta_max,0.2,20,20,gene_files[hdhsfs],
                                            'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                            'output1/', 'new2019_'+out_strs[hdhsfs])

            if hdhsfs==2:
                m.annotations_file='msigdb.txt'
            elif 'NCBI' in m.out_str:
                if int(m.gene_file[13:17]) in [3004,4419,4610,4974]:
                    m.annotations_file='ONTOLOGY_PARSER/annotation_file_human_combined_association_human_biological_process.txt'
                elif int(m.gene_file[13:17]) in [2969,3866]:
                    m.annotations_file='ONTOLOGY_PARSER/annotation_file_yeast_combined_association_yeast_biological_process.txt'
                elif int(m.gene_file[13:17]) in [1686]:
                    m.annotations_file='ONTOLOGY_PARSER/annotation_file_fly_combined_association_fly_biological_process.txt'
                elif int(m.gene_file[13:17]) in [3928]:
                    m.annotations_file='ONTOLOGY_PARSER/annotation_file_rat_combined_association_rat_biological_process.txt'
                elif int(m.gene_file[13:17]) in [3246,4423,4929,5092]:
                    m.annotations_file='ONTOLOGY_PARSER/annotation_file_mouse_combined_association_mouse_biological_process.txt'
                elif int(m.gene_file[13:17]) in [3216,3933]:
                    m.annotations_file='ONTOLOGY_PARSER/annotation_file_arabidopsis_combined_association_arabidopsis_biological_process.txt'
                m.gene_file='datasets/'+m.gene_file
                m.out_str+='_'+m.gene_file.split('_')[-1][:-4]

            #data collection
            C=[[]]*nsim
            MCMC_distr=[[]]*nsim
            STD_distr = [[]]*nsim
            for i in range(nsim):
                print(i)
                (C[i],MCMC_distr[i],STD_distr[i])=m.runMe(2)
            
            print(len(m.T),len(m.unique_genes),sum(m.sets_selected),m.nr_categories)

            
            if m.out_str=='DISEASE':
                m.out_str=m.gene_file.split('.')[0].split('datasets/disease_')[1]
        #save data
#            addon='_alpha_beta_max'+str(int(round(m.alpha_beta_max*1000)))
#            if m.annotations_file=='msigdb.txt':
#                f=open('saved_data/multiple_runs_new_msigdb_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','w+')
#            else:
#                f=open('saved_data/multiple_runs_new_go_'+m.out_str+'_s'+str(m.nr_categories)+'_cutoff'+str(m.cutoff)+'_topcutoff'+str(m.top_cutoff)+'_wk'+str(m.weight_kind)+'_ck'+str(m.category_kind)+'_belief'+str(m.belief)+'_burnin'+str(m.burnout)+'_steps'+str(m.MCMC_steps)+'_percentactive'+str(min_level)+addon+'_nsim'+str(nsim)+'.txt','w+')
#            pickle.dump([C,MCMC_distr], f)
#            f.close()
    return m

if __name__ == '__main__':
    m = main()
    
    
    
    upto = 30000
    x_time = 1000*(np.array(m.timer[1:(upto+1)])-np.array(m.timer)[:upto])
    x_which = np.array(m.which[:upto])
    x_accepted = np.array(m.accepted[:upto])
    for i in set(x_which):
        for j in set(x_accepted):
            indices = np.bitwise_and(x_which==i,x_accepted==j)
            print(i,j,sum(indices),np.mean(x_time[indices]))

