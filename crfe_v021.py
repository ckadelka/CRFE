'''This program implements the functional enrichment method CRFE,
described in C. Kadelka, M. Brandon, T. M. Murali, Concise Functional Enrichment
of Ranked Gene Lists. 
See www.math.vt.edu/people/claus89/crfe/ for help.'''

#v021:  cleaned up version, deleted nr_categories from argument list (always nr_categories==0, which leads to 2x speed up)
#v019:  re-added the possibility of specifying alpha and beta and not learning them during the MCMC by setting proportion_parameter_change=0
#       modified alpha_beta_max so that the highest false rate <= 0.5, that is alpha, beta <= (1+b)/(4b)
#v018:  added an additional column (mean position of perturbed genes) in output file
#v017:  corrected an error in initialize(), the belief was always set to 5 in v015&016
#v015:  include belief parameter among the set of learned parameters
#v014:  got rid of self.max_nr_alpha_plus_nr_beta in MCMC, propose_state sped up
#v013:  added additional information about threshold used and proportion perturbed genes in the output file,
#       only terms that show up at least once in C after the burnin period are reported in the output file,
#       the formula for the max value of alpha and beta given belief is corrected
#v011: corrected errors in the calculation of p_values in log_binomial_coefficient

import sys
import csv, math
if sys.version_info[0] < 3:
    import cPickle as pickle
else:
    import pickle
import random
import optparse
import datetime
import os
import time
import numpy as np
import pandas as pd

print(sys.argv[0])

which_python = sys.version_info[0]
            
class CRFE:
    ##############################
    ###  Initializing Methods  ###
    ##############################
    
    def __init__(self, repeats, lower_cutoff, upper_cutoff, belief, 
                 threshold, threshold_type, burnin_steps, MCMC_steps, 
                 alpha_beta_max, proportion_parameter_change, A, P, 
                 gene_file, annotations_file, output_folder, out_str='', 
                 seed=-1, max_belief=10, B=10,
                 LEARN_BELIEF=False,LEARN_PENALIZATION_PARAMETER=True,penalization_parameter=0.001,GET_INFO_ON_CURRENT_SETS=False,alpha=0.1,beta=0.25):        
        self.repeats = repeats
        self.lower_cutoff = lower_cutoff
        self.upper_cutoff = upper_cutoff
        self.belief = belief
        self.belief_alpha = belief
        self.belief_beta = belief
        self.initial_belief = belief
        self.threshold=threshold
        self.threshold_type=threshold_type
        self.burnin_steps = burnin_steps
        self.MCMC_steps = MCMC_steps
        self.alpha_beta_max = alpha_beta_max
        self.proportion_parameter_change=proportion_parameter_change
        self.A = A
        self.P = P
        self.B = B
        self.max_belief = max_belief #set max_belief == belief == 1 if you want to run MGSA as it was published. If max_belief > belief, then the max false rate is lower than it can be
        self.gene_file = gene_file
        self.annotations_file = annotations_file
        self.out_str = out_str
        self.output_folder = output_folder
        self.seed=seed
        self.LEARN_BELIEF = LEARN_BELIEF
        self.LEARN_PENALIZATION_PARAMETER = LEARN_PENALIZATION_PARAMETER
        self.penalization_parameter = penalization_parameter
        self.GET_INFO_ON_CURRENT_SETS = GET_INFO_ON_CURRENT_SETS
        
        self.alpha = alpha #this will be reset to closest allowable alpha by self.initialize if proportion_parameter_change > 0
        self.beta = beta #this will be reset to closest allowable beta by self.initialize if proportion_parameter_change > 0
        
        self.show_genes_in_output=False
        
        try:
            assert self.belief <= self.max_belief
        except AssertionError:
            print('\n!!! Warning:\nmax_belief must be as large or larger than belief. Now, max_belief == belief == %f)' % self.belief)
            self.max_belief = self.belief
            
        
    @staticmethod
    def uniq(list_with_duplicate_entries):
        temp = set(list_with_duplicate_entries)
        return list(temp)

    @staticmethod
    def discretize(vector,value):
        diff=[math.fabs(v-value) for v in vector]
        ind=sorted(range(len(vector)), reverse=False, key=lambda k: diff[k])
        return (vector[ind[0]],ind[0])

    def get_expression_data(self):
        """load gene_file and create an ordered (ranked) list of genes, and their expression levels"""        
        data = pd.read_csv(self.gene_file,sep='\t',header=None)
        n = data.shape[0]
        if data.shape[1]==1: #add fake levels if none are provided
            data['level'] = np.arange(n-1,-1,-1)/(n-1)
            self.LEVEL_LIST_PROVIDED = False
        else:
            data = data.sort_values(by=[data.columns[1]],ascending=False)
            self.LEVEL_LIST_PROVIDED = True
        self.level_list = np.array(data.iloc[:,1],dtype=float)
        self.genes_list = np.array(data.iloc[:,0])            


    #################################
    ###  Basic Algorirthm Methods ###
    #################################
                     
    def get_specific_rates(self, falserate, number, belief):
        '''This function is used to find the category-specific FPRs (alpha_i) or FNRs (beta_i). 
        Falserate is the desired overall value of alpha or beta'''
        
        falserate0 = 2.*falserate/(belief+1)
        return np.linspace(falserate0,belief*falserate0,number)

    #################################
    ###         MCMC Code         ###
    ################################# 
   
    def propose_MCMC_state(self,neighborhood_size):
        self.term_id_to_toggle, self.term_id_to_remove, self.term_id_to_add = -1,-1,-1
        
        r=random.random()
        if r > self.proportion_parameter_change: #add, delete, or switch a term
            proposal = int(random.random()*neighborhood_size)
        
            if proposal < self.n_terms: #add or delete a term
                self.term_id_to_toggle = proposal
                self.toggle_term(proposal)
                
                if self.LEARN_PENALIZATION_PARAMETER:
                    #learn self.penalization_parameter every time because it is easy
                    self.old_nr_penalization_parameter=self.nr_penalization_parameter
                    self.old_penalization_parameter=self.penalization_parameter
                    self.nr_penalization_parameter=min(max(1,self.number_of_selected_terms),self.P)	
                    self.penalization_parameter=self.nr_penalization_parameter*1./self.n_terms
            else: #switch two terms
                proposal -= self.n_terms
     
                selected_term_pos = proposal // self.number_of_unselected_terms 
                unselected_term_pos = proposal % self.number_of_unselected_terms + self.number_of_selected_terms
                
                self.term_id_to_remove = self.internal_list_of_terms_for_MCMC[selected_term_pos]
                self.term_id_to_add = self.internal_list_of_terms_for_MCMC[unselected_term_pos]
                
                self.remove_term(self.term_id_to_remove)
                self.add_term(self.term_id_to_add)
        else: #change one of the parameters
            if r < self.proportion_parameter_change/3 and self.LEARN_BELIEF or r < self.proportion_parameter_change/2 and self.LEARN_BELIEF==False: #change alpha
                self.old_nr_alpha = self.nr_alpha
                self.old_NP_log_alpha = self.NP_log_alpha
                self.nr_alpha=int(random.random()*self.A)
                self.alpha = self.possible_values_for_alpha_and_beta[self.nr_alpha]
                self.log_alphas=self.possible_log_alphas[self.nr_belief][self.nr_alpha]
                self.NP_log_alpha = np.sum(self.log_alphas[self.n_annotations_per_perturbed_gene==0])
                self.parameter_changed = 'alpha'
            elif self.LEARN_BELIEF==False or r < self.proportion_parameter_change*2/3: #change beta
                self.old_nr_beta = self.nr_beta
                self.old_EP_log_1_minus_beta=self.EP_log_1_minus_beta
                self.nr_beta=int(random.random()*self.A)
                self.beta = self.possible_values_for_alpha_and_beta[self.nr_beta]
                self.log_1_m_betas=self.possible_log_1_minus_betas[self.nr_belief][self.nr_beta]
                self.EP_log_1_minus_beta = np.dot(self.log_1_m_betas,self.n_annotations_per_perturbed_gene)
                self.parameter_changed = 'beta'
            else: #change belief parameter
                self.old_nr_belief = self.nr_belief
                self.old_EP_log_1_minus_beta=self.EP_log_1_minus_beta
                self.old_NP_log_alpha = self.NP_log_alpha
                self.nr_belief=int(random.random()*self.B)
                self.belief = self.possible_values_for_belief[self.nr_belief]
                self.log_alphas=self.possible_log_alphas[self.nr_belief][self.nr_alpha]
                self.log_1_m_betas=self.possible_log_1_minus_betas[self.nr_belief][self.nr_beta]             
                self.NP_log_alpha = np.sum(self.log_alphas[self.n_annotations_per_perturbed_gene==0])#np.dot(self.log_alphas,self.n_perturbed_genes_per_category-self.n_annotations_per_perturbed_gene)                    
                self.EP_log_1_minus_beta = np.dot(self.log_1_m_betas,self.n_annotations_per_perturbed_gene)
                self.parameter_changed = 'belief'
                

    def new_gene_annotated_to_selected_terms(self,gene):
        if gene < self.n_perturbed_genes: #gene is perturbed
            self.n_annotations_per_perturbed_gene[gene] = 1
            self.EP_log_1_minus_beta+=self.log_1_m_betas[gene]
            self.NP_log_alpha-=self.log_alphas[gene]
        else:#gene is unperturbed
            self.EU+=1
            self.NU-=1

    def gene_no_longer_annotated_to_selected_terms(self,gene):           
        if gene < self.n_perturbed_genes: #gene is perturbed
            self.n_annotations_per_perturbed_gene[gene] = 0
            self.EP_log_1_minus_beta-=self.log_1_m_betas[gene]
            self.NP_log_alpha+=self.log_alphas[gene]  
        else: #gene is perturbed
            self.EU-=1
            self.NU+=1


    def add_term(self,term_id_to_add):
        if self.TERM_IS_CURRENTLY_SELECTED[term_id_to_add]:
            return
        self.TERM_IS_CURRENTLY_SELECTED[term_id_to_add] = True
        
        #Go through all genes annotated to the term and increase the annotatation count
        for gene in self.T[term_id_to_add]: 
            self.n_selected_terms_gene_is_annotated_to[gene]+=1
            #check if gene was not annotated to any of the previously selected terms
            if self.n_selected_terms_gene_is_annotated_to[gene]==1:
                self.new_gene_annotated_to_selected_terms(gene)
                
        self.number_of_unselected_terms-=1
        if self.number_of_unselected_terms > 0: #otherwise, already sorted the way we want it	
            #internal sorting of the list that allows to quickly grab terms to be switched in propose_MCMC_state	
            pos = self.position_of_term_in_internal_list[term_id_to_add]	
            e0 = self.internal_list_of_terms_for_MCMC[self.number_of_selected_terms]	
            self.internal_list_of_terms_for_MCMC[pos] = e0	
            self.position_of_term_in_internal_list[e0] = pos	
            self.internal_list_of_terms_for_MCMC[self.number_of_selected_terms] = term_id_to_add	
            self.position_of_term_in_internal_list[term_id_to_add] = self.number_of_selected_terms	
        self.number_of_selected_terms += 1



    def remove_term(self,term_id_to_remove):
        if self.TERM_IS_CURRENTLY_SELECTED[term_id_to_remove]==False:
            return
        self.TERM_IS_CURRENTLY_SELECTED[term_id_to_remove] = False

        #Go through all genes annotated to the term and decrease the annotatation count
        for gene in self.T[term_id_to_remove]:
            self.n_selected_terms_gene_is_annotated_to[gene]-=1
            #check if gene is no longer annotated to any of the currently selected terms
            if (self.n_selected_terms_gene_is_annotated_to[gene] == 0):
                self.gene_no_longer_annotated_to_selected_terms(gene)
        
        self.number_of_selected_terms -= 1	
        if self.number_of_selected_terms > 0: #otherwise, already sorted the way we want it	
            #internal sorting of the list that allows to quickly grab terms to be switched in propose_MCMC_state
            pos = self.position_of_term_in_internal_list[term_id_to_remove]
            e1 = self.internal_list_of_terms_for_MCMC[self.number_of_selected_terms]
            self.internal_list_of_terms_for_MCMC[pos] = e1
            self.position_of_term_in_internal_list[e1] = pos
            self.internal_list_of_terms_for_MCMC[self.number_of_selected_terms] = term_id_to_remove
            self.position_of_term_in_internal_list[term_id_to_remove] = self.number_of_selected_terms
        self.number_of_unselected_terms+=1

    def toggle_term(self,term_id_to_toggle):        
        if self.TERM_IS_CURRENTLY_SELECTED[term_id_to_toggle]:
            self.remove_term(term_id_to_toggle)
        else:
            self.add_term(term_id_to_toggle)

    def undo_proposal(self):
        if self.term_id_to_toggle != -1: #addition or deletion of a term was not accepted
            self.toggle_term(self.term_id_to_toggle)
            if self.LEARN_PENALIZATION_PARAMETER:
                self.penalization_parameter=self.old_penalization_parameter
                self.nr_penalization_parameter=self.old_nr_penalization_parameter            
        elif self.term_id_to_remove != -1: #toggling of a term was not accepted
            self.add_term(self.term_id_to_remove)
            self.remove_term(self.term_id_to_add)            
        elif self.parameter_changed =='alpha': #change of alpha was not accepted
            self.nr_alpha = self.old_nr_alpha
            self.alpha = self.possible_values_for_alpha_and_beta[self.nr_alpha]
            self.log_alphas = self.possible_log_alphas[self.nr_belief][self.nr_alpha]
            self.NP_log_alpha = self.old_NP_log_alpha            
        elif self.parameter_changed =='beta': #change of beta was not accepted
            self.nr_beta = self.old_nr_beta
            self.beta = self.possible_values_for_alpha_and_beta[self.nr_beta]
            self.log_1_m_betas=self.possible_log_1_minus_betas[self.nr_belief][self.nr_beta]
            self.EP_log_1_minus_beta=self.old_EP_log_1_minus_beta
        else:
            self.nr_belief = self.old_nr_belief
            self.belief = self.possible_values_for_belief[self.nr_belief]
            self.log_alphas = self.possible_log_alphas[self.nr_belief][self.nr_alpha]
            self.log_1_m_betas=self.possible_log_1_minus_betas[self.nr_belief][self.nr_beta]
            self.NP_log_alpha = self.old_NP_log_alpha    
            self.EP_log_1_minus_beta=self.old_EP_log_1_minus_beta        

    def get_score(self):
        '''maximize [sum_i EP_i * (1-beta_i)] + [sum_i NP_i * alpha_i] + EU * beta + NU *(1-alpha)
        with respect to self.selected terms, alpha, beta and prob'''
        L = self.EP_log_1_minus_beta + self.NP_log_alpha + math.log(self.beta)*self.EU + math.log(1-self.alpha) * self.NU
        if self.penalization_parameter!=1:
            L += math.log(self.penalization_parameter) * self.number_of_selected_terms + math.log(1.0-self.penalization_parameter) * self.number_of_unselected_terms
        return L

    def initialize(self,initial_set_of_terms=[]):
        self.n_unperturbed_genes = self.n_genes-self.n_perturbed_genes
        
        self.n_annotations_per_perturbed_gene = np.zeros(self.n_perturbed_genes,dtype=int)#[0]*s
        self.n_selected_terms_gene_is_annotated_to=[0]*self.n_genes
        self.TERM_IS_CURRENTLY_SELECTED=[False]*self.n_terms
        self.number_of_unselected_terms=self.n_terms
        self.number_of_selected_terms=0
        
        if self.proportion_parameter_change>0:
            self.alpha_beta_max_upperbound = (1+self.max_belief) / (4*self.max_belief)
            if self.alpha_beta_max<=self.alpha_beta_max_upperbound:
                self.possible_values_for_alpha_and_beta = list(np.linspace(0,self.alpha_beta_max,self.A+1)[1:])
            else:
                #self.possible_values_for_alpha_and_beta = list(np.linspace(0,self.alpha_beta_max_upperbound,self.A+2)[1:-1]), old line, changed 2022/8/30
                print('\n!!! Warning:\nThe chosen upper bound (%f) for alpha and beta, the FPR and FNR, is higher than the maximal possible upper bound (%f) so that genes annotated to activated terms are never penalized more than those not annotated to activated terms (i.e. need to ensure 1-beta_k > alpha). Upper bound set to %f\n' % (self.alpha_beta_max,self.alpha_beta_max_upperbound,self.alpha_beta_max_upperbound))
                self.alpha_beta_max = self.alpha_beta_max_upperbound
                self.possible_values_for_alpha_and_beta = list(np.linspace(0,self.alpha_beta_max,self.A+1)[1:])
            #every time alpha or beta are changed, comparison to the following number ensures that no perturbed annotated gene is penalized more than an unperturbed, and vice versa for unperturbed not-annotated genes
            #self.max_nr_alpha_plus_nr_beta = int(self.alpha_beta_max_upperbound//self.possible_values_for_alpha_and_beta[0]-2) 

            self.possible_values_for_belief = list(np.linspace(1,self.max_belief,self.B))
            
            self.possible_log_alphas = []
            self.possible_log_1_minus_betas = []
            for i in range(self.B):
                possible_alphas_and_betas=np.array([self.get_specific_rates(self.possible_values_for_alpha_and_beta[i],self.n_perturbed_genes,self.belief) for i in range(self.A)])
                self.possible_log_alphas.append( [np.array(el) for el in np.log(possible_alphas_and_betas)])
                self.possible_log_1_minus_betas.append( [np.array(el) for el in np.log(1-possible_alphas_and_betas)] )
                
            (self.alpha,self.nr_alpha)=self.discretize(self.possible_values_for_alpha_and_beta,self.alpha)
            (self.beta,self.nr_beta)=self.discretize(self.possible_values_for_alpha_and_beta,self.beta)
            (self.belief,self.nr_belief)=self.discretize(self.possible_values_for_belief,self.initial_belief)
            
            self.log_alphas=self.possible_log_alphas[self.nr_belief][self.nr_alpha]
            self.log_1_m_betas=self.possible_log_1_minus_betas[self.nr_belief][self.nr_beta]
            
        else:
            alphas=self.get_specific_rates(self.alpha,self.n_perturbed_genes,self.belief_alpha) #xxxxx, add number of rates
            betas=self.get_specific_rates(self.beta,self.n_perturbed_genes,self.belief_beta) #xxxxx
            self.log_alphas = [math.log(el) for el in alphas]
            self.log_1_m_betas = [math.log(1-el) for el in betas]

        if self.LEARN_PENALIZATION_PARAMETER:
            self.penalization_parameter=min(max(len(initial_set_of_terms),1),self.P)*1./self.n_terms
            self.nr_penalization_parameter=min(max(len(initial_set_of_terms),1),self.P)-1
        
        #this list keeps track of which terms are currently selected (the first self.number_of_selected_terms)
        self.internal_list_of_terms_for_MCMC=list(range(self.n_terms))
        self.position_of_term_in_internal_list=list(range(self.n_terms))

        #Initially, no term is selected
        self.NP_log_alpha = sum(self.log_alphas)
        self.NU = self.n_unperturbed_genes
        self.EU,self.EP_log_1_minus_beta = 0,0

        for term in initial_set_of_terms:
            self.add_term(term)

        self.score=self.get_score()
        self.max_score= self.score
        self.max_score_step = -1
        self.max_score_alpha = self.alpha
        self.max_score_beta = self.beta
        self.max_score_prob = self.penalization_parameter
        
        self.best_set = initial_set_of_terms

        self.n_term_selected_in_MCMC = [0]*self.n_terms
        		        		
    def MCMC_alg(self,initial_set_of_terms=[],verbose=0):
        MCMC_start=time.time()
        
        self.initialize(initial_set_of_terms)          

        learned_params=[]
        
        self.step=1
        total_steps=self.burnin_steps + self.MCMC_steps
        alpha_count = [0]*self.A
        beta_count = [0]*self.A
        penalization_parameter_count = [0]*self.P
        belief_count = [0]*self.B

        neighborhood_size = self.n_terms + self.number_of_unselected_terms * self.number_of_selected_terms
        
        PRINT_AFTER_X_STEPS=(50000 if total_steps>500000 else total_steps/5)
        
        while self.step < total_steps:              
            self.propose_MCMC_state(neighborhood_size)
            new_score = self.get_score()
            if self.term_id_to_toggle!=-1:
                new_neighborhood_size = self.n_terms + self.number_of_unselected_terms * self.number_of_selected_terms
                log_accept_probability = (new_score - self.score) + math.log(neighborhood_size*1./new_neighborhood_size)
            else:
                log_accept_probability = new_score - self.score

            if math.log(random.random()) >= log_accept_probability:
                self.undo_proposal()
            else:
                self.score = new_score
                if self.term_id_to_toggle!=-1:
                    neighborhood_size = new_neighborhood_size
            
            #check if a new local optimum has been found
            if self.score > self.max_score:
                #Remember the score and the configuration
                self.max_score = self.score
                self.max_score_step = self.step
                self.max_score_alpha = self.alpha
                self.max_score_beta = self.beta
                self.max_score_prob = self.penalization_parameter
                self.best_set = self.internal_list_of_terms_for_MCMC[self.number_of_unselected_terms:self.n_terms]
                          
            if self.step>=self.burnin_steps:
                for i in range(self.number_of_selected_terms):
                    self.n_term_selected_in_MCMC[self.internal_list_of_terms_for_MCMC[i]] += 1
                if self.proportion_parameter_change>0:
                    alpha_count[self.nr_alpha] += 1
                    beta_count[self.nr_beta] += 1
                    if self.LEARN_PENALIZATION_PARAMETER:
                        penalization_parameter_count[self.nr_penalization_parameter-1] += 1   
                    belief_count[self.nr_belief] += 1 

            if self.GET_INFO_ON_CURRENT_SETS:
                if self.step==self.burnin_steps:                    
                    self.current_set_at_burnin = self.internal_list_of_terms_for_MCMC[:self.number_of_selected_terms]
    
                if self.step == total_steps - 1:                    
                    self.current_set_at_end = self.internal_list_of_terms_for_MCMC[:self.number_of_selected_terms]
                        
            self.step += 1
                
            if self.step%PRINT_AFTER_X_STEPS==0:
                if verbose>0:
                    print("%i out of %i steps completed; passed MCMC time %f" % (self.step,total_steps,time.time()-MCMC_start))
                    if verbose>1:
                        print(self.step, self.score, self.alpha,self.beta, self.belief, str(self.nr_penalization_parameter)+'/'+str(self.n_terms) if self.LEARN_PENALIZATION_PARAMETER else 0, self.EP_log_1_minus_beta,self.NP_log_alpha,self.EU,self.NU)
                learned_params.append([self.step,self.score,self.alpha,self.beta,self.belief,str(int(round(self.penalization_parameter*len(self.T))))+"/"+str(self.n_terms),self.penalization_parameter,self.number_of_selected_terms])
            
        #Get posterior probability distribution
        Final_Distribution = np.array(self.n_term_selected_in_MCMC)*1./self.MCMC_steps #Final_Distribution = x*1./self.MCMC_steps for x in self.n_term_selected_in_MCMC]  
        
        # Useful to return a best term set - can base this off of an arbitrary cutoff of appearances in the chain
        # Order the terms from highest frequency to lowest
        
        if self.proportion_parameter_change>0:            
            hist_alpha=([a*1./self.MCMC_steps for a in alpha_count],self.possible_values_for_alpha_and_beta)
            hist_beta=([b*1./self.MCMC_steps for b in beta_count],self.possible_values_for_alpha_and_beta)
            hist_prob=([p*1./self.MCMC_steps for p in penalization_parameter_count],range(1,self.P+1))
            hist_belief=([b*1./self.MCMC_steps for b in belief_count],self.possible_values_for_belief)
                              
        if self.proportion_parameter_change>0:
            return (self.max_score, Final_Distribution, learned_params, hist_alpha, hist_beta, hist_prob, hist_belief)
        return (self.max_score, Final_Distribution, learned_params,[],[],[],[])
 
    ########################################
    ###  Data analysis/plotting methods  ###
    ########################################

    def mean_gene_level(self,list_of_terms):
        """For a list of GO terms, C, this method returns the average prediction level of all genes annotated by a term"""
        output=[]
        for term in list_of_terms:
            helper=[self.levels[self.T[term][i]] for i in range(len(self.T[term]))]
            res=np.mean(helper) if helper!=[] else 0
            output.append(res)
        return output

    def mean_position_perturbed_genes(self,list_of_terms):
        """For a list of GO terms, C, this method returns the average position (in % of all perturbed genes) of all perturbed genes annotated by a term"""
        output = []
        for term in list_of_terms:
            if self.T[term][0]>=self.n_perturbed_genes: #if no perturbed gene is annotated to the term
                output.append(-1)
            else:
                index = sum(np.array(self.T[term]) < self.n_perturbed_genes)
                output.append(np.mean(self.T[term][:index]))
        return output

    def jaccard(self, liste,verbose=False):
        """Returns a table, in which all terms in liste are pairwise compared w.r.t their Jaccard index"""
        res=[[0 for i in range(len(liste))] for j in range(len(liste))]
        for i in range(len(liste)):
            for j in range(i+1,len(liste)):
                if i==j:
                    continue
                if self.T[liste[i]]!=[] or self.T[liste[j]]!=[]:
                    res[i][j]=len(set(self.T[liste[i]])&set(self.T[liste[j]]))*1./len(set(self.T[liste[i]])|set(self.T[liste[j]]))
                    res[j][i]=res[i][j]
                if verbose>0:
                    print(self.term_names[liste[i]],'&',self.term_names[liste[j]],res[i][j])
        if len(liste)==2:
            return res[0][1]
        return res

    #########################
    ###  P-value methods  ###
    #########################
    
    @staticmethod
    def hypergeometric_pmf(x, m, n, k):
        a = CRFE.log_binomial_coefficient(m, x)
        b = CRFE.log_binomial_coefficient(n, k-x)
        c = CRFE.log_binomial_coefficient(m+n, k)
        return np.exp(a+b-c)

    @staticmethod
    def log_binomial_coefficient(population, sample):
        if sample==0:
            return 0
        s = max(sample, population - sample)
        assert s <= population
        assert population > -1
        if s == population:
            return 0
        log_factor = 0
        for i in range(s+1, population + 1):
            log_factor += math.log(i)-math.log(i-s)
            return log_factor
    
    def p_value(self,n_annotated_perturbed_genes,n_annotated_genes,n_genes,n_perturbed_genes):
        s=0
        for i in range(n_annotated_perturbed_genes,min(n_perturbed_genes,n_annotated_genes)+1):
            s+=CRFE.hypergeometric_pmf(i,n_annotated_genes,n_genes-n_annotated_genes,n_perturbed_genes)
        return s
    
    def p_values_of_list(self, list_of_terms):
        ps = []
        
        perturbed_genes = set(range(self.n_perturbed_genes))
        
        for term in list_of_terms:
            annotated_genes = set(self.T[term])
            n_annotated_genes = len(annotated_genes)
            n_annotated_perturbed_genes = len(annotated_genes & perturbed_genes)
            ps.append(self.p_value(n_annotated_perturbed_genes,n_annotated_genes,self.n_genes,self.n_perturbed_genes))
        return ps
   
    ##################################################
    ###  Methods for initial clean-up of ontology  ###
    ##################################################

    def combine_equal_terms(self,verbose=0):
        #print "merge equal GO terms into one..."
        nr_terms = len(self.term_names)
        
        ind=self.n_genes_per_term.argsort()
        cc_ord=self.n_genes_per_term[ind]
        keep_these = np.zeros(0,dtype=int)
        delete_these = np.zeros(0,dtype=int)
        for i in range(nr_terms):
            for j in range(i+1,nr_terms):
                if cc_ord[i]<cc_ord[j]:
                    break
                elif set(self.T[ind[i]])==set(self.T[ind[j]]):
                    keep_these = np.append(keep_these,ind[i])
                    delete_these = np.append(delete_these,ind[j])

        term_should_be_kept = np.ones(nr_terms,dtype=bool)
        term_should_be_kept[delete_these] = False
        
        for i in range(keep_these.shape[0]):
            self.term_names[keep_these[i]]+=' & '+self.term_names[delete_these[i]]
            
        self.T = self.T[term_should_be_kept]
        self.term_names = self.term_names[term_should_be_kept]
        self.n_genes_per_term = self.n_genes_per_term[term_should_be_kept]
        
        if verbose>0:
            print("New combined nodes:\n")
            for name in self.term_names:
                if name.find('&')!=-1:
                    print(name)
    
    def delete_too_large_and_too_small_terms(self):
        if self.upper_cutoff>=self.lower_cutoff:
            indices_to_keep = np.bitwise_and(self.n_genes_per_term>=self.lower_cutoff,self.n_genes_per_term<=self.upper_cutoff)
        else:
            indices_to_keep = self.n_genes_per_term>=self.lower_cutoff
        self.T = self.T[indices_to_keep]
        self.term_names = self.term_names[indices_to_keep]
        self.n_genes_per_term = self.n_genes_per_term[indices_to_keep]

    def remove_genes_not_annotated_to_any_term(self,verbose=0):
        toberemoved=list(set(range(self.n_genes))-set.union(*[set([])] + [set(t) for t in self.T]))
        if toberemoved==[]:
            return
        toberemoved = np.array(toberemoved,dtype=int)
        delete_these_genes = np.zeros(self.n_genes,dtype=bool)
        delete_these_genes[toberemoved] = True
        if verbose>0:
            print("Deleted genes: "+', '.join(self.unique_genes[delete_these_genes]))
        cumsum = np.cumsum(delete_these_genes)
        for j in range(len(self.T)):
            for i in range(len(self.T[j])):
                self.T[j][i]-=cumsum[self.T[j][i]]

        self.unique_genes = self.unique_genes[~delete_these_genes]
        self.n_genes -= len(toberemoved)
        
    ############################
    ###  Load and Save Data  ###
    ############################
    
    def savevar(self, variable, v):
        """
            to save one of three important variables: T, unqiue_genes, term_names
        """
        if not os.path.exists('saved_data/'):
            os.makedirs('saved_data/')
        f=open('saved_data/save_'+self.out_str+'_'+str(self.lower_cutoff)+'to'+str(self.upper_cutoff)+'annotations_'+variable+'_py'+str(which_python)+'.txt','wb')
        pickle.dump(v, f)
        f.close()
    
    def loadvar(self, variable):
        """
            to load one of three important variables: T, unqiue_genes, term_names
        """
        if not os.path.exists('saved_data/'):
            os.makedirs('saved_data/')
        try:
            f=open('saved_data/save_'+self.out_str+'_'+str(self.lower_cutoff)+'to'+str(self.upper_cutoff)+'annotations_'+variable+'_py'+str(which_python)+'.txt','rb')
            print('saved_data/save_'+self.out_str+'_'+str(self.lower_cutoff)+'to'+str(self.upper_cutoff)+'annotations_'+variable+'_py'+str(which_python)+'.txt',f)
        except IOError:
            return None
        try:
            return pickle.load(f)
        except ValueError: #If the file exists but was written by Python 3 in pickle protocol 3 and is now loaded in Python 2.x, just recreate it and overwrite it
            return None
    
    def save_all(self, T, term_names, unique_genes):
        self.savevar('T', T)
        self.savevar('term_names', term_names)
        self.savevar('unique_genes', unique_genes)
    
    def load_annotation_data(self,verbose=0):
        """
            load all important variables from previously generated files, loads
            variables for full version (i.e. for lower_cutoff=1, upper_cutoff=0) and 
            then deletes terms that are too small or too large.
        """
        
        lower_cutoff=self.lower_cutoff
        upper_cutoff=self.upper_cutoff
        self.lower_cutoff=1
        self.upper_cutoff=0
        
        self.unique_genes=self.loadvar('unique_genes')
        self.T=self.loadvar('T')
        self.term_names=self.loadvar('term_names')
        
        if self.unique_genes is None or self.T is None or self.term_names is None:
            self.create_annotation_data(verbose)
        self.n_genes = len(self.unique_genes)
        
        self.lower_cutoff=lower_cutoff
        self.upper_cutoff=upper_cutoff

        self.n_genes_per_term = np.array(list(map(len,self.T)))        
        
        #remove terms that have too many or too little annotations
        self.delete_too_large_and_too_small_terms()
        self.n_terms = len(self.T)
        if self.P>0.5*self.n_terms:
            self.P = int(np.floor(0.5*self.n_terms))
        #remove genes that are no longer annotated by any term
        self.remove_genes_not_annotated_to_any_term(verbose)
        
        #sort genes in self.T and self.unique_genes to align with the ranking, e.g. gene 0 is the top ranked gene
        pass
    
    def create_annotation_data(self,verbose):
        """creates all the important variables from the annotations_file."""

        data = pd.read_csv(self.annotations_file,sep='\t',header=None,names=['term','genes'])
        self.term_names = np.array(data['term'])
        self.unique_genes = np.array(list(set.union(*[set([])] + [set(genes.split(' ')) for genes in data['genes']])))
        dict_genes = dict(zip(self.unique_genes,range(len(self.unique_genes))))
        dict_genes_list = dict(zip(self.genes_list,range(len(self.genes_list))))
        
        #create the annotations and delete genes in self.unique_genes but not in self.genes_list
        T = []
        for genes in data['genes']:
            T.append([])
            for gene in genes.split(' '):
                try: #check if gene in annotation data is also in expression data
                    dict_genes_list[gene]
                    T[-1].append(dict_genes[gene])
                except KeyError:
                    pass
            #T[-1]= np.sort(T[-1])
        self.T = np.array(T,dtype=object)	
        #self.T = np.array([np.sort(list(map(lambda x: dict_genes[x],genes.split(' ')))) for genes in data['genes']])
        
        self.n_genes_per_term = np.array(list(map(len,self.T)))
        self.n_genes = len(self.unique_genes)
        
        
        #remove terms that have too many or too little annotations
        self.delete_too_large_and_too_small_terms()
                
        #remove genes that are no longer annotated by any term
        self.remove_genes_not_annotated_to_any_term(verbose)
        
        self.savevar('T', self.T)
        self.savevar('term_names', self.term_names)
        self.savevar('unique_genes', self.unique_genes)
    
    ########################
    ###  Output Methods  ###
    ########################
    
    def get_filename(self,method,LONG=False):
        filename='%s_random%s_param_learning%s_%sto%sannotations' % (self.out_str,int(1000*method[1]),int(self.proportion_parameter_change*1000),self.lower_cutoff,self.upper_cutoff)
        if LONG:
            filename+='_belief%s_nr_perturbed%i' % (int(self.belief*1000),int(self.n_perturbed_genes))
        now = datetime.datetime.now()
        filename+=now.strftime("_%Y-%m-%d-%H-%M-%S")
        return filename
        
    def writeOut_MCMC(self, filename, learned_params, hist_alpha=[], hist_beta=[], hist_belief=[], hist_prob=[]):
        nr_terms=len(self.T)
        filename=self.output_folder+'mcmc_'+filename
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        f_out = open(filename+'.html', 'w')
        names=['Iteration step','Likelihood value','alpha (FPR)', 'beta (FNR)', 'belief','q (rational)','q (decimal)', '|C|']
        strings="<!DOCTYPE html>\n<html><head><title>Info about MCMC Process</title><link rel='stylesheet' href='../abc.css' type='text/css' media='screen'></head>\n<body>\n<h1>Additional info about the MCMC process</h1>"
        strings+="<h2>Learned parameter values</h2><table>"
        strings+="<tr valign=\"top\"><td>False positive rate (alpha):</td><td>%s</td></tr>" % (self.alpha)
        strings+="<tr valign=\"top\"><td>False negative rate (beta):</td><td>%s</td></tr>" % (self.beta)
        strings+="<tr valign=\"top\"><td>Belief parameter (beta):</td><td>%s</td></tr>" % (self.belief)
        strings+="<tr valign=\"top\"><td>Penalization parameter (q):</td><td>%s</td></tr>" % (self.penalization_parameter)
        strings+="<tr valign=\"top\"><td>Penalization per term:</td><td>%s</td></tr>" % (-math.log(self.penalization_parameter/(1-self.penalization_parameter)))
        strings+="</table>"
        
        if hist_alpha!=[]:
            strings+="<h2>Posterior Parameter Value Distribution</h2><table style=\"table-layout:fixed;\" cellspacing=\"0\" cellpadding=\"2\">"
            strings+="<tr style=\"font-weight:bold;\" valign=\"top\"><td>Value</td>"
            for i in range(len(hist_alpha[1])):
                strings+="<td style='width:50px;'>%s</td>" % (round(hist_alpha[1][i],3))
            strings+="</tr><tr style=\"background-color:#EEEEEE\" valign=\"top\"><td>alpha</td>"
            for i in range(len(hist_alpha[1])):
                strings+="<td style='width:50px;'>%s</td>" % (round(hist_alpha[0][i],3) if 1>hist_alpha[0][i]>0 else int(hist_alpha[0][i]))
            strings+="</tr><tr valign=\"top\"><td>beta</td>"
            for i in range(len(hist_beta[1])):
                strings+="<td style='width:50px;'>%s</td>" % (round(hist_beta[0][i],3) if 1>hist_beta[0][i]>0 else int(hist_beta[0][i]))
            strings+="</tr></table><p></p><table style=\"table-layout:fixed;\" cellspacing=\"0\">"
            strings+="<tr style=\"font-weight:bold;\" valign=\"top\"><td>q*%i</td>" % nr_terms
            for i in range(len(hist_belief[1])):
                strings+="<td>%s</td>" % (hist_belief[1][i])
            strings+="</tr><tr style=\"background-color:#EEEEEE\" valign=\"top\"><td>distr </td>"
            for i in range(len(hist_belief[1])):
                strings+="<td style='width:50px;'>%s</td>" % (round(hist_belief[0][i],3) if 1>hist_belief[0][i]>0 else int(hist_belief[0][i]))
            strings+="</tr></table><p></p><table style=\"table-layout:fixed;\" cellspacing=\"0\">"
            strings+="<tr style=\"font-weight:bold;\" valign=\"top\"><td>q*%i</td>" % nr_terms
            for i in range(len(hist_prob[1])):
                strings+="<td>%s</td>" % (hist_prob[1][i])
            strings+="</tr><tr style=\"background-color:#EEEEEE\" valign=\"top\"><td>distr </td>"
            for i in range(len(hist_prob[1])):
                strings+="<td style='width:50px;'>%s</td>" % (round(hist_prob[0][i],3) if 1>hist_prob[0][i]>0 else int(hist_prob[0][i]))
            strings+="</tr></table>"       
        
        strings+="<h2>Time course</h2><table cellspacing=\"0\"><tr valign=\"top\"><tr>"
        for i in range(len(names)):
            strings+="<td>%s</td>" % names[i]
        strings+="</tr>\n"
        for j in range(len(learned_params)):
            strings+="<tr style=\"background-color:%s;\">" % ("#EEEEEE" if j%2==0 else "#FFFFFF")
            for i in range(len(names)):
                strings+="<td style=\"width:150px;\">%s</td>" % (str(round(learned_params[j][i],3) if i in [1,2,3] else (round(learned_params[j][i],5) if i==6 else str(learned_params[j][i]))))
            strings+="</tr>\n"
        f_out.write(strings)
        f_out.write("</table>\n</body>\n</html>")
        f_out.flush()
        f_out.close() 
        
    def writeOut(self, filename, IS, mean_max_score, std_max_score, enriched, ps, mean_gene_level, mean_position_perturbed_genes, method, avg_posteriors, std_posteriors):
        filename=self.output_folder+'result_'+filename
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        text=["Creation Time:\t" + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n"]
        text.append("INPUT:")
        
        str1 = "File name:\t%s\nAnnotations file:\t%s\nGene file:\t%s\nlower limit annotations/process:\t%s\nupper limit annotations/process:\t%s\nthreshold:\t%s\nthreshold type:\t%s\nbelief:\t%s\nMethod:\t%s\nrepeats:\t%s\nburnin steps:\t%s\nMCMC steps:\t%s\nMaximal value for alpha and beta:\t%s\nMaximal possible value for alpha and beta:\t%s\nProbability parameter change:\t%s\nseed for RNG:\t%s\n" % (sys.argv[0],self.annotations_file, 
        self.gene_file, self.lower_cutoff, self.upper_cutoff if self.upper_cutoff>=self.lower_cutoff else "0 (i.e., None)",
        self.threshold,self.threshold_type,
        self.belief, method[0], self.repeats, self.burnin_steps, self.MCMC_steps, 
        self.alpha_beta_max if self.proportion_parameter_change>0 else 'fixed', self.alpha_beta_max_upperbound if self.proportion_parameter_change>0 else 'fixed',
        self.proportion_parameter_change, self.seed)
    
        text.append(str1)
        text.append("Initial set of processes:")
        if IS==[]:
            text.append('Empty set')
        else:
            for i in IS:
                text.append(str(i) + "\t" + self.term_names[i])
        text.append("\nINPUT ANALYSIS:")
        text.append("Total number of processes:\t"+ str(self.n_terms)) 
        text.append("Total number of genes:\t"+ str(self.n_genes))
        text.append("Number of perturbed genes:\t" +str(self.n_perturbed_genes))
        text.append("Proportion perturbed genes:\t" +str(self.n_perturbed_genes/self.n_genes))

        text.append("\nOUTPUT:")
        text.append("Learnt parameters:")
        str2 = "alpha (FPR):\t%s\nbeta (FNR):\t%s\nq (penalization parameter):\t%s\n" % ( self.alpha if self.proportion_parameter_change>0 else str(self.alpha) + ' (fixed)', self.beta if self.proportion_parameter_change>0 else str(self.beta) + ' (fixed)', self.penalization_parameter if self.LEARN_PENALIZATION_PARAMETER else str(self.penalization_parameter) + ' (fixed)')
        text.append(str2+"\n")
        text.append("Best log-likelihood value (mean +- standard deviation across replicates):\t"+str(mean_max_score)+' +- '+str(std_max_score)+'\n')
        text.append("Order\t#(Annotations)\t#(Annotations of perturbed genes)\tMean Posterior Probability\tStandard Deviation Posterior Probability\t"+("" if ps==[] else "p-value\t")+"Average Expression Level\tAverage position perturbed genes\tTerm Name" )
        if self.show_genes_in_output: 
            text.append("\tExplained Perturbed Genes") 
              
        set_perturbed_genes = set(self.genes_list[:self.n_perturbed_genes])
        out_list = [str(count+1) + "\t" + str(i[2]) + "\t" + str(len(set(self.T[i[0]])&set_perturbed_genes)) for (count,i) in zip(range(len(enriched)),enriched)]
        
        for i in range(len(out_list)):
            out_list[i]+="\t"+str(round(avg_posteriors[i],5))+"\t"+str(round(std_posteriors[i],5))

        if self.show_genes_in_output:
            ind=sorted(range(len(self.levels)), reverse=True, key=lambda k: self.levels[k])
            ind=ind[:self.n_perturbed_genes]                                
        for i in range(len(mean_gene_level)):
            if ps!=[]:
                out_list[i] += "\t" + str("%.3e" % ps[i])
            out_list[i]+="\t" + str(round(mean_gene_level[i],3)) + "\t" + str(round(mean_position_perturbed_genes[i],3)) + "\t" + str(enriched[i][1])
            if self.show_genes_in_output:
                out_list[i]+='\t'
                for j in ind:
                    if j in self.T[self.terms[i]]:
                        out_list[i]+=self.unique_genes[j]+" "
                out_list[i]=out_list[i][:-1]
            text.append(out_list[i])

        f_out = open(filename+".txt", 'w+')
        f_out.write("\n".join(text))
        f_out.close()

    def determine_method(self, parameter_initial_MCMC_set):
        if parameter_initial_MCMC_set>=1:
            initial_set_of_terms_description="a randomly chosen set of "+str(min(parameter_initial_MCMC_set,self.n_terms))+" processes"
        elif parameter_initial_MCMC_set>0:
            initial_set_of_terms_description="a randomly chosen set of processes (each process is in the set with probability p="+ str(parameter_initial_MCMC_set)+")"    
        else:       
            initial_set_of_terms_description="an empty set"
        text="MCMC"
        if self.proportion_parameter_change>0:
            text+=" with parameter learning"
        text+=" with "+initial_set_of_terms_description+" as initial set"
        return [text,parameter_initial_MCMC_set]

    #####################################
    ###  Main Function of this Class  ###
    #####################################        
    
    def runMe(self, verbose=False, parameter_initial_MCMC_set=0):        
        start=time.time()
        
        #load the expression data, initialize ordered list of genes (self.genes_list) and their expression levels (self.level_list)
        self.get_expression_data()
        
        if verbose:
            print('Checkpoint A, passed time',time.time()-start)
        
        #load/create the annotation data, initialize self.term_names, self.unique_genes, self.T, self.n_genes_per_term
        self.load_annotation_data()
        
        if verbose:
            print('Checkpoint B, passed time',time.time()-start)
        
        #If a gene is in expression data but not in annotation data, it is ignored.
        
        #up to v19, genes in self.unique_genes were sorted alphabetically, since v20 we sort by ranking in expression data
        
        pos_unique_genes = dict(zip(self.unique_genes,range(self.n_genes)))
        dict_new_old_index_unique_genes = dict()
        count=0
        self.levels = []
        self.unique_genes = []
        for gene,level in zip(self.genes_list,self.level_list): #self.genes_list is sorted already by expression value, top to bottom
            try:
                pos = pos_unique_genes[gene]
                dict_new_old_index_unique_genes.update({pos:count})
                count+=1
                self.levels.append(level)
                self.unique_genes.append(gene)
            except KeyError:
                continue
            
        for i in range(self.n_terms):
            self.T[i] = [dict_new_old_index_unique_genes[gene] for gene in self.T[i]]
        
        if self.threshold_type=='proportion':
            self.n_perturbed_genes = int(np.round(self.n_genes * self.threshold))
        else:
            self.n_perturbed_genes = sum(np.array(self.levels) >= self.threshold)
                           
        #set random seeds for both random number generators (both are used for performance boosts)
        if self.seed==-1:
            self.seed=np.random.randint(2**32 - 1)
        np.random.seed(self.seed)
        random.seed(self.seed)
        
        #choose random set of terms to start with
        initial_set_of_terms = []
        if parameter_initial_MCMC_set>0:
            if parameter_initial_MCMC_set>=1:
                initial_set_of_terms = np.random.choice(self.n_terms,int(min(self.n_terms,parameter_initial_MCMC_set)),replace=False)
            else:
                initial_set_of_terms = np.where(np.random.random(self.n_terms)<parameter_initial_MCMC_set)
        
        if verbose:
            print('Checkpoint C, passed time',time.time()-start)
        
        #Start of the actual algorithm
        if self.proportion_parameter_change>0:
            hist_alphas,hist_betas,hist_probs,hist_beliefs=[],[],[],[]
        max_scores,posteriors=[],[]
        if self.GET_INFO_ON_CURRENT_SETS:
            current_sets_at_burnin = []
            current_sets_at_end = []
        for i in range(self.repeats):
            (max_score, MCMC_Distr, learned_params, hist_alpha, hist_beta, hist_prob, hist_belief)=self.MCMC_alg(initial_set_of_terms=initial_set_of_terms,verbose = 2 if verbose else -1)
            max_scores.append(max_score)
            posteriors.append(MCMC_Distr)
            if self.proportion_parameter_change>0:
                hist_alphas.append(hist_alpha[0])
                hist_betas.append(hist_beta[0]) 
                hist_probs.append(hist_prob[0])
                hist_beliefs.append(hist_belief[0])
            if self.GET_INFO_ON_CURRENT_SETS:
                current_sets_at_burnin.append(self.current_set_at_burnin)
                current_sets_at_end.append(self.current_set_at_end)

        if self.proportion_parameter_change>0:
            mean_hist_alphas = np.mean(np.array(hist_alphas),0)
            mean_hist_betas = np.mean(np.array(hist_betas),0)
            mean_hist_probs = np.mean(np.array(hist_probs),0)
            mean_hist_beliefs = np.mean(np.array(hist_beliefs),0)
            self.alpha=self.possible_values_for_alpha_and_beta[np.argmax(mean_hist_alphas)]
            self.beta=self.possible_values_for_alpha_and_beta[np.argmax(mean_hist_betas)]
            self.belief=self.possible_values_for_belief[np.argmax(mean_hist_beliefs)]
            if self.LEARN_PENALIZATION_PARAMETER:
                self.penalization_parameter=(np.argmax(mean_hist_probs)+1)*1./self.n_terms
            
        avg_posteriors = np.mean(np.array(posteriors),0)
        std_posteriors = np.std(np.array(posteriors),0)                    
        self.terms=np.array(sorted(range(self.n_terms), reverse=True, key=lambda k: avg_posteriors[k]))
        avg_posteriors=avg_posteriors[self.terms]
        std_posteriors=std_posteriors[self.terms]
        
        #self.new_output(avg_posteriors,std_posteriors)
        
        enriched = []
        for i in range(len(self.terms)):
            #if avg_posteriors[i]<1./10 and i>20 or i>100: #output between 20 and 100 terms, stop when the posterior falls below 0.1
            if avg_posteriors[i]==0: #output between 20 and 100 terms, stop when the posterior falls below 0.1
                break
            enriched.append((self.terms[i], self.term_names[self.terms[i]], len(self.T[self.terms[i]])))
        
        if verbose:
            print('Checkpoint D, passed time',time.time()-start)
            
        ps = self.p_values_of_list(self.terms[:len(enriched)])        
        mean_gene_level=self.mean_gene_level(self.terms[:len(enriched)])
        mean_position_perturbed_genes=self.mean_position_perturbed_genes(self.terms[:len(enriched)])

        method=self.determine_method(parameter_initial_MCMC_set)
        filename=self.get_filename(method,LONG=True)                
                                                
        self.writeOut(filename, initial_set_of_terms, np.mean(max_scores), np.std(max_scores), enriched, ps, mean_gene_level, mean_position_perturbed_genes, method, avg_posteriors, std_posteriors)
        if self.proportion_parameter_change>0:        
            self.writeOut_MCMC(filename, learned_params, [mean_hist_alphas,self.possible_values_for_alpha_and_beta],[mean_hist_betas,self.possible_values_for_alpha_and_beta],[mean_hist_beliefs,self.possible_values_for_belief],[mean_hist_probs,np.arange(1,self.P+1)])

        if verbose:
            print('Checkpoint E, passed time',time.time()-start)    
                  
        if self.GET_INFO_ON_CURRENT_SETS:
            return (self.terms,avg_posteriors,std_posteriors,current_sets_at_burnin,current_sets_at_end)
        else:        
            return (self.terms,avg_posteriors,std_posteriors)

###############################################
###  Command Line Support via main function ###
###############################################

def main():
    '''
        This kicks off the process and parses the options from the arguments.
    '''
    p = optparse.OptionParser('This program implements the functional enrichment method CRFE, described in C Kadelka, M Brandon, TM Murali, Concise Functional Enrichment of Ranked Gene Lists. See www.math.vt.edu/people/claus89/crfe/ for help.')
            
    p.add_option('--annotations_file', '-a', default='human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv', help='Annotation data (csv or txt). Each row contains the term name followed by a tab, followed by all annotated genes separated by space. No header!')        
    p.add_option('--gene_file', '-g', default='bkz_gene_list.csv', help='Expression data (csv or txt). Each row contains the gene name, optionally followed by tab plus expression level. No header!')
    p.add_option('--threshold', '-t', default='0.3', help='If threshold_type==proportion, this describes the proportion of genes considered perturbed, otherwise it describes the threshold in the expression values to be used (default 0.3)')
    p.add_option('--threshold_type', '-T', default='proportion', help="Whether threshold should be interpreted as 1. a value, 2. proportion ('proportion') of perturbed genes, or 3. different levels of genes already provided in the expression data ('levels'), (default proportion)")
    p.add_option('--lower_cutoff', '-c', default ='20', help='Only terms with at least this many annotations are considered (default 5)')
    p.add_option('--upper_cutoff', '-C', default ='200', help='Only terms with at most this many annotations are considered. Use 0 to excluded an upper cutoff (default 200).')
    p.add_option('--belief', '-b', default='5', help='Belief for how much more active the highest gene set is compared to the least active gene set (default 5)')

    p.add_option('--repeats', '-n', default='1', help='Number of independent repeats of CRFE (default 1).')
    p.add_option('--burnin_steps', '-s', default = '100000', help='Length of (unrecorded) burnin period of MCMC simulation. Used to initialize the Markov chain (default 100000)')
    p.add_option('--MCMC_steps', '-S', default = '1000000', help='Length of (recorded) steps in the MCMC simulation after the (unrecorded) burnin period (default 1000000)')
    p.add_option('--alpha_beta_max', '-x', default = '1', help='Maximal value for the false positive and false negative rate that can be learned (default 1). Should not be changed. CRFE already automatically chooses the false poositive and negative rates such that (i) an explained perturbed gene is never penalized more than an unexplained perturbed gene, and that (ii) an unexplained unperturbed gene is never penalized more than an explained unperturbed gene.')   
    p.add_option('--diff_alpha_beta', '-X', default = '20', help='Number of different possible values for the false positive and false negative rate (default 20).')
    p.add_option('--diff_q', '-Q', default = '20', help='Number of different possible values for the term size penalization parameter q (default 20)')
    p.add_option('--probability_parameter_change', '-p', default = '0.2', help='Proportion of time parameters are changed in MCMC (default 0.2), if 0 no parameters are being learnt.')       
    p.add_option('--parameter_initial_MCMC_set','-r', default='0', help='When nonzero, the algorithm does not start with an empty set. If r>=1, r random processes are chosen, if 0<r<1, each process is in the initial set with probability r (default 0, i.e., empty set)')
   
    p.add_option('--output_folder', '-o', default='output/', help='Local reference to folder where output is saved (default: output/)')
    p.add_option('--identifier', '-i' ,default = '', help='Identifying string used to distinguish different experiments (empty default)')
    p.add_option('-v', action='store_true', dest='verbose', help='If True, more information is printed during the CRFE run (slows down CRFE).')
    p.add_option('--seed','-R', default='-1', help='If nonnegative, this seed is used too initialize the random number generators (default: -1, i.e. random seed)')
    options, arguments = p.parse_args()

    verbose = options.verbose
    
    try:
        lower_cutoff = int(options.lower_cutoff)
    except ValueError:
        lower_cutoff = 20
    try:
        upper_cutoff = int(options.upper_cutoff)
    except ValueError:
        upper_cutoff = 200

    try:
        repeats = int(options.repeats)
        if repeats<1:
            repeats=1
    except ValueError:
        repeats = 1    
            
    try:
        belief = int(options.belief)
    except ValueError:
        belief = 5
    try:
        threshold=float(options.threshold)
    except ValueError:
        threshold=0.3
    try:
        threshold_type=str(options.threshold_type)
    except ValueError:
        threshold_type='proportion'
           
    try:
        burnin_steps = int(options.burnin_steps)
    except ValueError:
        burnin_steps = 20000
    
    try:
        MCMC_steps = int(options.MCMC_steps)
    except ValueError:
        MCMC_steps = 1000000
    
    try:
        alpha_beta_max = float(options.alpha_beta_max)
    except ValueError:
        alpha_beta_max = 1
        
    try:
        proportion_parameter_change = float(options.probability_parameter_change)
        if proportion_parameter_change>1:
            proportion_parameter_change=0.2
        elif proportion_parameter_change<=0:
            proportion_parameter_change=0.01
    except ValueError:
        proportion_parameter_change = 0.2
    
    try:
        A = int(options.diff_alpha_beta)
    except ValueError:
        A = 20
        
    try:
        P = int(options.diff_q)
    except ValueError:
        P = 20

    try:
        parameter_initial_MCMC_set = float(options.parameter_initial_MCMC_set)
        if parameter_initial_MCMC_set>=1:
            parameter_initial_MCMC_set=int(parameter_initial_MCMC_set)
        if parameter_initial_MCMC_set<=0:
            parameter_initial_MCMC_set=0
    except ValueError:
        parameter_initial_MCMC_set = 0    
            
    try:
        seed = int(options.seed)
    except ValueError:
        seed = -1

    print(options.diff_alpha_beta)
    gene_file = options.gene_file
    annotations_file = options.annotations_file
    out_str = options.identifier
    output_folder = options.output_folder
    
    myGOAlgo = CRFE(repeats, lower_cutoff, upper_cutoff, belief, threshold, threshold_type, burnin_steps, MCMC_steps, alpha_beta_max, proportion_parameter_change, A, P, gene_file, annotations_file, output_folder, out_str, seed)
    myGOAlgo.runMe(verbose,parameter_initial_MCMC_set)
    return myGOAlgo
    
##Main Part

if __name__ == '__main__':
    main()

    # old data: doesn't work right now
    # lower_cutoff=20
    # upper_cutoff=200
    # gene_file = '3-cell-hep-d12_versus_2-cell-hep-d12_ID.txt'
    # annotation_file = 'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv'
    # identifier = 'LIVER'
    # repeats=1
    
    # new data: lung cancer RNA_seq
    lower_cutoff=20
    upper_cutoff=500
    belief=10
    burnin=50000
    steps=50000
    repeats=1
    threshold=0.3
    annotation_file = 'data/GOhuman_ft_named.txt'
    gene_file = 'data/GSE87340_tumor_normal_log2fc-overexp.txt'
    gene_file = 'data/GSE40419_tumor_normal_deseq2_log2fc-overexp.txt'
    identifier = 'real_hshdfsdf'
    
    m = CRFE(repeats, lower_cutoff, upper_cutoff, belief, 
                                        threshold, 'proportion', burnin, steps, 1,0,20,20, gene_file,
                                      annotation_file,'output/', identifier, seed=-1,
                                      LEARN_PENALIZATION_PARAMETER=True,penalization_parameter=0.001,GET_INFO_ON_CURRENT_SETS=True)
    (C,avgs,stds,first_sets,last_sets)=m.runMe(verbose=1, parameter_initial_MCMC_set=0)
    
    annotation_file = 'data/GOhuman_ft_named.txt'
    gene_file = 'adeno_symbol.txt'    
    identifier = 'adeno'    
    repeats = 1
    burnin=25000
    steps=25000
    m = CRFE(repeats, lower_cutoff, upper_cutoff, belief, 
                                        threshold, 'levels', burnin, steps, 1,0.2,20,20, gene_file,
                                      annotation_file,'output/', identifier, seed=-1,
                                      LEARN_PENALIZATION_PARAMETER=True,penalization_parameter=0.001,GET_INFO_ON_CURRENT_SETS=True)
    (C,avgs,stds,first_sets,last_sets)=m.runMe(verbose=1, parameter_initial_MCMC_set=0)    
    
    
X = set(range(m.n_perturbed_genes))#set(m.unique_genes[:m.n_perturbed_genes])
len_X = len(X)
U = set(range(m.n_perturbed_genes,m.n_genes))#set(m.unique_genes[m.n_perturbed_genes:])
len_U = len(U)
E = set()
EX,EU,qual = [],[],[]
for t in C[:50]:
    E = set.union(E,set(m.T[t]))
    EX.append(len(set.intersection(E,X)))
    EU.append(len(set.intersection(E,U)))
    try:
        qual.append((EX[-1]/len_X)/(EU[-1]/len_U))
    except ZeroDivisionError:
        qual.append(np.nan)
    
#        
#     within_jaccard_similarities = []
#     between_jaccard_similarities_first_sets = []
#     between_jaccard_similarities_last_sets = []
    
#     def jacc(array1,array2):
#         sarray1 = set(array1)
#         sarray2 = set(array2)
#         return len(sarray1.intersection(sarray2)) / len(sarray1.union(sarray2))
    
#     for i in range(repeats):
#         within_jaccard_similarities.append(jacc(first_sets[i],last_sets[i]))
#         for j in range(i+1,repeats):
#             between_jaccard_similarities_first_sets.append(jacc(first_sets[i],first_sets[j]))
#             between_jaccard_similarities_last_sets.append(jacc(last_sets[i],last_sets[j]))
           
#     import matplotlib.pyplot as plt
#     f,ax = plt.subplots()
#     ax.boxplot([within_jaccard_similarities],positions = [0])
#     ax.boxplot([between_jaccard_similarities_first_sets],positions = [1])
#     ax.boxplot([between_jaccard_similarities_last_sets],positions = [2])
#     ax.set_xticklabels(['within same repeat','between repeats\nat burnin','between repeats\nat end'])
#     ax.set_ylabel('Jaccard similarity of current set')




































#     lower_cutoff=20
#     upper_cutoff=500
#     belief=10
#     burnin=5000
#     steps=50000
#     repeats=1
#     threshold=0.1
#     annotation_file = 'data/GOhuman_ft_named.txt'
#     gene_file = 'data/GSE87340_tumor_normal_log2fc-overexp.txt'
#     identifier = 'real_an'
    
#     Cs = []
#     avgss = []
#     learned_alphas = []
#     learned_betas = []
#     proportions_perturbed = np.arange(0.05,0.51,0.1)
#     for threshold in proportions_perturbed:
#         m = CRFE(repeats, lower_cutoff, upper_cutoff, belief, 
#                                             threshold, 'proportion', burnin, steps, 1,0.2,20,20, gene_file,
#                                           annotation_file,'output/', identifier, seed=-1)
#         (C,avgs,stds)=m.runMe(verbose=1, parameter_initial_MCMC_set=0)
#         Cs.append(C)
#         avgss.append(avgs)
        
        
#     for i in range(len(Cs)):
#         print(sum(avgss[i]>0.5),sum(avgss[i]>0))
        
        
        
        
# import cProfile
# import pstats
# import io

# seed = 123456

# pr = cProfile.Profile()
# pr.enable()

# lower_cutoff=20
# upper_cutoff=500
# belief=10
# burnin=5000
# steps=50000
# repeats=1
# threshold=0.3
# annotation_file = 'data/GOhuman_ft_named.txt'
# gene_file = 'data/GSE87340_tumor_normal_log2fc-overexp.txt'
# identifier = 'real_an'
# m = CRFE(repeats, lower_cutoff, upper_cutoff, belief, 
#                                     threshold, 'proportion', burnin, steps, 1,0.2,20,20, gene_file,
#                                   annotation_file,'output/', identifier, seed=-1)

# my_result = m.runMe(verbose=1,parameter_initial_MCMC_set=0)

# pr.disable()
# s = io.StringIO()
# ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
# ps.print_stats()

# with open('runtime_cProfile3.txt', 'w+') as f:
#     f.write(s.getvalue())
