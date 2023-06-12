import numpy as np
import copy
import EMhard


class Bunch1:        # step_hard, step_soft, trial_hard, trial_soft
    def __init__(self, NUM_MUTATION, NUM_BLOCK, NUM_CLONE, K):
        self.mixture = np.zeros ( (NUM_BLOCK, NUM_CLONE), dtype = "float")
        self.mixture_record = np.zeros ( (K, NUM_BLOCK, NUM_CLONE), dtype = "float")
        self.membership = np.zeros ( (NUM_MUTATION), dtype = "int")
        self.membership_record = np.zeros ( (K, NUM_MUTATION), dtype = "int")
        self.membership_p = np.zeros ( (NUM_MUTATION, NUM_CLONE), dtype = "float")
        self.membership_p_record = np.zeros ( (K, NUM_MUTATION, NUM_CLONE), dtype = "float")
        self.membership_p_normalize = np.zeros ( (NUM_MUTATION, NUM_CLONE), dtype = "float")
        self.membership_p_normalize_record = np.zeros ( (K, NUM_MUTATION, NUM_CLONE), dtype = "float")
        self.likelihood = float("-inf")
        self.likelihood_record = np.array ([  float("-inf") ] * (K))
        self.stepindex = 0
        self.stepindex_record = np.zeros (K , dtype = "int")
        self.max_step_index = -1
        self.max_step_index_record = np.array ([-1] * (K))
        self.makeone_index = list (range (NUM_CLONE))
        self.makeone_index_record = [[]] * K
        self.fp_index = -1
        self.fp_index_record = np.array ([-1] * (K))
        self.fp_member_index = []
        self.fp_member_index_record = [[]] * K
        self.includefp = False
        self.includefp_record = np.array ( [False] * (K))
        self.fp_involuntary = False
        self.fp_involuntary_record = np.array ( [False] * (K))
        self.makeone_prenormalization = True
        self.makeone_prenormalization_record = np.zeros (K , dtype = "bool")

    def acc (self, mixture, membership, likelihood, membership_p, membership_p_normalize, makeone_index, fp_index, step_index, fp_member_index, includefp, fp_involuntary, makeone_prenormalization, max_step_index, K):
        self.mixture_record [K]= copy.deepcopy ( mixture )
        self.membership_record [K] = copy.deepcopy ( membership ) 
        self.likelihood_record [K] = likelihood
        self.membership_p_record [K] = copy.deepcopy ( membership_p )
        self.membership_p_normalize_record [K] = copy.deepcopy ( membership_p_normalize )
        self.makeone_index_record[K] = copy.deepcopy ( makeone_index )
        self.fp_index_record[K] = copy.deepcopy ( fp_index )
        self.stepindex = step_index
        self.stepindex_record[K] = step_index
        self.max_step_index = max_step_index
        self.max_step_index_record [K] = max_step_index
        self.fp_member_index = copy.deepcopy ( fp_member_index )
        self.fp_member_index_record [K] = copy.deepcopy ( fp_member_index )
        self.includefp = includefp
        self.includefp_record [K] = includefp
        self.fp_involuntary  = fp_involuntary 
        self.fp_involuntary_record [K] = fp_involuntary 
        self.makeone_prenormalization = makeone_prenormalization
        self.makeone_prenormalization_record [K] = makeone_prenormalization


    def find_max_likelihood_fp_voluntary (self, start, end):
        if start > end:
            return -1
        
        max, max_index = -9999999, -1
        for i in range (start, end + 1):
            if self.fp_involuntary_record[i] == False:   
                if self.makeone_prenormalization_record[i] == True:
                    if self.likelihood_record [i] > max:
                        max = self.likelihood_record [i]
                        max_index = i
        
        if max == -9999999:
            return -1
        return max_index
        
    def find_max_likelihood (self, start, end):
        try:
            i = np.argmax(self.likelihood_record [ start : end + 1 ]) + start
            return i
        except:
            return -1
        

    def copy (self, other, self_i, other_j):  # cluster_hard, cluster_soft
        self.mixture = copy.deepcopy ( other.mixture_record [ other_j ] )
        self.mixture_record [self_i] = copy.deepcopy ( other.mixture_record[other_j] )
        self.membership = copy.deepcopy ( other.membership_record [ other_j ] )
        self.membership_record [self_i] = copy.deepcopy ( other.membership_record[ other_j ] )
        self.membership_p = copy.deepcopy  ( other.membership_p_record[ other_j ] ) 
        self.membership_p_record [self_i] = copy.deepcopy ( other.membership_p_record[ other_j ] )
        self.membership_p_normalize_record [self_i] = copy.deepcopy ( other.membership_p_normalize_record[ other_j ] )
        self.likelihood = copy.deepcopy ( other.likelihood_record [ other_j ] )
        self.likelihood_record [ self_i ] = copy.deepcopy ( other.likelihood_record [other_j] )
        self.makeone_index = copy.deepcopy  ( other.makeone_index_record[ other_j ] )
        self.makeone_index_record [self_i] = copy.deepcopy ( other.makeone_index_record[ other_j ] )
        self.fp_index = other.fp_index_record[ other_j ]
        self.fp_index_record [self_i] = copy.deepcopy ( other.fp_index_record[ other_j ] )
        self.fp_member_index = other.fp_member_index_record [other_j ]
        self.fp_member_index_record [self_i] = copy.deepcopy ( other.fp_member_index_record [other_j ] )
        self.includefp =  other.includefp_record [other_j ]
        self.includefp_record [self_i] = other.includefp_record [other_j ]
        self.fp_involuntary =  other.fp_involuntary_record [other_j ]
        self.fp_involuntary_record [self_i] = other.fp_involuntary_record [other_j ]
        self.max_step_index_record[self_i] = other.max_step_index_record[other_j ]
        



# cluster_hard, cluster_soft
class Bunch2:
    def __init__(self, **kwargs):
        self.mixture_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.membership_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.membership_p_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.membership_p_normalize_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.likelihood_record = np.array ([  float("-inf") ] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))
        self.stepindex_record = np.array ([0] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))
        self.trialindex_record = np.array ([-1] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))
        self.makeone_index_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.fp_index_record = np.array ([-1] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))
        self.includefp_record = [False] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.fp_involuntary_record = [False] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.fp_member_index_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.max_step_index_record = np.array ([-1] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))

    def acc (self, mixture, membership, likelihood, membership_p, membership_p_normalize, step_index, trial_index, max_step_index, makeone_index, fp_index, includefp, fp_involuntary, fp_member_index, **kwargs):
        self.mixture = np.zeros ( (kwargs["NUM_BLOCK"], kwargs["NUM_CLONE"]), dtype = "float")
        self.mixture_record [kwargs["NUM_CLONE_NOMINAL"]] = copy.deepcopy  ( mixture ) 
        self.membership = np.zeros ( (kwargs["NUM_MUTATION"]), dtype = "int")
        self.membership_record [kwargs["NUM_CLONE_NOMINAL"]] = copy.deepcopy  ( membership )
        self.membership_p = np.zeros ( (kwargs["NUM_MUTATION"], kwargs["NUM_CLONE"]), dtype = "float")
        self.membership_p_record [kwargs["NUM_CLONE_NOMINAL"]] = copy.deepcopy  ( membership_p ) 
        self.membership_p_normalize = np.zeros ( (kwargs["NUM_MUTATION"], kwargs["NUM_CLONE"]), dtype = "float")
        self.membership_p_normalize_record [kwargs["NUM_CLONE_NOMINAL"]] = copy.deepcopy  ( membership_p_normalize )
        self.likelihood_record [kwargs["NUM_CLONE_NOMINAL"]] = likelihood
        self.stepindex = step_index
        self.stepindex_record [kwargs["NUM_CLONE_NOMINAL"]] = step_index
        self.trialindex= trial_index
        self.trialindex_record[kwargs["NUM_CLONE_NOMINAL"]]  = trial_index
        self.max_step_index = max_step_index
        self.max_step_index_record[kwargs["NUM_CLONE_NOMINAL"]]  = max_step_index
        self.makeone_index_record [kwargs["NUM_CLONE_NOMINAL"]] = copy.deepcopy  ( makeone_index )
        self.fp_index_record [kwargs["NUM_CLONE_NOMINAL"]] = fp_index
        self.includefp_record  [kwargs["NUM_CLONE_NOMINAL"]] = includefp
        self.fp_involuntary_record  [kwargs["NUM_CLONE_NOMINAL"]] = fp_involuntary
        self.fp_member_index_record  [kwargs["NUM_CLONE_NOMINAL"]] = copy.deepcopy  ( fp_member_index )

    def find_max_likelihood (self, start, end):
        i = np.argmax(self.likelihood_record [ start : end + 1]) + start
        return i

    def copy (self, other, self_i, other_j):
        other.mixture = copy.deepcopy  ( self.mixture_record [ self_i ] )
        other.mixture_record [other_j ] = copy.deepcopy  ( self.mixture_record [ self_i ] )
        other.likelihood = copy.deepcopy  ( self.likelihood_record [ self_i ] )
        other.likelihood_record [ other_j ] = copy.deepcopy  ( self.likelihood_record [ self_i ] )
        other.membership = copy.deepcopy  ( self.membership_record [ self_i ] )
        other.membership_record [other_j ] = copy.deepcopy  ( self.membership_record [ self_i ]  )
        other.membership_p_record [other_j ] = copy.deepcopy  ( self.membership_p_record [ self_i ] )
        #other.membership_p_normalize_record [other_j ] = self.membership_p_normalize_record [ self_i ] 
        other.makeone_index_record [ other_j ] = copy.deepcopy  ( self.makeone_index_record [ self_i ]  )
        other.fp_index = self.fp_index_record [ self_i ]
        other.fp_index_record [ other_j ] = copy.deepcopy  ( self.fp_index_record [ self_i ]  )
        other.stepindex =  self.stepindex
        other.includefp = self.includefp_record [ self_i ]
        other.includefp_record [ other_j ] = self.includefp_record [ self_i ]
        other.fp_involuntary = self.fp_involuntary_record [ self_i ]
        other.fp_involuntary_record [ other_j ] = self.fp_involuntary_record [ self_i ]
        other.fp_member_index = copy.deepcopy  ( self.fp_member_index_record [ self_i ] )
        other.fp_member_index_record [ other_j ] = copy.deepcopy  ( self.fp_member_index_record [ self_i ] )