import numpy as np
import copy
import pandas as pd
import EMhard


class Bunch1:        # step_hard, step_soft, trial_hard, trial_soft
    def __init__(self, NUM_MUTATION, NUM_BLOCK, NUM_CLONE, K):
        self.mixture = np.zeros ( (NUM_BLOCK, NUM_CLONE + 1), dtype = "float")
        self.mixture_record = np.zeros ( (K, NUM_BLOCK, NUM_CLONE + 1) , dtype = "float")
        self.initial_randomsample = []
        self.initial_randomsample_record = [[] for _ in range(K)]  # [[]] * K 
        self.initial_sampling = True
        self.membership = np.zeros ( (NUM_MUTATION), dtype = "int")
        self.membership_record = np.zeros ( (K, NUM_MUTATION), dtype = "int")
        self.membership_p = np.zeros ( (NUM_MUTATION, NUM_CLONE + 1) , dtype = "float")
        self.membership_p_record = np.zeros ( (K, NUM_MUTATION, NUM_CLONE + 1), dtype = "float")
        self.membership_p_normalize = np.zeros ( (NUM_MUTATION, NUM_CLONE + 1), dtype = "float")
        self.membership_p_normalize_record = np.zeros ( (K, NUM_MUTATION, NUM_CLONE + 1), dtype = "float")
        self.posterior = float("-inf")
        self.posterior_record = np.array ([  float("-inf") ] * (K))
        self.likelihood = float("-inf")
        self.likelihood_record = np.array ([  float("-inf") ] * (K))
        self.posterior_normalized = float("-inf")
        self.posterior_normalized_record = np.array ([  float("-inf") ] * (K))
        self.likelihood_normalized = float("-inf")
        self.likelihood_normalized_record = np.array ([  float("-inf") ] * (K))
        self.stepindex = 0
        self.stepindex_record = np.zeros (K , dtype = "int")
        self.max_step_index = -1
        self.max_step_index_record = np.array ([-1] * (K))
        self.makeone_index = list (range (NUM_CLONE))
        self.makeone_index_record =  [[] for _ in range(K)]   # [[]] * K 
        self.fp_index = -1
        self.fp_index_record = np.array ([-1] * (K))
        self.fp_member_index = []
        self.fp_member_index_record = [[] for _ in range(K)]
        self.includefp = False
        self.includefp_record = np.array ( [False] * (K))
        self.fp_involuntary = False
        self.fp_involuntary_record = np.array ( [False] * (K))
        self.tn_index = []
        self.tn_index_record =  [[] for _ in range(K)]
        self.lowvafcluster_index =   [[] for _ in range(NUM_BLOCK)]
        self.lowvafcluster_index_record =  [[] for _ in range(K)]
        self.lowvafcluster_cdf_list =   [[] for _ in range(NUM_BLOCK)]
        self.cluster_cdf_list =   [[] for _ in range(NUM_BLOCK)]
        self.checkall_strict = True
        self.checkall_strict_record = np.zeros (K , dtype = "bool")
        self.checkall_lenient = True
        self.checkall_lenient_record = np.zeros (K , dtype = "bool")
        self.makeone_experience = False
        self.hard_or_soft = "hard"
        self.less_than_min_cluster_size = False
        self.df = pd.DataFrame ( columns = ["NUM_BLOCK", "TRIAL", "STEP", "STEP_TOTAL", "POSTERIOR", "LIKELIHOOD", "HARD_SOFT"] ) 

    def acc (self, mixture, membership, posterior, likelihood, membership_p, membership_p_normalize, makeone_index, tn_index, fp_index, step_index, fp_member_index, lowvafcluster_index, includefp, fp_involuntary, checkall_strict, checkall_lenient, max_step_index, K):
        self.mixture_record [K]= copy.deepcopy ( mixture )
        self.membership_record [K] = copy.deepcopy ( membership ) 
        self.posterior_record [K] = posterior
        self.likelihood_record [K] = likelihood
        self.membership_p_record [K] = copy.deepcopy ( membership_p )
        self.membership_p_normalize_record [K] = copy.deepcopy ( membership_p_normalize )
        self.makeone_index_record[K] = copy.deepcopy ( makeone_index )
        self.tn_index_record[K] = copy.deepcopy ( tn_index )
        self.fp_index_record[K] = copy.deepcopy ( fp_index )
        self.stepindex = step_index
        self.stepindex_record[K] = step_index
        self.max_step_index = max_step_index
        self.max_step_index_record [K] = max_step_index
        self.fp_member_index = copy.deepcopy ( fp_member_index )
        self.fp_member_index_record [K] = copy.deepcopy ( fp_member_index )
        self.lowvafcluster_index_record [K] = copy.deepcopy ( lowvafcluster_index )
        self.includefp = includefp
        self.includefp_record [K] = includefp
        self.fp_involuntary  = fp_involuntary 
        self.fp_involuntary_record [K] = fp_involuntary 
        self.checkall_strict = checkall_strict
        self.checkall_strict_record [K] = checkall_strict
        self.checkall_lenient = checkall_lenient
        self.checkall_lenient_record [K] = checkall_lenient


    def find_max_likelihood_nocondition (self, start, end):         # For 
        if start > end:
            return -1, False
        
        max, max_index = -9999999, -1
        for i in range (start, end + 1):
            if self.likelihood_record [i] > max:
                max = self.likelihood_record [i]
                max_index = i
        
        if max == -9999999:
            return -1, False
        
        return max_index, self.checkall_strict_record [max_index]
    
    def find_max_likelihood_strictfirst_fromtheend (self, start, end):   # step
        if start > end:
            return -1, False
        
        max, max_index = -9999999, -1
        for i in range (end, start - 1, -1):
            if self.checkall_strict_record[i] == True:   # 뒤에서부터 돌면서 strict를 return
                return i, True

        return end, True             # lenient라도 맨 마지막을 선택

    def find_max_likelihood_strictfirst (self, start, end):   # trial
        if start > end:
            return -1, False
        
        max, max_index = -9999999, -1
        for i in range (start, end):
            if self.checkall_strict_record[i] == True:   # Qualified trial only
                if self.likelihood_record [i] > max:
                    max = self.likelihood_record [i]
                    max_index = i
        if max_index != -1 :   # checkall(strict) == True를 찾으면 그걸 보내줌
            return max_index, True
        
        # checkall(lenient) == True라도 찾자
        max, max_index = -9999999, -1
        for i in range (start, end ):
            if self.checkall_lenient_record[i] == True:   # Qualified trial only
                if self.likelihood_record [i] > max:
                    max = self.likelihood_record [i]
                    max_index = i
        if max_index != -1:
            return max_index, False
        
         # 아예 모든게 -inf, -inf, -inf일때
        return -1, False
            
            

    def copy (self, other, self_i, other_j):  # step_soft <- cluster_hard
        self.mixture = copy.deepcopy ( other.mixture_record [ other_j ] )
        self.mixture_record [self_i] = copy.deepcopy ( other.mixture_record[other_j] )
        self.membership = copy.deepcopy ( other.membership_record [ other_j ] )
        self.membership_record [self_i] = copy.deepcopy ( other.membership_record[ other_j ] )
        self.membership_p = copy.deepcopy  ( other.membership_p_record[ other_j ] ) 
        self.membership_p_record [self_i] = copy.deepcopy ( other.membership_p_record[ other_j ] )
        self.membership_p_normalize_record [self_i] = copy.deepcopy ( other.membership_p_normalize_record[ other_j ] )
        self.posterior = copy.deepcopy ( other.posterior_record [ other_j ] )
        self.posterior_record [ self_i ] = copy.deepcopy ( other.posterior_record [other_j] )
        self.likelihood = copy.deepcopy ( other.likelihood_record [ other_j ] )
        self.likelihood_record [ self_i ] = copy.deepcopy ( other.likelihood_record [other_j] )
        self.makeone_index = copy.deepcopy  ( other.makeone_index_record[ other_j ] )
        self.makeone_index_record [self_i] = copy.deepcopy ( other.makeone_index_record[ other_j ] )
        self.tn_index = copy.deepcopy  ( other.tn_index_record[ other_j ] )
        self.tn_index_record [self_i] = copy.deepcopy ( other.tn_index_record[ other_j ] )
        self.fp_index = other.fp_index_record[ other_j ]
        self.fp_index_record [self_i] = copy.deepcopy ( other.fp_index_record[ other_j ] )
        self.fp_member_index = other.fp_member_index_record [other_j ]
        self.fp_member_index_record [self_i] = copy.deepcopy ( other.fp_member_index_record [other_j ] )
        self.includefp =  other.includefp_record [other_j ]
        self.includefp_record [self_i] = other.includefp_record [other_j ]
        self.fp_involuntary =  other.fp_involuntary_record [other_j ]
        self.fp_involuntary_record [self_i] = other.fp_involuntary_record [other_j ]
        self.max_step_index_record[self_i] = other.max_step_index_record[other_j ]
        self.checkall_strict = other.checkall_strict_record[other_j]
        self.checkall_strict_record[self_i] = other.checkall_strict_record[other_j]
        self.checkall_lenient = other.checkall_lenient_record[other_j]
        self.checkall_lenient_record[self_i] = other.checkall_lenient_record[other_j]


        



# cluster_hard, cluster_soft
class Bunch2:
    def __init__(self, **kwargs):
        self.mixture_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.membership_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.membership_p_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.membership_p_normalize_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.posterior_record = np.array ([  float("-inf") ] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))
        self.likelihood_record = np.array ([  float("-inf") ] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))
        self.stepindex_record = np.array ([0] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))
        self.trialindex_record = np.array ([-1] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))
        self.makeone_index_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1) 
        self.checkall_strict_record = [False] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.checkall_lenient_record = [False] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.fp_index_record = np.array ([-1] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))
        self.includefp_record = [False] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.fp_involuntary_record = [False] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.fp_member_index_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.tn_index_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.max_step_index_record = np.array ([-1] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))
        self.df = pd.DataFrame ( columns =  ["NUM_BLOCK", "TRIAL", "STEP", "STEP_TOTAL", "POSTERIOR", "LIKELIHOOD", "HARD_SOFT", "NUM_CLONE"] ) 

    def acc (self, mixture, membership, posterior, likelihood, membership_p, membership_p_normalize, step_index, trial_index, max_step_index, makeone_index, trial_checkall_strict, trial_checkall_lenient, tn_index, fp_index, includefp, fp_involuntary, fp_member_index, **kwargs):
        self.mixture = np.zeros ( (kwargs["NUM_BLOCK"], kwargs["NUM_CLONE"] + 1), dtype = "float")
        self.mixture_record [kwargs["NUM_CLONE_ITER"]] = copy.deepcopy  ( mixture ) 
        self.membership = np.zeros ( (kwargs["NUM_MUTATION"]), dtype = "int")
        self.membership_record [kwargs["NUM_CLONE_ITER"]] = copy.deepcopy  ( membership )
        self.membership_p = np.zeros ( (kwargs["NUM_MUTATION"], kwargs["NUM_CLONE"] + 1), dtype = "float")
        self.membership_p_record [kwargs["NUM_CLONE_ITER"]] = copy.deepcopy  ( membership_p ) 
        self.membership_p_normalize = np.zeros ( (kwargs["NUM_MUTATION"], kwargs["NUM_CLONE"] + 1), dtype = "float")
        self.membership_p_normalize_record [kwargs["NUM_CLONE_ITER"]] = copy.deepcopy  ( membership_p_normalize )
        self.posterior_record [kwargs["NUM_CLONE_ITER"]] = posterior
        self.likelihood_record [kwargs["NUM_CLONE_ITER"]] = likelihood
        self.stepindex = step_index
        self.stepindex_record [kwargs["NUM_CLONE_ITER"]] = step_index
        self.trialindex= trial_index
        self.trialindex_record[kwargs["NUM_CLONE_ITER"]]  = trial_index
        self.max_step_index = max_step_index
        self.max_step_index_record[kwargs["NUM_CLONE_ITER"]]  = max_step_index
        self.makeone_index_record [kwargs["NUM_CLONE_ITER"]] = copy.deepcopy  ( makeone_index )
        self.checkall_strict_record [kwargs["NUM_CLONE_ITER"]] = trial_checkall_strict
        self.checkall_lenient_record [kwargs["NUM_CLONE_ITER"]] = trial_checkall_lenient
        self.fp_index_record [kwargs["NUM_CLONE_ITER"]] = fp_index
        self.includefp_record  [kwargs["NUM_CLONE_ITER"]] = includefp
        self.fp_involuntary_record  [kwargs["NUM_CLONE_ITER"]] = fp_involuntary
        self.fp_member_index_record  [kwargs["NUM_CLONE_ITER"]] = copy.deepcopy  ( fp_member_index )
        self.tn_index_record [kwargs["NUM_CLONE_ITER"]] = copy.deepcopy  ( tn_index )

    def find_max_likelihood (self, start, end):
        i = np.argmax(self.likelihood_record [ start : end + 1]) + start
        return i

    def copy (self, other, self_i, other_j):
        other.mixture = copy.deepcopy  ( self.mixture_record [ self_i ] )
        other.mixture_record [other_j ] = copy.deepcopy  ( self.mixture_record [ self_i ] )
        other.posterior = copy.deepcopy  ( self.posterior_record [ self_i ] )
        other.posterior_record [ other_j ] = copy.deepcopy  ( self.posterior_record [ self_i ] )
        other.likelihood = copy.deepcopy  ( self.likelihood_record [ self_i ] )
        other.likelihood_record [ other_j ] = copy.deepcopy  ( self.likelihood_record [ self_i ] )
        other.membership = copy.deepcopy  ( self.membership_record [ self_i ] )
        other.membership_record [other_j ] = copy.deepcopy  ( self.membership_record [ self_i ]  )
        other.membership_p_record [other_j ] = copy.deepcopy  ( self.membership_p_record [ self_i ] )
        #other.membership_p_normalize_record [other_j ] = self.membership_p_normalize_record [ self_i ] 
        other.makeone_index_record [ other_j ] = copy.deepcopy  ( self.makeone_index_record [ self_i ]  )
        other.checkall_strict_record [ other_j ] = copy.deepcopy  ( self.checkall_strict_record [ self_i ]  )
        other.tn_index_record [ other_j ] = copy.deepcopy  ( self.tn_index_record [ self_i ]  )
        other.fp_index = self.fp_index_record [ self_i ]
        other.fp_index_record [ other_j ] = copy.deepcopy  ( self.fp_index_record [ self_i ]  )
        other.stepindex =  self.stepindex
        other.includefp = self.includefp_record [ self_i ]
        other.includefp_record [ other_j ] = self.includefp_record [ self_i ]
        other.fp_involuntary = self.fp_involuntary_record [ self_i ]
        other.fp_involuntary_record [ other_j ] = self.fp_involuntary_record [ self_i ]
        other.fp_member_index = copy.deepcopy  ( self.fp_member_index_record [ self_i ] )
        other.fp_member_index_record [ other_j ] = copy.deepcopy  ( self.fp_member_index_record [ self_i ] )
        other.checkall_strict_record [ other_j ] = self.checkall_strict_record [self_i]
        other.checkall_lenient_record [ other_j ] = self.checkall_lenient_record [self_i]