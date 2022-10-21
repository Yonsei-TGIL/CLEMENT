import numpy as np
import copy

# step : 각 step에서의 값을 저장
# trial : 각 trial에서의 best 값을 저장
class Bunch1:
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
        self.outlier_index = []
        self.outlier_index_record = [[]] * K
        self.includeoutlier = False
        self.includeoutlier_record = np.array ( [False] * (K))

    def acc (self, mixture, membership, likelihood, membership_p, membership_p_normalize, makeone_index, fp_index, step_index, outlier_index, includeoutlier, max_step_index, K):
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
        self.outlier_index = copy.deepcopy ( outlier_index )
        self.outlier_index_record [K] = copy.deepcopy ( outlier_index )
        self.includeoutlier = includeoutlier
        self.includeoutlier_record [K] = includeoutlier


    def find_max_likelihood (self, start, end):
        try:
            i = np.argmax(self.likelihood_record [ start : end + 1 ]) + start
            return i
        except:
            return -1
        

    def copy (self, other, self_i, other_j):
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
        self.outlier_index = other.outlier_index_record [other_j ]
        self.outlier_index_record [self_i] = copy.deepcopy ( other.outlier_index_record [other_j ] )
        self.includeoutlier =  other.includeoutlier_record [other_j ]
        self.includeoutlier_record [self_i] = other.includeoutlier_record [other_j ]
        self.max_step_index_record[self_i] = other.max_step_index_record[other_j ]



# self :  각 NUM_CLONE에서의 best 값을 저장  (record가 중요)
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
        self.includeoutlier_record = [False] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.outlier_index_record = [[]] * (kwargs["NUM_CLONE_TRIAL_END"] + 1)
        self.max_step_index_record = np.array ([-1] * (kwargs["NUM_CLONE_TRIAL_END"] + 1))

    def acc (self, mixture, membership, likelihood, membership_p, membership_p_normalize, step_index, trial_index, max_step_index, makeone_index, fp_index, includeoutlier, outlier_index, **kwargs):
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
        self.includeoutlier_record  [kwargs["NUM_CLONE_NOMINAL"]] = includeoutlier
        self.outlier_index_record  [kwargs["NUM_CLONE_NOMINAL"]] = copy.deepcopy  ( outlier_index )

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
        other.includeoutlier = self.includeoutlier_record [ self_i ]
        other.includeoutlier_record [ other_j ] = self.includeoutlier_record [ self_i ]
        other.outlier_index = copy.deepcopy  ( self.outlier_index_record [ self_i ] )
        other.outlier_index_record [ other_j ] = copy.deepcopy  ( self.outlier_index_record [ self_i ] )