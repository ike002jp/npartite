#!/usr/bin/env python
# encoding: utf-8

'''
'''

#--------------------------------------------------
# modules
#--------------------------------------------------
from __future__ import division

import math
try:
    import numpy
except ImportError:
    # for PyPy environment
    import numpypy as numpy

#--------------------------------------------------
# static variables
#--------------------------------------------------
RADIX = 2

#--------------------------------------------------
# public functions
#--------------------------------------------------
def calculate_nmi(ind_list1, ind_list2):
    flatten_list1 = _flatten(ind_list1)
    flatten_list2 = _flatten(ind_list2)
    return _calculate_nmi(flatten_list1, flatten_list2)

#--------------------------------------------------
# private functions
#--------------------------------------------------
def _flatten(multi_ind_list):
    flatten_list = []
    
    max_ind = 0
    map_dic = {}
    for part, _ind_list in enumerate(multi_ind_list):
        for ind in _ind_list:
            key = "%d-%d" % (part, ind)
            if key not in map_dic:
                map_dic[key] = max_ind
                max_ind += 1
                
            map_ind = map_dic[key]
            flatten_list.append(map_ind)
    
    return flatten_list


def _calculate_nmi(ind_list1, ind_list2):
    element_num1 = len(ind_list1)
    element_num2 = len(ind_list2)
    
    if element_num1 != element_num2:
        msg = 'The number of elements of the two arguments should be equal'
        raise ValueError, msg
        
    N = _make_confusion_matrix(ind_list1, ind_list2)
    molecule = _calculate_molecule(N)
    denominator = _calculate_denominator(N)
    if molecule == 0:
        nmi_val = 0
    else:
        nmi_val = molecule / denominator
    
    return float(nmi_val)

def _make_confusion_matrix(ind_list1, ind_list2):
    element_set_from_com_1 = _convert_to_dic(ind_list1)
    element_set_from_com_2 = _convert_to_dic(ind_list2)
    class_num1 = len(element_set_from_com_1)
    class_num2 = len(element_set_from_com_2)
    
    N = numpy.zeros((class_num1, class_num2))
    for row in xrange(class_num1):
        element_set1 = element_set_from_com_1[row]
        
        for col in xrange(class_num2):
            element_set2 = element_set_from_com_2[col]
            N[row, col] = len( element_set1 & element_set2 )
            
    return N
    
def _calculate_molecule(N):
    class_num1, class_num2 = N.shape
    n = N.sum()
    rtn = 0
    for i in xrange(class_num1):
        N_i_sum = N[i,:].sum()
        for j in xrange(class_num2):
            N_ij = N[i,j]
            N_j_sum = N[:,j].sum()
            if not N_ij == 0:
                log_term = math.log( ((N_ij * n) / (N_i_sum * N_j_sum)), RADIX)
                #log_term = math.log10( (N_ij * n) / (N_i_sum * N_j_sum) )
                rtn += N_ij * log_term
    return -2 * rtn
    
def _calculate_denominator(N):
    class_num1, class_num2 = N.shape
    n = N.sum()
    
    term1 = 0
    for i in xrange(class_num1):
        N_i_sum = N[i,:].sum()
        if not N_i_sum == 0:
            term1 += N_i_sum * math.log( (N_i_sum / n), RADIX)
            #term1 += N_i_sum * math.log10(N_i_sum / n)
        
    term2 = 0
    for j in xrange(class_num2):
        N_j_sum = N[:,j].sum()
        if not N_j_sum == 0:
            term2 += N_j_sum * math.log( (N_j_sum / n), RADIX)
            #term2 += N_j_sum * math.log10(N_j_sum / n)
    return term1 + term2

def _convert_to_dic(ind_list):
    # make mapping dic
    map_dic = {}
    max_ind = 0
    for class_ind in ind_list:
        if class_ind not in map_dic.keys():
            map_dic[class_ind] = max_ind
            max_ind += 1
    
    # make dic
    dic = {}
    for element_ind, class_ind in enumerate(ind_list):
        map_class_ind = map_dic[class_ind]
        dic.setdefault(map_class_ind, set()).add(element_ind)
        
    return dic


#--------------------------------------------------
# test
#--------------------------------------------------
if __name__ == "__main__":
    # test 1
    # answer : 1.0
    ind_list1 = range(5)
    ind_list2 = range(5)
    nmi = calculate_nmi(ind_list1, ind_list2)
    answer = 1.0
    tn = answer == round(nmi, 13)
    print "    %.5f, %s" % (nmi, tn)
    
    # test 2
    # answer : 1.0
    ind_list1 = [1, 1, 2, 2, 3]
    ind_list2 = [1, 1, 2, 2, 4]
    nmi = calculate_nmi(ind_list1, ind_list2)
    answer = 1.0
    tn = answer == round(nmi, 13)
    print "    %.5f, %s" % (nmi, tn)
    
    # test 3
    # answer : 0.737175493807
    ind_list1 = [1, 1, 2, 2, 3]
    ind_list2 = [1, 1, 2, 4, 4]
    nmi = calculate_nmi(ind_list1, ind_list2)
    answer = 0.737175493807
    tn = answer == round(nmi, 13)
    print "    %.5f, %s" % (nmi, tn)


