#!/usr/bin/env python
# encoding: utf-8


#--------------------------------------------------
# modules
#--------------------------------------------------
from __future__ import division

from maker import SyntheticNetworkMaker

#--------------------------------------------------
# static variables
#--------------------------------------------------


#--------------------------------------------------
# private classes
#--------------------------------------------------

#--------------------------------------------------
# public classes
#--------------------------------------------------
class SimpleCaseMaker(SyntheticNetworkMaker):

    def __init__(self, corres_egnum, noise_ratio, 
                 is_scale_free=False, community_size=10):
        name = 'Simple Case'
        comsize = community_size
        community_index = [ [0] * comsize + [1] * comsize,
                            [0] * comsize + [1] * comsize,
                            [0] * comsize + [1] * comsize ]
        correct_corres_set = set([ (0, 0, 0), (1, 1, 1) ])
        egnum = corres_egnum * 2
        SyntheticNetworkMaker.__init__(self, name, 
                                       community_index, correct_corres_set, 
                                       egnum, noise_ratio, is_scale_free)

class OverlappingCaseMaker(SyntheticNetworkMaker):

    def __init__(self, corres_egnum, noise_ratio, 
                 is_scale_free=False, community_size=10):
        name = 'Overlapping Case'
        comsize = community_size
        community_index = [ [0] * comsize + [1] * comsize,
                            [0] * comsize,
                            [0] * comsize + [1] * comsize ]
        correct_corres_set = set([ (0, 0, 0), (1, 0, 1) ])
        egnum = corres_egnum * 2
        SyntheticNetworkMaker.__init__(self, name, 
                                       community_index, correct_corres_set, 
                                       egnum, noise_ratio, is_scale_free)

class ContradictiveCaseMaker(SyntheticNetworkMaker):

    def __init__(self, corres_egnum, noise_ratio, 
                 is_scale_free=False, community_size=10):
        name = 'Contradictive Case'
        comsize = community_size
        community_index = [ [0] * comsize + [1] * comsize,
                            [0] * comsize + [1] * comsize,
                            [0] * comsize + [1] * comsize ]
        correct_corres_set = set([ (0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0) ])
        egnum = corres_egnum * 4
        SyntheticNetworkMaker.__init__(self, name, 
                                       community_index, correct_corres_set, 
                                       egnum, noise_ratio, is_scale_free)

class SimplePlusCaseMaker(SyntheticNetworkMaker):

    def __init__(self, corres_egnum, noise_ratio, 
                 is_scale_free=False, community_size=10):
        name = 'Simpe Plus Case'
        comsize = community_size
        community_index = [ [0] * comsize + [1] * comsize,
                            [0] * comsize + [1] * comsize,
                            [0] * comsize + [1] * comsize ]
        correct_corres_set = set([ (0, 0, 0), (0, 1, 0), (1, 1, 1) ])
        egnum = corres_egnum * 3
        SyntheticNetworkMaker.__init__(self, name, 
                                       community_index, correct_corres_set, 
                                       egnum, noise_ratio, is_scale_free)

class Simple3CaseMaker(SyntheticNetworkMaker):

    def __init__(self, corres_egnum, noise_ratio, 
                 is_scale_free=False, community_size=10):
        name = 'Simple 3 Case'
        comsize = community_size
        community_index = [ [0] * comsize + [1] * comsize + [2] * comsize,
                            [0] * comsize + [1] * comsize + [2] * comsize,
                            [0] * comsize + [1] * comsize + [2] * comsize ]
        correct_corres_set = set([ (0, 0, 0), (1, 1, 1), (2, 2, 2) ])
        egnum = corres_egnum * 3
        SyntheticNetworkMaker.__init__(self, name, 
                                       community_index, correct_corres_set, 
                                       egnum, noise_ratio, is_scale_free)

class Simple3PlusCaseMaker(SyntheticNetworkMaker):

    def __init__(self, corres_egnum, noise_ratio, 
                 is_scale_free=False, community_size=10):
        name = 'Simple 3 Plus Case'
        comsize = community_size
        community_index = [ [0] * comsize + [1] * comsize + [2] * comsize,
                            [0] * comsize + [1] * comsize + [2] * comsize,
                            [0] * comsize + [1] * comsize + [2] * comsize ]
        correct_corres_set = set([ (0, 0, 0), (0, 1, 0), 
                                   (1, 1, 1), (1, 2, 1), (2, 2, 2) ])
        egnum = corres_egnum * 5
        SyntheticNetworkMaker.__init__(self, name, 
                                       community_index, correct_corres_set, 
                                       egnum, noise_ratio, is_scale_free)


def _make_and_print(syn_maker, warn):
    syn_network = syn_maker.make(warn=warn)

    com_labels = syn_network.correct_community_labels()
    eglist = syn_network.edge_list()
    eglist.sort()

    print syn_network.name()
    print syn_network.name_with_config()
    print eglist
    print com_labels
    print ''


if __name__ == '__main__':
    
    print 'typical network usage'
    corres_egnum = 30
    noise_ratio = 0.1
    is_scale_free = False
    warn = False
    syn_maker = SimpleCaseMaker(corres_egnum, noise_ratio, is_scale_free)
    _make_and_print(syn_maker, warn)

    syn_maker = OverlappingCaseMaker(corres_egnum, noise_ratio, is_scale_free)
    _make_and_print(syn_maker, warn)

    syn_maker = ContradictiveCaseMaker(corres_egnum, noise_ratio, is_scale_free)
    _make_and_print(syn_maker, warn)

    syn_maker = SimplePlusCaseMaker(corres_egnum, noise_ratio, is_scale_free)
    _make_and_print(syn_maker, warn)

    syn_maker = Simple3CaseMaker(corres_egnum, noise_ratio, is_scale_free)
    _make_and_print(syn_maker, warn)

    syn_maker = Simple3PlusCaseMaker(corres_egnum, noise_ratio, is_scale_free)
    _make_and_print(syn_maker, warn)


