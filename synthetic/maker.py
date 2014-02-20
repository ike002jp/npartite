#!/usr/bin/env python
# encoding: utf-8


#--------------------------------------------------
# modules
#--------------------------------------------------
from __future__ import division
import itertools
import copy
import random
import os.path
import pickle
import math

#--------------------------------------------------
# static variables
#--------------------------------------------------


#--------------------------------------------------
# private classes
#--------------------------------------------------

#--------------------------------------------------
# public classes
#--------------------------------------------------
class SyntheticNetwork(object):

    def __init__(self, name, config, eglist, has_problems):
        self._name = name
        self._config = config
        self._eglist = eglist
        self.has_problems = has_problems

    def name_with_config(self):
        conf = self._config
        return '%s(vrt-%d eg-%d noise-%.3f scale-%s)' % \
            (self._name, 
             conf.vrtnum, conf.egnum, 
             conf.noise_ratio, 
             conf.is_scale_free)

    def name(self):
        return '%s' % (self._name)

    def edge_list(self):
        return self._eglist

    def correct_community_labels(self):
        return self._config.community_index

    def is_same(self, syn_network):
        egset1 = set(self._eglist)
        egset2 = set( syn_network.get_edge_list() )
        return egset1 == egset2

    def save_as_edge_list(self, place, 
                          file_name_format='%(name)s.pickle',
                          name_dict={}, overwrite=False):

        # check the directory
        if not os.path.isdir(place):
            msg = '"%s" is not directory path.\n' % place
            msg += 'The third argument "place" should be '
            msg += 'a string of directory path.' 
            raise IOError(msg)

        # make file name
        _name_dict = {'name' : self._name}
        _name_dict.update(name_dict)
        file_name = file_name_format % _name_dict

        # check the path
        save_file_path = os.path.join(place, file_name)
        if (not overwrite) and os.path.isfile(save_file_path):
            msg = 'The file "%s" has already existed.' % save_file_path
            raise IOError(msg)

        # save
        f = open(save_file_path, 'w')
        pickle.dump(self._eglist, f)
        f.close()


class SyntheticNetworkMaker(object):

    def __init__(self, name, community_index, correct_corres_set, 
                       egnum, noise_ratio, is_scale_free):

        self._name = name
        self._config = _NetworkConfiguration(
                            community_index, correct_corres_set, 
                            egnum, noise_ratio, is_scale_free)
        self._making_status = _MakingStatus(self._config)

    def make(self, warn=False, until_succeed=True):

        if until_succeed:
            eglist = self._make_edge_list_until_succeed(warn)
            successful = True
        if not successful:
            eglist, successful = self._make_edge_list(warn)

        name = self._name
        config = self._config
        has_problems = not successful
        syn_network = SyntheticNetwork(name, config, eglist, has_problems)
        return syn_network

    def make_many(self, quantity, log=False,  warn=False):
        num = 0
        syn_networks = []

        while num < quantity:
            new_syn_network = self.make(warn=warn)

            # check that the new one does not exist 
            is_new = True
            for syn_network in syn_networks:
                if new_syn_network.is_same(syn_network):
                    is_new = False
                    break

            if is_new:
                syn_networks.append(new_syn_network)
                num += 1

            if log:
                self._print_proccess_log(self, num, quantity)

        return syn_networks

    def _make_edge_list(self, warn=False):
        self._making_status.init()

        eglist = []
        self._add_correct_edges(eglist)
        self._add_noise_edges(eglist)

        successful = not self._has_problems(eglist, warn)
        return eglist, successful

    def _make_edge_list_until_succeed(self, warn=False):
        eglist, successful = self._make_edge_list(warn)
        if not successful:
            return self._make_edge_list_until_succeed(warn)
        return eglist

    def _add_correct_edges(self, eglist):
        config = self._config
        added_egnum = 0
        while added_egnum < config.correct_egnum:
            is_added = self._add_edge(eglist, config.correct_corres_set)
            if is_added:
                added_egnum += 1

    def _add_noise_edges(self, eglist):
        config = self._config
        added_egnum = 0
        while added_egnum < config.noise_egnum:
            is_added = self._add_edge(eglist, config.noise_corres_set)
            if is_added:
                added_egnum += 1

    def _add_edge(self, eglist, corres_set):
        target_corres = self._choose_corres(corres_set)
        eg = self._choose_edge(target_corres)
        if eg in eglist:
            return False

        # add edge
        eglist.append(eg)

        # update status
        status = self._making_status
        for part, vrt in enumerate(eg):
            status.remove_vrt_from_zero_degree(part, vrt)
            status.plus_degree(part, vrt)

        return True

    def _choose_corres(self, corres_set):
        status = self._making_status
        
        if status.exist_zero_degree_vrt():
            return self._choose_corres_having_vrt_with_zerodeg(corres_set)
        else:
            correses = list(corres_set)
            return random.choice(correses)
        
    def _choose_corres_having_vrt_with_zerodeg(self, corres_set):
        max_vrtnum_with_zerodeg = 0
        candidate_correses = []
        for corres in corres_set:
            vrtnum_with_zerodeg = self._count_vrt_with_zerodeg_in_corres(corres)
                
            if vrtnum_with_zerodeg > max_vrtnum_with_zerodeg:
                max_vrtnum_with_zerodeg = vrtnum_with_zerodeg
                candidate_correses = [corres]
            elif vrtnum_with_zerodeg == max_vrtnum_with_zerodeg:
                candidate_correses.append(corres)
                
        return random.choice( candidate_correses )

    def _count_vrt_with_zerodeg_in_corres(self, corres):
        zero_degree_vrtnum = 0
        for part, com in enumerate(corres):
            zero_degree_vrtnum += self._count_vrt_with_zerodeg_in_com(part, com)
        
        return zero_degree_vrtnum
    
    def _count_vrt_with_zerodeg_in_com(self, part, com):
        config = self._config
        status = self._making_status
        
        memberset = config.memberset_from_com[part][com]
        zero_memberset = status.get_zero_degree_vrtset(part)
        return len( memberset & zero_memberset )
    
    def _choose_edge(self, corres):
        config = self._config

        # vertices which will be chose as an edge
        vrtset_list = self._make_vrtset_list_for_choosing_edge(corres)
        
        if config.is_scale_free:
            eg = [self._choose_vrt_in_proportion_to_degree(part, vrtset) 
                  for part, vrtset in enumerate(vrtset_list)]
        else:
            eg = [random.choice( list(vrtset) ) 
                  for vrtset in vrtset_list]
        
        return tuple(eg)
    
    def _make_vrtset_list_for_choosing_edge(self, corres):
        config = self._config
        status = self._making_status

        if not status.exist_zero_degree_vrt():
            vrtset_list = self._make_vrtset_list_from_corres(corres)
            return vrtset_list
        
        vrtset_list = []
        for part, com in enumerate(corres):
            vrtset_in_com = config.memberset_from_com[part][com]
            zero_vrtset = status.get_zero_degree_vrtset(part)
            vrtset = vrtset_in_com & zero_vrtset
            if len(vrtset) == 0:
                vrtset = vrtset_in_com
            vrtset_list.append(vrtset)
            
        return vrtset_list
    
    def _make_vrtset_list_from_corres(self, corres):
        config = self._config
        vrtset_list = [ set( config.memberset_from_com[part][com] )
                        for part, com in enumerate(corres) ]
        return vrtset_list

    def _choose_vrt_in_proportion_to_degree(self, part, vrtset):
        status = self._making_status
        
        degree_vrt_list = []
        for vrt in vrtset:
            degree = status.get_degree(part, vrt)
            if degree == 0:
                degree = 1
            degree_vrt_list.extend([vrt] * degree)
            
        random.shuffle(degree_vrt_list)
        return random.choice(degree_vrt_list)
    
    def _has_problems(self, edge_list, warn):
        
        is_duplicated = ( len(edge_list) != len( set(edge_list) ) )
        exist_zero_degree_vrt = self._making_status.exist_zero_degree_vrt()
        if is_duplicated and warn:
            print 'Edge list contains same edges.'

        if exist_zero_degree_vrt and warn:
            print 'There are some vertices whose degree is zero.'

        return (is_duplicated or exist_zero_degree_vrt)

    def _print_proccess_log(self, current_number, total_number):
        percent_step = 20
        percent_list = range(0, 100, percent_step)

        number_step = int( math.ceil(total_number * percent_step / 100) )
        target_number_list = range(0, total_number, number_step)
        target_number_list = target_number_list[:len(percent_list)]
        if current_number in target_number_list:
            ind = target_number_list.index(current_number)
            percent = percent_list[ind]
            print '%d%% has finished (%d networks have generated).' % (
                        percent, current_number)

class _NetworkConfiguration(object):
    def __init__(self, community_index,  correct_corres_set, 
                       egnum, noise_ratio, is_scale_free):
        self.partnum = len(community_index)

        # community
        self.community_index = community_index
        self.community_set = [ set(_community_index)
                               for _community_index in community_index]
        self.memberset_from_com = [{} for _ in xrange(self.partnum)]
        for part, _community_index in enumerate(self.community_index):
            for vrt_ind, com in enumerate(_community_index):
                self.memberset_from_com[part].setdefault(com,
                                                         set()).add(vrt_ind)

        # corresondency
        self.correct_corres_set = correct_corres_set
        self.noise_corres_set = self._make_noise_corres_set()

        # vertex
        self.vrtnum_list = [len(_community_index) 
                            for _community_index in community_index ]
        self.vrtnum = sum(self.vrtnum_list)
        self.vertex_index = [set( range(vrtnum) ) 
                             for vrtnum in self.vrtnum_list]

        # edge
        self.egnum = egnum
        self.noise_egnum = int(egnum * noise_ratio)
        self.correct_egnum = egnum - self.noise_egnum

        # 
        self.degree = egnum / (self.vrtnum / 3)
        self.noise_ratio = noise_ratio
        self.is_scale_free = is_scale_free

    def _make_noise_corres_set(self):
        all_corres_set = set( itertools.product(*self.community_set) )
        return all_corres_set - self.correct_corres_set


class _MakingStatus(object):

    def __init__(self, config):
        self._config = config
        self._zero_degree_vertex = None
        self._degree_from_vertex = None

    def init(self):
        config = self._config
        self._zero_degree_vertex = copy.deepcopy(config.vertex_index)
        self._degree_from_vertex = [{} for _ in xrange(config.partnum)]

    def get_zero_degree_vrtset(self, part):
        return self._zero_degree_vertex[part]

    def get_degree(self, part, vrt):
        return self._degree_from_vertex[part].setdefault(vrt, 0)

    def plus_degree(self, part, vrt):
        degree = self._degree_from_vertex[part].setdefault(vrt, 0)
        self._degree_from_vertex[part][vrt] = degree + 1
        #return self._degree_from_vertex[part][vrt]

    def remove_vrt_from_zero_degree(self, part, vrt):
        self._zero_degree_vertex[part].discard(vrt)

    def exist_zero_degree_vrt(self):
        zero_vrtnum = 0
        for _zero_degree_vertex in self._zero_degree_vertex:
            zero_vrtnum += len(_zero_degree_vertex)

        if zero_vrtnum == 0:
            return False
        return True

if __name__ == '__main__':
    partnum = 3
    community_index = [ [0] * 5 + [1] * 5
                        for _ in xrange(partnum)]
    correct_corres_set = set([ (0, 0, 0), (1, 1, 1) ])
    egnum = 30
    noise_ratio = 0.1
    #is_scale_free = False
    is_scale_free = True
    syn_maker = SyntheticNetworkMaker(
                    'Test Network',
                    community_index, correct_corres_set, 
                    egnum, noise_ratio, is_scale_free)
    warn = True
    syn_network = syn_maker.make(warn=warn)
    edge_list = syn_network.edge_list()
    edge_list.sort()
    com_labels = syn_network.correct_community_labels()
    print 'basic usage'
    print syn_network.name()
    print syn_network.name_with_config()
    print com_labels
    print edge_list
    print syn_network.has_problems
    print ''
    
    '''
    print 'save network'
    place = './'
    syn_network.save_as_edge_list(place, overwrite=False)
    '''

