#!/usr/bin/env python
# encoding: utf-8


#--------------------------------------------------
# modules
#--------------------------------------------------
from __future__ import division
import itertools
import copy
import random

#--------------------------------------------------
# static variables
#--------------------------------------------------


#--------------------------------------------------
# private classes
#--------------------------------------------------
class _BasicInformation(object):

    def __init__(self, edges):

        if not isinstance(edges, tuple):
            raise ValueError("Edges should be tuple.")

        self._edges = edges

        # parameters
        _n_part_vertex_list = zip(*edges)
        self._vertex_num_list = [len( set(vertex_list) ) 
                                    for vertex_list in _n_part_vertex_list]
        self._total_vertex_num = sum(self._vertex_num_list)
        self._edge_num = len(edges)
        self._partnum = len(edges[0])
        self._adj_egset_from_vrt = \
            self._make_adjacency_between_vertex_and_edge(edges)

        self._degrees = [{} for _ in range(self._partnum)]
        for part, adj_egset_from_vrt in enumerate(self._adj_egset_from_vrt):
            for vrt, adj_egset in adj_egset_from_vrt.iteritems():
                self._degrees[part][vrt] = len(adj_egset)

    def _make_adjacency_between_vertex_and_edge(self, edges):
        between_vrt_eg = [{} for _ in range(self._partnum)]
        for eg_ind, eg in enumerate(edges):
            for part, vertex in enumerate(eg):
                between_vrt_eg[part].setdefault(vertex, set()).add(eg_ind)
            
        return tuple(between_vrt_eg)

    def edgenum(self):
        return self._edge_num

    def partnum(self):
        return self._partnum

    def get_vertex_num_list(self):
        return self._vertex_num_list

    def get_edges_as_tuple(self):
        return self._edges

    def get_edge(self, eg_ind):
        return self._edges[eg_ind]

    def adj_egset_to_vrt(self, part, vrt):
        return self._adj_egset_from_vrt[part][vrt]

    def adj_row_egset_to_vrt(self, part, vrt):
        row_egset = set()
        for eg_ind in self._adj_egset_from_vrt[part][vrt]:
            row_egset.add( self.get_edge(eg_ind) )
        return row_egset

    def adj_egset_to_eg(self, eg_ind):
        eg = self._edges[eg_ind]
        adj_egset = set()
        for part, vrt in enumerate(eg):
            adj_egset.update(self._adj_egset_from_vrt[part][vrt])

        return adj_egset

    def degree(self, part, vrt):
        return self._degrees[part][vrt]
    

class _CommunityStatus(object):

    def __init__(self, basic_info):
        self._basic_info = basic_info
        self._partnum = basic_info.partnum()

        _vrtnum_list = self._basic_info.get_vertex_num_list()
        self._com_labels = [ [None] * _vrtnum 
                                  for _vrtnum in _vrtnum_list ]
 
        self._memberset_in_com = None
        self._membernum_in_com = None
        self._egset_from_com = None
        self._egnum_from_com = None
        self._corresset_from_com =None
        self._egset_from_corres = None
        self._egnum_from_corres = None

    # accessor : whole community structures
    def count_all_comnums(self):
        return sum([len(__memberset_in_com) 
                    for __memberset_in_com in self._memberset_in_com])

    def get_com_labels(self):
        return self._com_labels

    # accessor : each community
    def com_of_vrt(self, part=None, vrt=None):
        return self._com_labels[part][vrt]

    def memberset_in_com(self, part, com):
        #return self._memberset_in_com[part].get(com, set())
        return self._memberset_in_com[part].setdefault(com, set())

    def membernum_in_com(self, part, com):
        #return self._membernum_in_com[part].get(com, 0)
        return self._membernum_in_com[part].setdefault(com, 0)

    def egset_from_com(self, part, com):
        #return self._egset_from_com[part].get(com, set())
        return self._egset_from_com[part].setdefault(com, set())

    def egnum_from_com(self, part, com):
        #return self._egnum_from_com[part].get(com, 0)
        return self._egnum_from_com[part].setdefault(com, 0)

    def corresset_from_com(self, part, com):
        #return self._corresset_from_com[part].get(com, set())
        return self._corresset_from_com[part].setdefault(com, set())

    def iter_com_combinations(self, part, combi_num=2):
        coms = self._memberset_in_com[part].keys()
        return itertools.combinations(coms, combi_num)

    # accessor : correspondency
    def iter_corres_egnum(self):
        return self._egnum_from_corres.iteritems()

    def egset_from_corres(self, corres):
        return self._egset_from_corres.get(corres, set())
        #return self._egset_from_corres.setdefault(corres, set())

    def egnum_from_corres(self, corres):
        return self._egnum_from_corres.get(corres, 0)
        #return self._egnum_from_corres.setdefault(corres, 0)

    # manipulator : community index
    def assign_unique_com_labels(self):
        vertex_num_list = self._basic_info.get_vertex_num_list()
        com_labels = [range(num) for num in vertex_num_list]
        self.set_com_labels(com_labels)

    def set_com_labels(self, com_labels):
        self._com_labels = com_labels
        self._update_community_structure()

    # inner manipulator : used in this class
    def _del_com(self, part, com):
        del self._memberset_in_com[part][com]
        del self._membernum_in_com[part][com]
        del self._egset_from_com[part][com]
        del self._egnum_from_com[part][com]
        del self._corresset_from_com[part][com]

    def _del_corres(self, corres):
        del self._egset_from_corres[corres]
        del self._egnum_from_corres[corres]
        for part, com in enumerate(corres):
            if com in self._corresset_from_com[part]:
                self._corresset_from_com[part][com].discard(corres)

    def _assign_com(self, part, vrt, com):
        self._com_labels[part][vrt] = com

    # maniurator : merge, move
    def merge_coms(self, part, com1, com2):
        moved_vrtset = self._memberset_in_com[part][com1]
        partnum = self._partnum
        moved_vrts_info = [{} for _ in range(partnum)]
        for vrt in moved_vrtset:
            moved_vrts_info[part][vrt] = [com1, com2]

        moving_diff_info = self.diff_of_moving_vrts(moved_vrts_info)
        self.update_com_with_diff_info(moving_diff_info)
        
    def merge_coms_tentatively(self, part, com1, com2):
        moved_vrtset = self._memberset_in_com[part][com1]
        partnum = self._partnum
        moved_vrts_info = [{} for _ in range(partnum)]
        self._rollback_moved_vrts_info = [{} for _ in range(partnum)]
        for vrt in moved_vrtset:
            moved_vrts_info[part][vrt] = [com1, com2]
            self._rollback_moved_vrts_info[part][vrt] = [com2, com1]

        moving_diff_info = self.diff_of_moving_vrts(moved_vrts_info)
        self.update_com_with_diff_info(moving_diff_info)
        
    def move_vrts_tentatively(self, moved_vrts_info):
        partnum = self._partnum
        self._rollback_moved_vrts_info = [{} for _ in range(partnum)]
        for part, _moved_vrts_info in enumerate(moved_vrts_info):
            for vrt, (prev_com, next_com) in _moved_vrts_info.iteritems():
                self._rollback_moved_vrts_info[part][vrt] = [next_com, prev_com]

        moving_diff_info = self.diff_of_moving_vrts(moved_vrts_info)
        self.update_com_with_diff_info(moving_diff_info)

    def rollback_merging_coms(self):
        moving_diff_info = self.diff_of_moving_vrts(
                                            self._rollback_moved_vrts_info)
        self.update_com_with_diff_info(moving_diff_info)
        self._rollback_moved_vrts_info = None
    
    def rollback_moving_vrts(self):
        moving_diff_info = self.diff_of_moving_vrts(
                                            self._rollback_moved_vrts_info)
        self.update_com_with_diff_info(moving_diff_info)
        self._rollback_moved_vrts_info = None


    # status updater : whole community structure
    def _update_community_structure(self):
        edge_list = self._basic_info.get_edges_as_tuple()
        partnum = self._partnum

        egset_from_corres = {}
        egnum_from_corres = {}
        egset_from_com = [{} for _ in range(partnum)]
        egnum_from_com = [{} for _ in range(partnum)]
        memberset_in_com = [{} for _ in range(partnum)]
        membernum_in_com = [{} for _ in range(partnum)]
        corresset_from_com = [{} for _ in range(partnum)]

        for eg in edge_list:
            # correspondency
            corres = tuple( [self.com_of_vrt(part, vertex_index)
                             for part, vertex_index 
                             in enumerate(eg)] )
            egset_from_corres.setdefault(corres, set()).add(eg)
            egnum_from_corres[corres] = egnum_from_corres.get(corres, 0) + 1
            
            # community in each part
            for part, com in enumerate(corres):
                egset_from_com[part].setdefault(com, set()).add(eg)
                egnum_from_com[part][com] = egnum_from_com[part].get(com, 0) + 1
                memberset_in_com[part].setdefault(com, set())
                if eg[part] not in memberset_in_com[part][com]:
                    memberset_in_com[part][com].add(eg[part])
                    membernum_in_com[part][com] = \
                        membernum_in_com[part].get(com, 0) + 1
                corresset_from_com[part].setdefault(com, set()).add(corres)
                
        self._memberset_in_com = tuple(memberset_in_com)
        self._membernum_in_com = tuple(membernum_in_com)
        self._egset_from_com = tuple(egset_from_com)
        self._egnum_from_com = tuple(egnum_from_com)
        self._corresset_from_com = tuple(corresset_from_com)
        self._egset_from_corres = egset_from_corres
        self._egnum_from_corres = egnum_from_corres
        
    def update_com_with_diff_info(self, moving_diff_info):
        moved_vrts_info = moving_diff_info['moved_vrts_info']
        new_egset_from_com = moving_diff_info['new_egset_from_com']
        new_egnum_from_com = moving_diff_info['new_egnum_from_com']
        new_egset_from_corres = moving_diff_info['new_egset_from_corres']
        new_egnum_from_corres = moving_diff_info['new_egnum_from_corres']
        added_corresset = moving_diff_info['added_corresset']
 
        # vertex indexes, community member
        for part, _moved_vrt in enumerate(moved_vrts_info):
            for vrt, (prev_com, next_com) in _moved_vrt.iteritems():
                self._com_labels[part][vrt] = next_com
                self._memberset_in_com[part][prev_com].discard(vrt)
                self._memberset_in_com[part].setdefault(next_com, 
                                                          set()).add(vrt)
                self._membernum_in_com[part][prev_com] -= 1
                self._membernum_in_com[part].setdefault(next_com, 0)
                self._membernum_in_com[part][next_com] += 1

        # each community information
        for part, _new_egnum_from_com in enumerate(new_egnum_from_com):
            for com, egnum in _new_egnum_from_com.iteritems():
                if egnum == 0:
                    self._del_com(part, com)
                else:
                    _new_egset = new_egset_from_com[part][com]
                    self._egset_from_com[part][com] = _new_egset
                    self._egnum_from_com[part][com] = egnum
                
        # correspondence information
        for corres, egnum in new_egnum_from_corres.iteritems():
            if egnum == 0:
                self._del_corres(corres)
            else:
                _egset = new_egset_from_corres[corres]
                self._egset_from_corres[corres] = _egset
                self._egnum_from_corres[corres] = egnum
                
                if corres in added_corresset:
                    for part, com in enumerate(corres):
                        self._corresset_from_com[part].setdefault(com, set()).add(corres)

        # adjacent correspondency to community
        #for corres in added_corresset:
        #    for part, com in enumerate(corres):
        #        self._corresset_from_com[part].setdefault(com, set()).add(corres)

    # information for the tentative move
    def diff_of_merging_coms(self, part, com1, com2):
        moved_vrtset = self._memberset_in_com[part][com1]
        partnum = self._partnum
        moved_vrts_info = [{} for _ in range(partnum)]
        for vrt in moved_vrtset:
            moved_vrts_info[part][vrt] = [com1, com2]

        return self.diff_of_moving_vrts(moved_vrts_info)

    def diff_of_moving_vrts(self, moved_vrts_info):

        # communities and edges which are affected by the moves
        # moved_vrts_info, prev_com == next_com is the bad case.
        partnum = self._partnum
        moved_egset = set()
        affected_corresset = set()
        new_egset_from_com = [{} for _ in range(partnum)]
        new_egnum_from_com = [{} for _ in range(partnum)]
        #checked_coms = [set() for _ in range(partnum)]
        for part, _moved_vrt in enumerate(moved_vrts_info):
            for vrt, (prev_com, next_com) in _moved_vrt.iteritems():
                if prev_com not in new_egset_from_com[part]:
                    affected_corresset |= self.corresset_from_com(part, prev_com)
                if next_com not in new_egset_from_com[part]:
                    affected_corresset |= self.corresset_from_com(part, next_com)

                _moved_egset = self._basic_info.adj_row_egset_to_vrt(part, vrt)
                moved_egset |= _moved_egset

                new_egset_from_com[part].setdefault(
                    prev_com,
                    self.egset_from_com(part, prev_com).copy()
                    ).difference_update(_moved_egset)
                new_egset_from_com[part].setdefault(
                    next_com,
                    self.egset_from_com(part, next_com).copy()
                    ).update(_moved_egset)

                egnum = len(_moved_egset)
                egnum_from_com = new_egnum_from_com[part].setdefault(
                                    prev_com,
                                    self.egnum_from_com(part, prev_com)
                                    )
                new_egnum_from_com[part][prev_com] = egnum_from_com - egnum 
                egnum_from_com = new_egnum_from_com[part].setdefault(
                                    next_com,
                                    self.egnum_from_com(part, next_com)
                                    )
                new_egnum_from_com[part][next_com] = egnum_from_com + egnum

        # all correspondences which are affected by the move
        new_egnum_from_corres = {}
        new_egset_from_corres = {}
         
        added_corresset = set()
        for eg in moved_egset:
            prev_corres = [None] * partnum
            next_corres = [None] * partnum
            for part, vrt in enumerate(eg):
                _moved_vrt = moved_vrts_info[part]
                if vrt in _moved_vrt:
                    prev_com, next_com = _moved_vrt[vrt]
                    prev_corres[part] = prev_com
                    next_corres[part] = next_com
                else:
                    com = self.com_of_vrt(part, vrt)
                    prev_corres[part] = com
                    next_corres[part] = com

            prev_corres = tuple(prev_corres)
            if prev_corres not in new_egset_from_corres:
                new_egset_from_corres[prev_corres] = self.egset_from_corres(prev_corres).copy()
                new_egnum_from_corres[prev_corres] = self.egnum_from_corres(prev_corres)
 
            next_corres = tuple(next_corres)
            if next_corres not in new_egset_from_corres:
                added_corresset.add(next_corres)
                new_egset_from_corres[next_corres] = self.egset_from_corres(next_corres).copy()
                new_egnum_from_corres[next_corres] = self.egnum_from_corres(next_corres)

            new_egset_from_corres[prev_corres].remove(eg)
            new_egnum_from_corres[prev_corres] -= 1
            new_egset_from_corres[next_corres].add(eg)
            new_egnum_from_corres[next_corres] += 1

        for corres in affected_corresset:
            if corres not in new_egset_from_corres:
                new_egset_from_corres[corres] = self.egset_from_corres(corres).copy()
                new_egnum_from_corres[corres] = self.egnum_from_corres(corres)
 
        moving_diff_info = {
            'moved_vrts_info' : moved_vrts_info,
            'new_egset_from_com' : new_egset_from_com,
            'new_egnum_from_com' : new_egnum_from_com,
            'new_egset_from_corres' : new_egset_from_corres,
            'new_egnum_from_corres' : new_egnum_from_corres,
            'added_corresset' : added_corresset,
        }
        #print moving_diff_info
        #print affected_com
        #print moved_egset
        #print affected_corresset
        #print self._egset_from_corres
        #print self._corresset_from_com
        #exit()
        return moving_diff_info

    def minimal_diff_of_moving_vrts(self, moved_vrts_info):

        # communities and edges which are affected by the moves
        # moved_vrts_info, prev_com == next_com is the bad case.
        partnum = self._partnum
        moved_egset = set()
        affected_corresset = set()
        checked_com = [set() for _ in range(partnum)]
        new_egnum_from_com = [{} for _ in range(partnum)]
        #checked_coms = [set() for _ in range(partnum)]
        for part, _moved_vrt in enumerate(moved_vrts_info):
            for vrt, (prev_com, next_com) in _moved_vrt.iteritems():
                if prev_com not in checked_com[part]:
                    affected_corresset |= self.corresset_from_com(part, prev_com)
                if next_com not in checked_com[part]:
                    affected_corresset |= self.corresset_from_com(part, next_com)

                _moved_egset = self._basic_info.adj_row_egset_to_vrt(part, vrt)
                moved_egset |= _moved_egset

                #egnum = len(_moved_egset)
                egnum = self._basic_info.degree(part, vrt)
                egnum_from_com = new_egnum_from_com[part].setdefault(
                                    prev_com,
                                    self.egnum_from_com(part, prev_com)
                                    )
                new_egnum_from_com[part][prev_com] = egnum_from_com - egnum 
                egnum_from_com = new_egnum_from_com[part].setdefault(
                                    next_com,
                                    self.egnum_from_com(part, next_com)
                                    )
                new_egnum_from_com[part][next_com] = egnum_from_com + egnum

        # all correspondences which are affected by the move
        new_egnum_from_corres = {}
        for eg in moved_egset:
            prev_corres = [None] * partnum
            next_corres = [None] * partnum
            for part, vrt in enumerate(eg):
                _moved_vrt = moved_vrts_info[part]
                if vrt in _moved_vrt:
                    prev_com, next_com = _moved_vrt[vrt]
                    prev_corres[part] = prev_com
                    next_corres[part] = next_com
                else:
                    com = self.com_of_vrt(part, vrt)
                    prev_corres[part] = com
                    next_corres[part] = com

            prev_corres = tuple(prev_corres)
            if prev_corres not in new_egnum_from_corres:
                new_egnum_from_corres[prev_corres] = self.egnum_from_corres(prev_corres)
 
            next_corres = tuple(next_corres)
            if next_corres not in new_egnum_from_corres:
                new_egnum_from_corres[next_corres] = self.egnum_from_corres(next_corres)

            new_egnum_from_corres[prev_corres] -= 1
            new_egnum_from_corres[next_corres] += 1

        for corres in affected_corresset:
            if corres not in new_egnum_from_corres:
                new_egnum_from_corres[corres] = self.egnum_from_corres(corres)
 
        moving_diff_info = {
            #'moved_vrts_info' : moved_vrts_info,
            #'new_egset_from_com' : new_egset_from_com,
            'new_egnum_from_com' : new_egnum_from_com,
            #'new_egset_from_corres' : new_egset_from_corres,
            'new_egnum_from_corres' : new_egnum_from_corres,
            #'added_corresset' : added_corresset,
        }
        return moving_diff_info

class _EdgeStatus(object):

    def __init__(self, basic_info):
        self._basic_info = basic_info
        self._partnum = basic_info.partnum()
        self._egcl_labels = [None] * self._basic_info.edgenum()

        self._memberset_from_egcl = None
        self._membernum_from_egcl = None

    def get_egcl_labels(self):
        return self._egcl_labels

    def egcl_label(self, eg_ind):
        return self._egcl_labels[eg_ind]

    def egcl_size(self, egcl_label):
        return self._membernum_from_egcl[egcl_label]

    def memberset_from_egcl(self, egcl_label):
        return self._memberset_from_egcl[egcl_label]

    def iter_adj_egcl_to_vrt(self, part, vrt):
        adj_egclset = set()
        adj_egcls = [self.egcl_label(eg_ind)
                     for eg_ind 
                     in self._basic_info.adj_egset_to_vrt(part, vrt)]

        for adj_egcl in set(adj_egcls):
            yield adj_egcl

    def assign_unique_egcl_labels(self):
        egcl_labels = range( self._basic_info.edgenum() )
        self.set_egcl_labels(egcl_labels)
        return self._egcl_labels

    def set_egcl_labels(self, egcl_labels):
        self._egcl_labels = egcl_labels
        self._update_egcl_structure()

    def calculate_egclnum(self):
        return len(self._memberset_from_egcl)

    def iter_egcl_combination(self, combi_num=2):
        #iter_egcls = self._memberset_from_egcl.iterkeys
        #return itertools.combinations(iter_egcls(), combi_num)
        egcls = self._memberset_from_egcl.keys()
        return itertools.combinations(egcls, combi_num)

    def _update_egcl_structure(self):
        memberset_from_egcl = {}
        membernum_from_egcl = {}
        for eg_ind, egcl in enumerate(self._egcl_labels):
            memberset_from_egcl.setdefault(egcl, set()).add(eg_ind)
            membernum_from_egcl.setdefault(egcl, 0)
            membernum_from_egcl[egcl] += 1

        self._memberset_from_egcl = memberset_from_egcl
        self._membernum_from_egcl = membernum_from_egcl
        #print self._egcl_labels
        #print self._memberset_from_egcl
        #print self._membernum_from_egcl
        #print ''

    def merge_egcls(self, egcl1, egcl2):
        self.merge_egcls_tentatively(egcl1, egcl2)

    def merge_egcls_tentatively(self, egcl1, egcl2):
        assign_info = {}
        self._rollback_merge_info = None
        
        # delete egcl1
        egcl1_members = self._memberset_from_egcl[egcl1]
        egcl1_num = self._membernum_from_egcl[egcl1]
        del self._memberset_from_egcl[egcl1]
        del self._membernum_from_egcl[egcl1]

        # merge egcl1 to egcls
        self._memberset_from_egcl[egcl2].update(egcl1_members)
        self._membernum_from_egcl[egcl2] += egcl1_num

        # update labels of edges
        for eg_ind in egcl1_members:
            self._egcl_labels[eg_ind] = egcl2

        # store the information for rollback
        self._rollback_info_of_merging_egcls = \
            (egcl1, egcl1_members, egcl1_num, egcl2)

    def rollback_merging_egcls(self):
        egcl1, egcl1_members, egcl1_num, egcl2 = self._rollback_info_of_merging_egcls

        # split the merged egdge cluster
        self._memberset_from_egcl[egcl2].difference_update(egcl1_members)
        self._membernum_from_egcl[egcl2] -= egcl1_num
        self._memberset_from_egcl[egcl1] = egcl1_members
        self._membernum_from_egcl[egcl1] = egcl1_num

        # assign the labels to edges which are the members of splited cluster
        for eg_ind in egcl1_members:
            self._egcl_labels[eg_ind] = egcl1


class _HierarchicalEdgeStatus(_EdgeStatus):

    def __init__(self, basic_info):
        self._basic_info = basic_info
        self._partnum = basic_info.partnum()

        self._egcl_labels = [None] * self._basic_info.edgenum()

        self._memberset_from_egcl = None
        self._membernum_from_egcl = None
        self._adj_egclset_from_egcl = None
        self._adj_egclset_from_vrt = None
        self._egcl_label_from_egcl = None

        #self._egclset_from_egcl_label = None
        self._membernum_from_egcl_label = None
        self._memberset_from_egcl_label = None

    def egcl_label(self, egcl):
        return self._egcl_label_from_egcl[egcl]

    def egcl_size(self, egcl_label):
        return self._membernum_from_egcl_label[egcl_label]

    def memberset_from_egcl_label(self, egcl_label):
        return self._memberset_from_egcl_label[egcl_label]

    def iter_adj_egcl_to_vrt(self, part, vrt):
        return [self._egcl_label_from_egcl[adj_egcl] 
                    for adj_egcl
                    in self._adj_egclset_from_vrt[part][vrt] ]

    def egcls_randomly(self):
        egcls = self._membernum_from_egcl.keys()
        random.shuffle(egcls)
        return egcls

    def adj_egclset(self, egcl):
        return self._adj_egclset_from_egcl[egcl]

    def merge_egcls_hierarchically(self):
        self._update_egcl_structure()

    def calculate_egclnum(self):
        return len(self._egcls_from_egcl_label)

    def _update_egcl_structure(self):
        basic_info = self._basic_info
        partnum = self._partnum

        memberset_from_egcl = {}
        membernum_from_egcl = {}
        adj_egclset_from_egcl = {}
        adj_egclset_from_vrt = [{} for _ in xrange(partnum)]
        egcl_label_from_egcl = {}

        #egclset_from_egcl_label = {}
        memberset_from_egcl_label = {}
        membernum_from_egcl_label = {}

        for eg_ind, egcl in enumerate(self._egcl_labels):
            memberset_from_egcl.setdefault(egcl, set()).add(eg_ind)
            membernum_from_egcl.setdefault(egcl, 0)
            membernum_from_egcl[egcl] += 1

            for adj_eg_ind in basic_info.adj_egset_to_eg(eg_ind):
                adj_egcl = self._egcl_labels[adj_eg_ind]
                adj_egclset_from_egcl.setdefault(egcl, set()).add(adj_egcl)

            for part, vrt in enumerate(basic_info.get_edge(eg_ind)):
                for adj_eg_ind in basic_info.adj_egset_to_vrt(part, vrt):
                    adj_egcl = self._egcl_labels[adj_eg_ind]
                    adj_egclset_from_vrt[part].setdefault(
                        vrt, set() ).add(adj_egcl)

            egcl_label_from_egcl.setdefault(egcl, egcl)

            #egclset_from_egcl_label.setdefault(egcl, set([egcl]))
            memberset_from_egcl_label.setdefault(egcl, set()).add(eg_ind)
            membernum_from_egcl_label.setdefault(egcl, 0)
            membernum_from_egcl_label[egcl] += 1

        #for egcl_label, egclset in egclset_from_egcl_label.iteritems():
        #    membernum = 0
        #    memberset = set()
        #    for egcl in egclset:
        #        memberset.update( memberset_from_egcl[egcl] )
        #        membernum += membernum_from_egcl[egcl]
        #    membernum_from_egcl_label[egcl_label] = membernum
        #    memberset_from_egcl_label[egcl_label] = memberset

        self._memberset_from_egcl = memberset_from_egcl
        self._membernum_from_egcl = membernum_from_egcl
        self._adj_egclset_from_egcl = adj_egclset_from_egcl
        self._adj_egclset_from_vrt = adj_egclset_from_vrt
        self._egcl_label_from_egcl = egcl_label_from_egcl

        #self._egclset_from_egcl_label = egclset_from_egcl_label
        self._membernum_from_egcl_label = membernum_from_egcl_label
        self._memberset_from_egcl_label = memberset_from_egcl_label

        #print self._egcl_labels
        #print self._memberset_from_egcl
        #print self._membernum_from_egcl
        #print ''

    #def move_egcl(self, egcl, next_egcl_label):
    #    self.move_egcl_tentatively(egcl, next_egcl_label)

    def move_egcl_tentatively(self, egcl, next_egcl_label):
        memberset = self._memberset_from_egcl[egcl]
        membernum = self._membernum_from_egcl[egcl]

        for eg_ind in memberset:
            self._egcl_labels[eg_ind] = next_egcl_label

        prev_egcl_label = self._egcl_label_from_egcl[egcl]
        self._egcl_label_from_egcl[egcl] = next_egcl_label

        # update
        #self._egclset_from_egcl_label[prev_egcl_label].discard(egcl)
        #self._egclset_from_egcl_label[next_egcl_label].add(egcl)

        self._membernum_from_egcl_label[prev_egcl_label] -= membernum
        self._membernum_from_egcl_label[next_egcl_label] += membernum

        self._memberset_from_egcl_label[prev_egcl_label].difference_update(memberset)
        self._memberset_from_egcl_label[next_egcl_label].update(memberset)

        # store the information for rollback
        self._rollback_info_of_moving_egcl = \
            (egcl, prev_egcl_label, next_egcl_label)

    def rollback_moving_egcl(self):
        egcl, prev_egcl_label, next_egcl_label = \
            self._rollback_info_of_moving_egcl

        memberset = self._memberset_from_egcl[egcl]
        membernum = self._membernum_from_egcl[egcl]
        for eg_ind in memberset:
            self._egcl_labels[eg_ind] = prev_egcl_label

        self._egcl_label_from_egcl[egcl] = prev_egcl_label

        # update
        #self._egclset_from_egcl_label[prev_egcl_label].add(egcl)
        #self._egclset_from_egcl_label[next_egcl_label].remove(egcl)

        self._membernum_from_egcl_label[prev_egcl_label] += membernum
        self._membernum_from_egcl_label[next_egcl_label] -= membernum

        self._memberset_from_egcl_label[prev_egcl_label].update(memberset)
        self._memberset_from_egcl_label[next_egcl_label].difference_update(memberset)

        # the information for rollback
        self._rollback_info_of_moving_egcl = None
        
        
#--------------------------------------------------
# public classes
#--------------------------------------------------
class NetworkStatus(object):

    def __init__(self, edge_list):
        edge_list = self._transform_edges_to_tuple(edge_list)
        self.basic = _BasicInformation(edge_list)

    def _transform_edges_to_tuple(self, edges):
        if isinstance(edges[0], tuple):
            tuple_edges = tuple(edges)
        else:
            tuple_edges = tuple( map(tuple, edges) )

        return tuple_edges

    def add_com(self):
        self.com = _CommunityStatus(self.basic)

    def add_egcl(self, is_hierarchical=False):
        if is_hierarchical:
            self.egcl = _HierarchicalEdgeStatus(self.basic)
        else:
            self.egcl = _EdgeStatus(self.basic)


class ComEgclSynchronalManeger(object):

    def __init__(self, status):
        self._status = status
        self._partnum = status.basic.partnum()

    def assign_unique_egcl_labels(self):
        status = self._status
        # assign labels
        status.egcl.assign_unique_egcl_labels()

        # update communities
        self._assign_coms_based_on_all_egcls()

    def merge_egcls_tentatively(self, egcl1, egcl2):
        status = self._status
        # merge egcls
        status.egcl.merge_egcls_tentatively(egcl1, egcl2)

        # update communities
        target_egset = status.egcl.memberset_from_egcl(egcl2)
        self._assign_coms_based_on_egcls_tentatively(target_egset)
        
    def rollback_merging_egcls(self):
        self._status.egcl.rollback_merging_egcls()
        self._status.com.rollback_moving_vrts()

    def merge_egcls(self, egcl1, egcl2):
        self.merge_egcls_tentatively(egcl1, egcl2)

    def move_egcl_tentatively(self, egcl, next_egcl_label):
        status = self._status

        # move egcl
        prev_egcl_label = status.egcl.egcl_label(egcl)
        status.egcl.move_egcl_tentatively(egcl, next_egcl_label)

        # update communities
        target_egset = ( status.egcl.memberset_from_egcl_label(prev_egcl_label) |
                         status.egcl.memberset_from_egcl_label(next_egcl_label) )
        self._assign_coms_based_on_egcls_tentatively(target_egset)

    def rollback_moving_egcl(self):
        self._status.egcl.rollback_moving_egcl()
        self._status.com.rollback_moving_vrts()

    def move_egcl(self, egcl, next_egcl_label):
        self.move_egcl_tentatively(egcl, next_egcl_label)

    def _assign_coms_based_on_all_egcls(self):
        status = self._status
        com_labels = status.com.get_com_labels()
        partnum = self._partnum

        checked_vrtset = [set() for _ in xrange(partnum)]
        for eg in status.basic.get_edges_as_tuple():
            for part, vrt in enumerate(eg):
                if vrt in checked_vrtset[part]:
                    continue
                checked_vrtset[part].add(vrt)

                com = self._com_based_on_egcl(part, vrt)
                com_labels[part][vrt] = com

        status.com.set_com_labels(com_labels)

    def _com_based_on_egcl(self, part, vrt):
        egcl_status = self._status.egcl
        max_size = -1
        max_egcl = None
        for egcl in egcl_status.iter_adj_egcl_to_vrt(part, vrt):
        #for adj_egcl in egcl_status._adj_egclset_from_vrt[part][vrt]:
            #egcl = egcl_status.egcl_label(adj_egcl)
            size = egcl_status.egcl_size(egcl)
            #size = egcl_status._membernum_from_egcl_label[egcl]

            if size > max_size:
                max_size = size
                max_egcl = egcl
            elif size == max_size and max_egcl < egcl:
                max_egcl = egcl

        return max_egcl

    def _assign_coms_based_on_egcls_tentatively(self, target_egset):
        moved_vrts_info = self._moved_vrts_info_by_assignment(target_egset)

        # move the vertices tentatively
        self._status.com.move_vrts_tentatively(moved_vrts_info)

    def _moved_vrts_info_by_assignment(self, target_egset):
        status = self._status
        com_labels = status.com.get_com_labels()
        partnum = self._partnum

        # information of the moves of the vertices
        moved_vrts_info = [{} for _ in xrange(partnum)]
        checked_vrtset = [set() for _ in xrange(partnum)]
        for eg_ind in target_egset:
            eg = status.basic.get_edge(eg_ind)
            for part, vrt in enumerate(eg):
                # the vertex has already checked or not
                if vrt in checked_vrtset[part]:
                    continue
                checked_vrtset[part].add(vrt)

                next_com = self._com_based_on_egcl(part, vrt)
                prev_com = com_labels[part][vrt]
                if prev_com != next_com:
                    moved_vrts_info[part][vrt] = [prev_com, next_com]

        return moved_vrts_info
 
    def diff_of_moving_egcl(self, egcl, next_egcl_label):
        status = self._status

        # move egcl
        prev_egcl_label = status.egcl.egcl_label(egcl)
        status.egcl.move_egcl_tentatively(egcl, next_egcl_label)

        # moved vrtices information
        target_egset = ( status.egcl.memberset_from_egcl_label(prev_egcl_label) |
                         status.egcl.memberset_from_egcl_label(next_egcl_label) )
        moved_vrts_info = self._moved_vrts_info_by_assignment(target_egset)
        moving_diff_info = status.com.diff_of_moving_vrts(moved_vrts_info)

        # rollback
        status.egcl.rollback_moving_egcl()

        return moving_diff_info

    def minimal_diff_of_moving_egcl(self, egcl, next_egcl_label):
        status = self._status

        # move egcl
        prev_egcl_label = status.egcl.egcl_label(egcl)
        status.egcl.move_egcl_tentatively(egcl, next_egcl_label)

        # moved vrtices information
        target_egset = ( status.egcl.memberset_from_egcl_label(prev_egcl_label) |
                         status.egcl.memberset_from_egcl_label(next_egcl_label) )
        moved_vrts_info = self._moved_vrts_info_by_assignment(target_egset)
        moving_diff_info = status.com.minimal_diff_of_moving_vrts(moved_vrts_info)

        # rollback
        status.egcl.rollback_moving_egcl()

        return moving_diff_info

    '''
    def diff_of_moving_egcl_opt(self, egcl, prev_egcl_label, next_egcl_label):
        status = self._status
        partnum = self._partnum
        com_labels = status.com.get_com_labels()

        # update egcl size
        status.egcl._egcl_label_from_egcl[egcl] = next_egcl_label

        memberset = status.egcl._memberset_from_egcl[egcl]
        status.egcl._memberset_from_egcl_label[prev_egcl_label].difference_update(memberset)
        status.egcl._memberset_from_egcl_label[next_egcl_label].update(memberset)

        membernum = status.egcl._membernum_from_egcl[egcl]
        status.egcl._membernum_from_egcl_label[prev_egcl_label] -= membernum
        status.egcl._membernum_from_egcl_label[next_egcl_label] += membernum


        # moved vrts info
        target_egset = ( status.egcl.memberset_from_egcl_label(prev_egcl_label) |
                         status.egcl.memberset_from_egcl_label(next_egcl_label) )
        moved_vrts_info = [{} for _ in xrange(partnum)]
        checked_vrtset = [set() for _ in xrange(partnum)]
        for eg_ind in target_egset:
            eg = status.basic.get_edge(eg_ind)
            for part, vrt in enumerate(eg):
                # the vertex has already checked or not
                if vrt in checked_vrtset[part]:
                    continue
                checked_vrtset[part].add(vrt)

                # next com
                max_size = -1
                for _egcl in status.egcl.iter_adj_egcl_to_vrt(part, vrt):
                    size = status.egcl.egcl_size(_egcl)
                    if size > max_size:
                        max_size = size
                        max_egcl = _egcl
                next_com = max_egcl

                #next_com = self._com_based_on_egcl(part, vrt)
                prev_com = com_labels[part][vrt]
                if prev_com != next_com:
                    moved_vrts_info[part][vrt] = [prev_com, next_com]

        moving_diff_info = status.com.diff_of_moving_vrts(moved_vrts_info)

        # rollback
        status.egcl._egcl_label_from_egcl[egcl] = prev_egcl_label
        status.egcl._membernum_from_egcl_label[prev_egcl_label] += membernum
        status.egcl._membernum_from_egcl_label[next_egcl_label] -= membernum
        status.egcl._memberset_from_egcl_label[prev_egcl_label].update(memberset)
        status.egcl._memberset_from_egcl_label[next_egcl_label].difference_update(memberset)

        return moving_diff_info
    '''

if __name__ == 'main':
    pass
