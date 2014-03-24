#!/usr/bin/env python
# encoding: utf-8


#--------------------------------------------------
# modules
#--------------------------------------------------
from __future__ import division
import itertools
import copy
import random
from itertools import izip

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
        self._vertex_num_list = [len( set(vertex_list) ) 
                                 for vertex_list in zip(*edges)]
        self._total_vertex_num = sum(self._vertex_num_list)

        self._edge_num = len(edges)

        self._partnum = len(edges[0])
        self._partlist = range(self._partnum)

        self._adj_egset_from_vrt = \
            self._make_adjacency_between_vertex_and_edge(edges)

        self._degrees = [{} for _ in self._partlist]
        for part, adj_egset_from_vrt in enumerate(self._adj_egset_from_vrt):
            for vrt, adj_egset in adj_egset_from_vrt.iteritems():
                self._degrees[part][vrt] = len(adj_egset)

    def _make_adjacency_between_vertex_and_edge(self, edges):
        between_vrt_eg = [{} for _ in self._partlist]
        for eg_ind, eg in enumerate(edges):
            for part, vertex in izip(self._partlist, eg):
                between_vrt_eg[part].setdefault(vertex, set()).add(eg_ind)
            
        return tuple(between_vrt_eg)

    def edgenum(self):
        return self._edge_num

    def partnum(self):
        return self._partnum

    def partlist(self):
        return self._partlist

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
        self._basic = basic_info
        self._partnum = basic_info.partnum()
        self._partlist = basic_info.partlist()

        _vrtnum_list = self._basic.get_vertex_num_list()
        self._com_labels = [ [None] * _vrtnum 
                                  for _vrtnum in _vrtnum_list ]
 
        self._memberset_in_com = None
        self._membernum_in_com = None
        self._egset_from_com = None
        self._egnum_from_com = None

        self._egset_from_corres = None
        self._egnum_from_corres = None

    # whole community structures
    def count_all_comnums(self):
        return sum([len(_memberset_in_com) 
                    for _memberset_in_com in self._memberset_in_com])

    def com_labels(self):
        return self._com_labels

    # each community
    def com_of_vrt(self, part, vrt):
        return self._com_labels[part][vrt]

    def memberset_in_com(self, part, com):
        return self._memberset_in_com[part].get(com, set())
        #return self._memberset_in_com[part].setdefault(com, set())

    def membernum_in_com(self, part, com):
        return self._membernum_in_com[part].get(com, 0)
        #return self._membernum_in_com[part].setdefault(com, 0)

    def egset_from_com(self, part, com):
        return self._egset_from_com[part].get(com, set())
        #return self._egset_from_com[part].setdefault(com, set())

    def egnum_from_com(self, part, com):
        return self._egnum_from_com[part].get(com, 0)
        #return self._egnum_from_com[part].setdefault(com, 0)

    def iter_com_combinations(self, part, combi_num=2):
        coms = self._memberset_in_com[part].keys()
        return itertools.combinations(coms, combi_num)

    # correspondency
    def iter_corres_egnum(self):
        return self._egnum_from_corres.iteritems()

    def egset_from_corres(self, corres):
        return self._egset_from_corres.get(corres, set())
        #return self._egset_from_corres.setdefault(corres, set())

    def egnum_from_corres(self, corres):
        return self._egnum_from_corres.get(corres, 0)
        #return self._egnum_from_corres.setdefault(corres, 0)

    def _corres_from_eg(self, eg):
        corres = [self.com_of_vrt(part, vrt)
                  for part, vrt in izip(self._partlist, eg)]
        return tuple(corres)

    # community index
    def assign_unique_com_labels(self):
        vertex_num_list = self._basic.get_vertex_num_list()
        com_labels = [range(num) for num in vertex_num_list]
        self.set_com_labels(com_labels)

    def set_com_labels(self, com_labels):
        self._com_labels = com_labels
        self._update_community_structure()

    # used in this class
    def _del_com(self, part, com):
        del self._memberset_in_com[part][com]
        del self._membernum_in_com[part][com]
        del self._egset_from_com[part][com]
        del self._egnum_from_com[part][com]

    def _del_corres(self, corres):
        del self._egset_from_corres[corres]
        del self._egnum_from_corres[corres]

    # merge, move
    def merge_coms(self, part, com1, com2):
        partlist = self._partlist
        
        moved_vrtset = self._memberset_in_com[part][com1]
        moved_vrts_info = [{} for _ in partlist]
        for vrt in moved_vrtset:
            moved_vrts_info[part][vrt] = [com1, com2]

        moving_diff_info = self.diff_of_moving_vrts(moved_vrts_info)
        self.update_com_with_diff_info(moving_diff_info)
        
    def merge_coms_tentatively(self, part, com1, com2):
        partlist = self._partlist
        
        moved_vrtset = self._memberset_in_com[part][com1]
        moved_vrts_info = [{} for _ in partlist]
        self._rollback_moved_vrts_info = [{} for _ in partlist]

        for vrt in moved_vrtset:
            moved_vrts_info[part][vrt] = [com1, com2]
            self._rollback_moved_vrts_info[part][vrt] = [com2, com1]

        moving_diff_info = self.diff_of_moving_vrts(moved_vrts_info)
        self.update_com_with_diff_info(moving_diff_info)
        
    def move_vrts_tentatively(self, moved_vrts_info):

        self._rollback_moved_vrts_info = [
            { vrt : [next_com, prev_com] for 
                vrt, (prev_com, next_com) in _moved_vrts_info.iteritems() } 
            for _moved_vrts_info in moved_vrts_info ]

        moving_diff_info = self.diff_of_moving_vrts(moved_vrts_info)
        self.update_com_with_diff_info(moving_diff_info)

    def rollback_merging_coms(self):
        self.rollback_moving_vrts()

    def rollback_moving_vrts(self):
        moving_diff_info = self.diff_of_moving_vrts(
                                            self._rollback_moved_vrts_info)
        self.update_com_with_diff_info(moving_diff_info)
        self._rollback_moved_vrts_info = None


    # whole community structure
    def _update_community_structure(self):
        edge_list = self._basic.get_edges_as_tuple()
        partlist = self._partlist

        memberset_in_com = [{} for _ in partlist]
        membernum_in_com = [{} for _ in partlist]
        egset_from_com = [{} for _ in partlist]
        egnum_from_com = [{} for _ in partlist]

        egset_from_corres = {}
        egnum_from_corres = {}

        for eg in edge_list:
            # correspondency
            corres = self._corres_from_eg(eg)
            if corres not in egset_from_corres:
                egset_from_corres[corres] = set()
                egnum_from_corres[corres] = 0
            egset_from_corres[corres].add(eg)
            egnum_from_corres[corres] += 1
            
            # community in each part
            for part, vrt, com in izip(partlist, eg, corres):
                if com not in memberset_in_com[part]:
                    memberset_in_com[part][com] = set()
                    membernum_in_com[part][com] = 0
                    egset_from_com[part][com] = set()
                    egnum_from_com[part][com] = 0

                if vrt not in memberset_in_com[part][com]:
                    memberset_in_com[part][com].add(vrt)
                    membernum_in_com[part][com] += 1
                egset_from_com[part][com].add(eg)
                egnum_from_com[part][com] += 1
                
        self._memberset_in_com = tuple(memberset_in_com)
        self._membernum_in_com = tuple(membernum_in_com)
        self._egset_from_com = tuple(egset_from_com)
        self._egnum_from_com = tuple(egnum_from_com)

        self._egset_from_corres = egset_from_corres
        self._egnum_from_corres = egnum_from_corres

    def update_com_with_diff_info(self, moving_diff_info):
        partlist = self._partlist

        moved_vrts_info = moving_diff_info['moved_vrts_info']

        new_egset_from_com = moving_diff_info['new_egset_from_com']
        new_egnum_from_com = moving_diff_info['new_egnum_from_com']

        new_egset_from_corres = moving_diff_info['new_egset_from_corres']
        new_egnum_from_corres = moving_diff_info['new_egnum_from_corres']

        # community labels, members in community
        for part, _moved_vrt in izip(partlist, moved_vrts_info):
            for vrt, (prev_com, next_com) in _moved_vrt.iteritems():
                # labels
                self._com_labels[part][vrt] = next_com

                # prev com
                self._memberset_in_com[part][prev_com].remove(vrt)
                self._membernum_in_com[part][prev_com] -= 1

                # next com
                if next_com not in self._memberset_in_com[part]:
                    self._memberset_in_com[part][next_com] = set()
                    self._membernum_in_com[part][next_com] = 0
                self._memberset_in_com[part][next_com].add(vrt)
                self._membernum_in_com[part][next_com] += 1

        # edges from each community
        for part, _new_egnum_from_com in izip(partlist, new_egnum_from_com):
            for com, egnum in _new_egnum_from_com.iteritems():
                if egnum == 0:
                    self._del_com(part, com)
                else:
                    _new_egset = new_egset_from_com[part][com]
                    self._egset_from_com[part][com] = _new_egset
                    self._egnum_from_com[part][com] = egnum
                
        # correspondence
        for corres, egnum in new_egnum_from_corres.iteritems():
            if egnum == 0:
                self._del_corres(corres)
            else:
                _egset = new_egset_from_corres[corres]
                self._egset_from_corres[corres] = _egset
                self._egnum_from_corres[corres] = egnum
                
    def update_com_with_minimal_diff_info(self, moving_diff_info):
        partlist = self._partlist

        moved_vrts_info = moving_diff_info['moved_vrts_info']
        new_egnum_from_com = moving_diff_info['new_egnum_from_com']
        new_egnum_from_corres = moving_diff_info['new_egnum_from_corres']
        changed_corresset = moving_diff_info['changed_corresset']

        new_egset_from_com = [{} for _ in partlist]
        new_egset_from_corres = {}

        # community labels, members in community
        for part, _moved_vrt in izip(partlist, moved_vrts_info):
            for vrt, (prev_com, next_com) in _moved_vrt.iteritems():
                # labels
                self._com_labels[part][vrt] = next_com

                # prev com
                self._memberset_in_com[part][prev_com].remove(vrt)
                self._membernum_in_com[part][prev_com] -= 1

                # next com
                if next_com not in self._memberset_in_com[part]:
                    self._memberset_in_com[part][next_com] = set()
                    self._membernum_in_com[part][next_com] = 0
                self._memberset_in_com[part][next_com].add(vrt)
                self._membernum_in_com[part][next_com] += 1



                _adj_egset = self._basic.adj_row_egset_to_vrt(part, vrt)
                egnum = self._basic.degree(part, vrt)

                new_egset_from_com[part].setdefault(
                    prev_com, self.egset_from_com(part, prev_com)
                    ).difference_update(_adj_egset)

                new_egset_from_com[part].setdefault(
                    next_com, self.egset_from_com(part, next_com)
                    ).update(_adj_egset)

       # edges from each community
        for part, _new_egnum_from_com in izip(partlist, new_egnum_from_com):
            for com, egnum in _new_egnum_from_com.iteritems():
                if egnum == 0:
                    self._del_com(part, com)
                else:
                    _new_egset = new_egset_from_com[part][com]
                    self._egset_from_com[part][com] = _new_egset
                    self._egnum_from_com[part][com] = egnum


        for eg, (prev_corres, next_corres) in changed_corresset.iteritems():
            new_egset_from_corres.setdefault(
                prev_corres, self.egset_from_corres(prev_corres)
                ).remove(eg)
 
            new_egset_from_corres.setdefault(\
                next_corres, self.egset_from_corres(next_corres)
                ).add(eg)

        # correspondence
        for corres, egnum in new_egnum_from_corres.iteritems():
            if egnum == 0:
                self._del_corres(corres)
            else:
                _egset = new_egset_from_corres[corres]
                self._egset_from_corres[corres] = _egset
                self._egnum_from_corres[corres] = egnum
 
    # diff for merging and moving
    def diff_of_merging_coms(self, part, com1, com2):
        moved_vrtset = self._memberset_in_com[part][com1]
        moved_vrts_info = [{} for _ in self._partlist]
        for vrt in moved_vrtset:
            moved_vrts_info[part][vrt] = [com1, com2]

        return self.diff_of_moving_vrts(moved_vrts_info)

    def diff_of_moving_vrts(self, moved_vrts_info):
        partlist = self._partlist

        # edges from communities
        # BAD CASE : prev_com == next_com in moved_vrts_info
        affected_egset = set()
        new_egset_from_com = [{} for _ in partlist]
        new_egnum_from_com = [{} for _ in partlist]
        for part, _moved_vrt in izip(partlist, moved_vrts_info):
            for vrt, (prev_com, next_com) in _moved_vrt.iteritems():

                _adj_egset = self._basic.adj_row_egset_to_vrt(part, vrt)
                egnum = self._basic.degree(part, vrt)
                affected_egset |= _adj_egset

                if prev_com not in new_egset_from_com[part]:
                    new_egset_from_com[part][prev_com] = \
                        self.egset_from_com(part, prev_com).copy()
                    new_egnum_from_com[part][prev_com] = \
                        self.egnum_from_com(part, prev_com)

                new_egset_from_com[part][prev_com].difference_update(
                    _adj_egset)
                new_egnum_from_com[part][prev_com] -= egnum 

                if next_com not in new_egset_from_com[part]:
                    new_egset_from_com[part][next_com] = \
                        self.egset_from_com(part, next_com).copy()
                    new_egnum_from_com[part][next_com] = \
                        self.egnum_from_com(part, next_com)

                new_egset_from_com[part][next_com].update(_adj_egset)
                new_egnum_from_com[part][next_com] += egnum

        # correspondences which are affected by the move
        new_egnum_from_corres = {}
        new_egset_from_corres = {}
        #added_corresset = set()
        for eg in affected_egset:
            prev_corres = self._corres_from_eg(eg)
            next_corres = [_moved_vrt.get(
                            vrt, [None, self.com_of_vrt(part, vrt)])[1] for
                           part, vrt, _moved_vrt in
                           izip(partlist, eg, moved_vrts_info)]
            next_corres = tuple(next_corres)

            if prev_corres not in new_egset_from_corres:
                new_egset_from_corres[prev_corres] = \
                    self.egset_from_corres(prev_corres).copy()
                new_egnum_from_corres[prev_corres] = \
                    self.egnum_from_corres(prev_corres)
 
            new_egset_from_corres[prev_corres].remove(eg)
            new_egnum_from_corres[prev_corres] -= 1

            if next_corres not in new_egset_from_corres:
                #added_corresset.add(next_corres)
                new_egset_from_corres[next_corres] = \
                    self.egset_from_corres(next_corres).copy()
                new_egnum_from_corres[next_corres] = \
                    self.egnum_from_corres(next_corres)

            new_egset_from_corres[next_corres].add(eg)
            new_egnum_from_corres[next_corres] += 1

        moving_diff_info = {
            'moved_vrts_info' : moved_vrts_info,
            'new_egset_from_com' : new_egset_from_com,
            'new_egnum_from_com' : new_egnum_from_com,
            'new_egset_from_corres' : new_egset_from_corres,
            'new_egnum_from_corres' : new_egnum_from_corres,
        }
        return moving_diff_info

    def minimal_diff_of_moving_vrts(self, moved_vrts_info):
        partlist = self._partlist

        # edges from communities
        # BAD CASE : prev_com == next_com in moved_vrts_info
        affected_egset = set()
        new_egnum_from_com = [{} for _ in partlist]
        for part, _moved_vrt in izip(partlist, moved_vrts_info):
            for vrt, (prev_com, next_com) in _moved_vrt.iteritems():

                _adj_egset = self._basic.adj_row_egset_to_vrt(part, vrt)
                egnum = self._basic.degree(part, vrt)
                affected_egset |= _adj_egset

                new_egnum_from_com[part].setdefault(
                    prev_com, self.egnum_from_com(part, prev_com))
                new_egnum_from_com[part][prev_com] -= egnum 

                new_egnum_from_com[part].setdefault(
                    next_com, self.egnum_from_com(part, next_com))
                new_egnum_from_com[part][next_com] += egnum

        # correspondences which are affected by the move
        new_egnum_from_corres = {}
        changed_corresset = {}
        for eg in affected_egset:
            prev_corres = self._corres_from_eg(eg)
            next_corres = [_moved_vrt.get(
                            vrt, [None, self.com_of_vrt(part, vrt)])[1] for
                           part, vrt, _moved_vrt in
                           izip(partlist, eg, moved_vrts_info)]
            next_corres = tuple(next_corres)
            changed_corresset[eg] = (prev_corres, next_corres)

            new_egnum_from_corres.setdefault(
                prev_corres, self.egnum_from_corres(prev_corres))
            new_egnum_from_corres[prev_corres] -= 1

            new_egnum_from_corres.setdefault(
                next_corres, self.egnum_from_corres(next_corres))
            new_egnum_from_corres[next_corres] += 1

        moving_diff_info = {
            'moved_vrts_info' : moved_vrts_info,
            'new_egnum_from_com' : new_egnum_from_com,
            'new_egnum_from_corres' : new_egnum_from_corres,
            'changed_corresset' : changed_corresset,
        }
        return moving_diff_info

class _EdgeStatus(object):

    def __init__(self, basic_info):
        self._basic = basic_info
        self._partnum = basic_info.partnum()
        self._edgenum = basic_info.edgenum()

        self._memberset_from_egcl = None
        self._membernum_from_egcl = None
        self._egcl_labels = None

    def egcl_size(self, egcl_label):
        return self._membernum_from_egcl[egcl_label]

    def memberset_from_egcl(self, egcl_label):
        return self._memberset_from_egcl[egcl_label]

    def adj_egclset_to_vrt(self, part, vrt):
        adj_egcls = [self._egcl_labels[eg_ind]
                     for eg_ind 
                     in self._basic.adj_egset_to_vrt(part, vrt)]

        return set(adj_egcls)

    def assign_unique_egcl_labels(self):
        egcl_labels = range(self._edgenum)
        self._set_egcl_labels(egcl_labels)

    def _set_egcl_labels(self, egcl_labels):
        self._egcl_labels = egcl_labels
        self._update_egcl_structure()

    def count_egclnum(self):
        return len(self._memberset_from_egcl)

    def iter_egcl_combination(self, combi_num=2):
        #iter_egcls = self._memberset_from_egcl.iterkeys
        #return itertools.combinations(iter_egcls(), combi_num)
        egcls = self._memberset_from_egcl.keys()
        return itertools.combinations(egcls, combi_num)

    def _update_egcl_structure(self):
        memberset_from_egcl = {}
        membernum_from_egcl = {}
        for eg_ind, egcl in izip(xrange(self._edgenum), self._egcl_labels):
            memberset_from_egcl.setdefault(egcl, set()).add(eg_ind)
            membernum_from_egcl.setdefault(egcl, 0)
            membernum_from_egcl[egcl] += 1

        self._memberset_from_egcl = memberset_from_egcl
        self._membernum_from_egcl = membernum_from_egcl

    def merge_egcls(self, egcl1, egcl2):
        self.merge_egcls_tentatively(egcl1, egcl2)

    def merge_egcls_tentatively(self, egcl1, egcl2):

        # delete egcl1
        egcl1_members = self._memberset_from_egcl[egcl1]
        egcl1_num = self._membernum_from_egcl[egcl1]
        del self._memberset_from_egcl[egcl1]
        del self._membernum_from_egcl[egcl1]

        # merge egcl1 to egcl2
        self._memberset_from_egcl[egcl2].update(egcl1_members)
        self._membernum_from_egcl[egcl2] += egcl1_num

        # update labels of edges
        for eg_ind in egcl1_members:
            self._egcl_labels[eg_ind] = egcl2

        # store the information for rollback
        self._rollback_info_of_merging_egcls = \
            (egcl1, egcl1_members, egcl1_num, egcl2)

    def rollback_merging_egcls(self):
        egcl1, egcl1_members, egcl1_num, egcl2 = \
            self._rollback_info_of_merging_egcls

        # split the merged egdge cluster
        self._memberset_from_egcl[egcl2].difference_update(egcl1_members)
        self._membernum_from_egcl[egcl2] -= egcl1_num
        self._memberset_from_egcl[egcl1] = egcl1_members
        self._membernum_from_egcl[egcl1] = egcl1_num

        # assign the labels to edges which are the members of splited cluster
        for eg_ind in egcl1_members:
            self._egcl_labels[eg_ind] = egcl1

class _HierarchicalEdgeStatus(object):

    def __init__(self, basic_info):
        self._basic = basic_info
        self._partnum = basic_info.partnum()
        self._partlist = basic_info.partlist()
        self._edgenum = basic_info.edgenum()

        self._adj_egclset_to_vrt = None
        self._egcl_of_edge = [None] * self._edgenum

        self._label_of_egcl = None
        self._egset_of_egcl = None
        self._egnum_of_egcl = None
        self._adj_egclset_to_egcl = None

        self._egclset_of_label = None
        self._egnum_of_label = None

    def adj_egclset_to_vrt(self, part, vrt):
        return self._adj_egclset_to_vrt[part][vrt]

    def egcls_randomly(self):
        egcls = self._egnum_of_egcl.keys()
        random.shuffle(egcls)
        return egcls

    def label_of_egcl(self, egcl):
        return self._label_of_egcl[egcl]

    def egset_of_egcl(self, egcl):
        return self._egset_of_egcl[egcl]

    def adj_egclset_to_egcl(self, egcl):
        return self._adj_egclset_to_egcl[egcl]

    def egnum_of_label(self, egcl_label):
        return self._egnum_of_label[egcl_label]

    def egclset_of_label(self, egcl_label):
        return self._egclset_of_label[egcl_label]

    def assign_unique_egcl_labels(self):
        egcls_of_edge = range(self._edgenum)
        self._set_egcl_labels(egcls_of_edge)

    def _set_egcl_labels(self, egcls_of_edge):
        self._egcl_of_edge = egcls_of_edge
        self._update_egcl_structure()

    def merge_egcls_hierarchically(self):
        self._update_egcl_structure()

    def _update_egcl_structure(self):
        basic = self._basic
        partlist = self._partlist

        egset_of_egcl = {}
        egnum_of_egcl = {}
        label_of_egcl = {}
        adj_egclset_to_egcl = {}

        adj_egclset_to_vrt = [{} for _ in partlist]

        egclset_of_label = {}
        egnum_of_label = {}

        for eg_ind, egcl in enumerate(self._egcl_of_edge):
            egset_of_egcl.setdefault(egcl, set()).add(eg_ind)
            egnum_of_egcl.setdefault(egcl, 0)
            egnum_of_egcl[egcl] += 1

            # adj_egclset_to_egcl
            for adj_eg_ind in basic.adj_egset_to_eg(eg_ind):
                adj_egcl = self._egcl_of_edge[adj_eg_ind]
                adj_egclset_to_egcl.setdefault(egcl, set()).add(adj_egcl)

            # adj_egclset_to_vrt
            for part, vrt in izip(partlist, basic.get_edge(eg_ind)):
                adj_egclset_to_vrt[part].setdefault(vrt, set()).add(egcl)

            # next_egcl
            label_of_egcl[egcl] = egcl
            egclset_of_label.setdefault(egcl, set()).add(egcl)
            egnum_of_label.setdefault(egcl, 0)
            egnum_of_label[egcl] += 1

        self._egset_of_egcl = egset_of_egcl
        self._egnum_of_egcl = egnum_of_egcl
        self._label_of_egcl = label_of_egcl
        self._adj_egclset_to_egcl = adj_egclset_to_egcl

        self._adj_egclset_to_vrt = adj_egclset_to_vrt

        self._egclset_of_label = egclset_of_label
        self._egnum_of_label = egnum_of_label

    def move_egcl(self, egcl, next_label):
        self.move_egcl_tentatively(egcl, next_label)

        # edges
        egset = self._egset_of_egcl[egcl]
        for eg_ind in egset:
            self._egcl_of_edge[eg_ind] = next_label

    def move_egcl_tentatively(self, egcl, next_label):
        egnum = self._egnum_of_egcl[egcl]

        # egcls
        prev_label = self._label_of_egcl[egcl]
        self._label_of_egcl[egcl] = next_label

        self._egclset_of_label[prev_label].remove(egcl)
        self._egclset_of_label[next_label].add(egcl)

        # egnum
        self._egnum_of_label[prev_label] -= egnum
        self._egnum_of_label[next_label] += egnum

        # store the information for rollback
        self._rollback_info_of_moving_egcl = \
            (egcl, prev_label, next_label)

    def rollback_moving_egcl(self):
        egcl, prev_label, next_label = \
            self._rollback_info_of_moving_egcl

        egnum = self._egnum_of_egcl[egcl]

        # egcls
        self._label_of_egcl[egcl] = prev_label

        self._egclset_of_label[prev_label].add(egcl)
        self._egclset_of_label[next_label].remove(egcl)

        # egnum
        self._egnum_of_label[prev_label] += egnum
        self._egnum_of_label[next_label] -= egnum

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

    def add_egcl(self):
        self.egcl = _EdgeStatus(self.basic)

    def add_hiegcl(self):
        self.hiegcl = _HierarchicalEdgeStatus(self.basic)

class ComEgclSynchronalManeger(object):

    def __init__(self, status):
        self._status = status
        self._partnum = status.basic.partnum()
        self._partlist = status.basic.partlist()

    def assign_unique_egcl_labels(self):
        status = self._status
        # assign labels
        status.egcl.assign_unique_egcl_labels()

        # update communities
        self._assign_coms_based_on_all_egcls()

    def merge_egcls(self, egcl1, egcl2):
        self.merge_egcls_tentatively(egcl1, egcl2)

    def merge_egcls_tentatively(self, egcl1, egcl2):
        status = self._status
        # merge egcls
        status.egcl.merge_egcls_tentatively(egcl1, egcl2)

        # update communities
        target_egset = status.egcl.memberset_from_egcl(egcl2)
        moved_vrts_info = self._moved_vrts_info_by_assignment(target_egset)
        self._status.com.move_vrts_tentatively(moved_vrts_info)
        
    def rollback_merging_egcls(self):
        self._status.egcl.rollback_merging_egcls()
        self._status.com.rollback_moving_vrts()

    def _assign_coms_based_on_all_egcls(self):
        status = self._status
        partlist = self._partlist
        com_labels = status.com.com_labels()

        checked_vrtset = [set() for _ in partlist]
        for eg in status.basic.get_edges_as_tuple():
            for part, vrt in izip(partlist, eg):
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
        for egcl in egcl_status.adj_egclset_to_vrt(part, vrt):
            size = egcl_status.egcl_size(egcl)

            if size > max_size:
                max_size = size
                max_egcl = egcl
            elif size == max_size and max_egcl < egcl:
                max_egcl = egcl

        return max_egcl

    def _moved_vrts_info_by_assignment(self, target_egset):
        status = self._status
        com_labels = status.com.com_labels()
        partlist = self._partlist

        # information of the moves of the vertices
        moved_vrts_info = [{} for _ in partlist]
        checked_vrtset = [set() for _ in partlist]
        for eg_ind in target_egset:
            eg = status.basic.get_edge(eg_ind)
            for part, vrt in izip(partlist, eg):
                # the vertex has already checked or not
                if vrt in checked_vrtset[part]:
                    continue
                checked_vrtset[part].add(vrt)

                next_com = self._com_based_on_egcl(part, vrt)
                prev_com = com_labels[part][vrt]
                if prev_com != next_com:
                    moved_vrts_info[part][vrt] = [prev_com, next_com]

        return moved_vrts_info
 
class ComHiegclSynchronalManeger(object):

    def __init__(self, status):
        self._status = status
        self._partnum = status.basic.partnum()
        self._partlist = status.basic.partlist()

    def assign_unique_egcl_labels(self):
        status = self._status

        # assign labels to edges
        status.hiegcl.assign_unique_egcl_labels()

        # update communities
        self._assign_coms_based_on_all_egcls()

    def move_egcl(self, egcl, next_label):
        self.move_egcl_tentatively(egcl, next_label)

    def move_egcl_tentatively(self, egcl, next_label):
        status = self._status

        # move egcl
        prev_label = status.hiegcl.label_of_egcl(egcl)
        status.hiegcl.move_egcl_tentatively(egcl, next_label)

        # update communities
        target_egclset = ( status.hiegcl.egclset_of_label(prev_label) |
                           status.hiegcl.egclset_of_label(next_label) )
        target_egset = set()
        for egcl in target_egclset:
            target_egset |= status.hiegcl.egset_of_egcl(egcl)
        moved_vrts_info = self._moved_vrts_info_by_assignment(target_egset)
        self._status.com.move_vrts_tentatively(moved_vrts_info)

    def rollback_moving_egcl(self):
        self._status.hiegcl.rollback_moving_egcl()
        self._status.com.rollback_moving_vrts()

    def _assign_coms_based_on_all_egcls(self):
        status = self._status
        com_labels = status.com.com_labels()
        partnum = self._partnum
        partlist = self._partlist

        checked_vrtset = [set() for _ in partlist]
        for eg in status.basic.get_edges_as_tuple():
            for part, vrt in enumerate(eg):
                if vrt in checked_vrtset[part]:
                    continue
                checked_vrtset[part].add(vrt)

                com = self._com_based_on_egcl(part, vrt)
                com_labels[part][vrt] = com

        status.com.set_com_labels(com_labels)

    def _com_based_on_egcl(self, part, vrt):
        hiegcl = self._status.hiegcl
        max_size = -1
        max_label = None
        for egcl in hiegcl.adj_egclset_to_vrt(part, vrt):
            egcl_label = hiegcl.label_of_egcl(egcl)
            size = hiegcl.egnum_of_label(egcl_label)

            if size > max_size:
                max_size = size
                max_label = egcl_label
            elif size == max_size and max_label < egcl_label:
                max_label = egcl_label

        return max_label

    def _moved_vrts_info_by_assignment(self, target_egset):
        status = self._status
        com_labels = status.com.com_labels()
        com_based_on_egcl = self._com_based_on_egcl
        partnum = self._partnum
        partlist = self._partlist

        # information of the moves of the vertices
        moved_vrts_info = [{} for _ in partlist]
        checked_vrtset = [set() for _ in partlist]
        for eg_ind in target_egset:
            eg = status.basic.get_edge(eg_ind)
            for part, vrt in izip(partlist, eg):
                # the vertex has already checked or not
                if vrt in checked_vrtset[part]:
                    continue
                checked_vrtset[part].add(vrt)

                next_com = com_based_on_egcl(part, vrt)
                prev_com = com_labels[part][vrt]
                if prev_com != next_com:
                    moved_vrts_info[part][vrt] = [prev_com, next_com]

        return moved_vrts_info
 
    def diff_of_moving_egcl(self, egcl, next_label):
        status = self._status

        # move egcl
        prev_label = status.hiegcl.label_of_egcl(egcl)
        status.hiegcl.move_egcl_tentatively(egcl, next_label)

        # moved vrtices information
        target_egclset = ( status.hiegcl.egclset_of_label(prev_label) |
                           status.hiegcl.egclset_of_label(next_label) )
        target_egset = set()
        for egcl in target_egclset:
            target_egset |= status.hiegcl.egset_of_egcl(egcl)
        moved_vrts_info = self._moved_vrts_info_by_assignment(target_egset)
        #moving_diff_info = status.com.diff_of_moving_vrts(moved_vrts_info)
        moving_diff_info = \
            status.com.minimal_diff_of_moving_vrts(moved_vrts_info)

        # rollback
        status.hiegcl.rollback_moving_egcl()

        return moving_diff_info

if __name__ == 'main':
    pass

