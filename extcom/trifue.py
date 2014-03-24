#!/usr/bin/env python
# encoding: utf-8


#--------------------------------------------------
# modules
#--------------------------------------------------
from __future__ import division

from random import Random
import copy

import modularity as newmod

class _WrappingOptimization(object):
    
    def __init__(self):
        pass
        
    def name(self):
        return self.__class__.__name__
    
    def start(self, modularity, edge_list):
        
        if isinstance(edge_list[0], tuple):
            self.edges = edge_list
        else:
            self.edges = tuple([tuple(e) for e in edge_list])
        xyz = zip(*edge_list)
        self.num_x = max(xyz[0]) + 1
        self.num_y = max(xyz[1]) + 1
        self.num_z = max(xyz[2]) + 1
        
        if isinstance(modularity, newmod.NeubauerModularity):
            self.mod_proto = NeubauerTripartiteModularityOld(self.edges)
            
        ans_mod, mod_list, ans_ind = self.optimize()
        return ans_mod, ans_ind, mod_list

class TriFastUnfoldingForEdges(_WrappingOptimization):
    
    def optimize(self):
        
        # make some dictionaries for convenience
        eg_set_from_egcl, adj_egcl_set_from_egcl, adj_egcl_set_from_vrt_ind = \
            self._make_dictionaries_initially()
            
        # make community indexes of vertices
        vrt_ind = self._make_vrt_ind_initially()
        
        # current modularity
        mod_val = self.mod_proto.get_modularity(vrt_ind)
        mod_list = []
        
        # optimize modularity hierarchically
        while True:
            # optimize modularity with current hierarchical structure 
            new_mod_val, new_mod_list, new_vrt_ind, egcl_label_from_egcl = \
                self._optimize_modularity(mod_val,
                                          eg_set_from_egcl,
                                          adj_egcl_set_from_egcl,
                                          adj_egcl_set_from_vrt_ind,
                                          vrt_ind)
            
            # check the increase of modularity
            if new_mod_val <= mod_val:
                ans_ind = vrt_ind
                ans_mod = mod_val
                break
            
            # revise max value of modularity and community indexes  
            mod_val = new_mod_val
            mod_list.append(new_mod_list)
            vrt_ind = new_vrt_ind
            
            # merge edge clusters
            eg_set_from_egcl = \
                self._make_egcl_and_eg_member(eg_set_from_egcl, 
                                              egcl_label_from_egcl)
            # revise data structures
            adj_egcl_set_from_egcl = \
                self._make_ajacency_between_edge_cluster(eg_set_from_egcl)
            adj_egcl_set_from_vrt_ind = \
                self._make_adjacency_between_vertex_edge_cluster(
                                                         eg_set_from_egcl)
            
        ans_ind = ( self._adjast_index(ans_ind[0]),
                    self._adjast_index(ans_ind[1]),
                    self._adjast_index(ans_ind[2]) )
        
        return (ans_mod, mod_list, ans_ind)
    
    def _make_dictionaries_initially(self):
        
        # key   : index of an edges cluster
        # value : set of members of the cluster
        eg_set_from_egcl = dict( [ (ind, set([e]) )
                                 for ind, e in enumerate(self.edges)] )
        
        # key   : index of an edges cluster
        # value : set of indexes of edge clusters which are 
        #         adjacent to the edge cluster of the key
        adj_egcl_set_from_egcl = \
            self._make_ajacency_between_edge_cluster(eg_set_from_egcl)
        
        # key   : index of a vertex
        # value : set of indexes of edge clusters which are 
        #         adjacent to the vertex of the key
        adj_egcl_set_from_vrt_ind = \
            self._make_adjacency_between_vertex_edge_cluster(eg_set_from_egcl)
        
        return eg_set_from_egcl, adj_egcl_set_from_egcl, adj_egcl_set_from_vrt_ind
    
    def _make_ajacency_between_edge_cluster(self, eg_set_from_egcl):
        adj_egcl_set_from_egcl = {}
        
        between_vrt_egcl = self._make_adjacency_between_vertex_edge_cluster(eg_set_from_egcl)
        adj_egcl_set_from_vrt_x = between_vrt_egcl['x']
        adj_egcl_set_from_vrt_y = between_vrt_egcl['y']
        adj_egcl_set_from_vrt_z = between_vrt_egcl['z']
        
        for cls_ind, eg_list in eg_set_from_egcl.iteritems():
            for (_x, _y, _z) in eg_list:
                adj_egcl_set_from_egcl[cls_ind] = (adj_egcl_set_from_vrt_x[_x] | 
                                                   adj_egcl_set_from_vrt_y[_y] | 
                                                   adj_egcl_set_from_vrt_z[_z])
            adj_egcl_set_from_egcl[cls_ind].discard(cls_ind)
        
        return adj_egcl_set_from_egcl
    
    def _make_adjacency_between_vertex_edge_cluster(self, eg_set_from_egcl):
        between_vrt_egcl = { 'x' : {},
                             'y' : {},
                             'z' : {} }
        for egcl_ind, eg_set in eg_set_from_egcl.iteritems():
            for (_x, _y, _z) in eg_set:
                between_vrt_egcl['x'].setdefault(_x, set()).add(egcl_ind)
                between_vrt_egcl['y'].setdefault(_y, set()).add(egcl_ind)
                between_vrt_egcl['z'].setdefault(_z, set()).add(egcl_ind)
            
        return between_vrt_egcl
    
    def _make_vrt_ind_initially(self):
        cur_ind_x = [-1] * self.num_x
        cur_ind_y = [-1] * self.num_y
        cur_ind_z = [-1] * self.num_z
        for egcl_ind, (_x, _y, _z) in enumerate(self.edges):
            if cur_ind_x[_x] < egcl_ind:
                cur_ind_x[_x] = egcl_ind
            if cur_ind_y[_y] < egcl_ind:
                cur_ind_y[_y] = egcl_ind
            if cur_ind_z[_z] < egcl_ind:
                cur_ind_z[_z] = egcl_ind
                
        ind = (cur_ind_x, cur_ind_y, cur_ind_z)
        return ind
    
    def _make_egcl_and_eg_member(self, eg_set_from_egcl, egcl_label_from_egcl):
        _eg_set_from_egcl = {}
        for egcl_ind, new_egcl_ind in egcl_label_from_egcl.iteritems():
            eg_set = eg_set_from_egcl[egcl_ind]
            _eg_set_from_egcl.setdefault(new_egcl_ind, set()).update(eg_set)
            
        return _eg_set_from_egcl
    
    def _adjast_index(self, index):
        dic = {}
        num = 0
        for ind in index:
            if ind not in dic.keys():
                dic.setdefault(ind, num)
                num += 1
                
        clean_ind = index[:]
        for ith, ind in enumerate(index):
            clean_ind[ith] = dic[ind]
            
        return clean_ind

    def _optimize_modularity(self, mod_val, eg_set_from_egcl, 
                             adj_egcl_set_from_egcl, adj_egcl_set_from_vrt_ind, 
                             vrt_ind):
        
        # prepare for calculation
        random = Random() # used for shuffling the edge clusters
        calculate_delta = self.mod_proto.calculate_delta_value # to refer to the function directly
        reflect_the_moves = self.mod_proto.reflect_the_moves # same reason as the above
        get_vrt_ind = self._get_vrt_ind # same reason as the above
        x_vrt_ind, y_vrt_ind, z_vrt_ind = vrt_ind
        mod_list = [mod_val]
        
        # make some data structures for convenience
        egcl_ind_list = eg_set_from_egcl.keys()
        egcl_label_from_egcl = dict( [ (egcl_ind, egcl_ind) 
                                    for egcl_ind in egcl_ind_list] )
        
        egcl_set_from_egcl_label = dict( [ (egcl_ind, set([egcl_ind])) 
                                        for egcl_ind in egcl_ind_list] )
        
        egcl_size_from_egcl = dict( [ (egcl_ind, len(eg_set)) 
                                     for (egcl_ind, eg_set) 
                                     in eg_set_from_egcl.iteritems()] )
        
        egcl_size_from_egcl_label = dict( [ (egcl_ind, len(eg_set)) 
                                           for (egcl_ind, eg_set) 
                                           in eg_set_from_egcl.iteritems()] )
        
        # optimize modularity with some moves
        prev_mod_val = -1
        current_mod_val = mod_val
        while prev_mod_val < current_mod_val:
            prev_mod_val = current_mod_val
            
            # take an edge cluster randomly
            random.shuffle(egcl_ind_list)
            for egcl_ind in egcl_ind_list:
                egcl_label = egcl_label_from_egcl[egcl_ind]
                egcl_size = egcl_size_from_egcl[egcl_ind]
                adj_egcl_set = adj_egcl_set_from_egcl[egcl_ind]
                
                # move the edge cluster to its adjacent cluster temporary
                max_delta_mod = -1  # the highest modularity value derived by the move
                for adj_egcl_ind in adj_egcl_set:
                    adj_egcl_label = egcl_label_from_egcl[adj_egcl_ind]
                    if egcl_label == adj_egcl_label:
                        continue
                    
                    # move tentatively 
                    egcl_label_from_egcl[egcl_ind] = adj_egcl_label
                    egcl_set_from_egcl_label[egcl_label].discard(egcl_ind)
                    egcl_size_from_egcl_label[egcl_label] -= egcl_size
                    egcl_set_from_egcl_label[adj_egcl_label].add(egcl_ind)
                    egcl_size_from_egcl_label[adj_egcl_label] += egcl_size
                    
                    # assign vertices label based on edge cluster
                    _ind, _moved_vrt = get_vrt_ind(
                                                    eg_set_from_egcl,
                                                    egcl_label_from_egcl,
                                                    egcl_set_from_egcl_label,
                                                    egcl_size_from_egcl_label,
                                                    adj_egcl_set_from_vrt_ind,
                                                    egcl_label,
                                                    adj_egcl_label,
                                                    vrt_ind)
                    
                    # calculate modularity
                    #_mod_val = self.mod_proto.get_modularity(_ind) # calculate modularity naively
                    _delta_mod, _reflection_info = calculate_delta(_moved_vrt) # calculate the differences
                    
                    # check the increase of modularity
                    if _delta_mod > max_delta_mod:
                        max_delta_mod = _delta_mod
                        max_adj_egcl_label = adj_egcl_label
                        max_moved_vrt = _moved_vrt
                        max_reflection_info = _reflection_info
                    
                    # restore the states before the move
                    egcl_label_from_egcl[egcl_ind] = egcl_label
                    egcl_set_from_egcl_label[egcl_label].add(egcl_ind)
                    egcl_size_from_egcl_label[egcl_label] += egcl_size
                    egcl_set_from_egcl_label[adj_egcl_label].discard(egcl_ind)
                    egcl_size_from_egcl_label[adj_egcl_label] -= egcl_size
                    for _ind, (_prev, _next) in _moved_vrt['x'].iteritems():
                        x_vrt_ind[_ind] = _prev
                    for _ind, (_prev, _next) in _moved_vrt['y'].iteritems():
                        y_vrt_ind[_ind] = _prev
                    for _ind, (_prev, _next) in _moved_vrt['z'].iteritems():
                        z_vrt_ind[_ind] = _prev
                
                # tentative moves to its adjacent are finished
                # check the increase of modularity derived by the move
                if max_delta_mod > 0:
                    mod_val += max_delta_mod
                    mod_list.append(mod_val)
                    
                    # move the cluster to its adjacency with max increase
                    egcl_label_from_egcl[egcl_ind] = max_adj_egcl_label
                    egcl_set_from_egcl_label[egcl_label].discard(egcl_ind)
                    egcl_size_from_egcl_label[egcl_label] -= egcl_size
                    egcl_set_from_egcl_label[max_adj_egcl_label].add(egcl_ind)
                    egcl_size_from_egcl_label[max_adj_egcl_label] += egcl_size
                    
                    # assign vertex labels
                    for ind, (_prev, _next) in max_moved_vrt['x'].iteritems():
                        x_vrt_ind[ind] = _next
                    for ind, (_prev, _next) in max_moved_vrt['y'].iteritems():
                        y_vrt_ind[ind] = _next
                    for ind, (_prev, _next) in max_moved_vrt['z'].iteritems():
                        z_vrt_ind[ind] = _next
                        
                    # notify modularity that the cluster has moved
                    reflect_the_moves(max_moved_vrt, max_reflection_info)
                
            current_mod_val = mod_val
            
        return mod_val, mod_list, vrt_ind, egcl_label_from_egcl
    
    def _get_vrt_ind(self, eg_set_from_egcl, egcl_label_from_egcl, 
                     egcl_set_from_egcl_label, egcl_size_from_egcl_label,
                     adj_egcl_set_from_vrt_ind,
                     egcl_label, adj_egcl_label, vrt_ind):
        
        # detect target edges which will move
        target_egcl_set = (egcl_set_from_egcl_label[egcl_label] | 
                           egcl_set_from_egcl_label[adj_egcl_label])
        target_eg_set = set()
        for egcl_ind in target_egcl_set:
            target_eg_set |= eg_set_from_egcl[egcl_ind]
        
        adj_egcl_set_from_x_vrt_ind = adj_egcl_set_from_vrt_ind['x']
        adj_egcl_set_from_y_vrt_ind = adj_egcl_set_from_vrt_ind['y']
        adj_egcl_set_from_z_vrt_ind = adj_egcl_set_from_vrt_ind['z']
        _x_ind_list, _y_ind_list, _z_ind_list = vrt_ind
        x_moved_vrt, y_moved_vrt, z_moved_vrt = {}, {}, {}
        x_assigned_vrt, y_assigned_vrt, z_assigned_vrt = set(), set(), set()
        
        # assign vertex label
        for (_x, _y, _z) in target_eg_set:
            # for x vertices
            if _x not in x_assigned_vrt:
                x_assigned_vrt.add(_x)
                max_size = 0 
                for egcl_ind in adj_egcl_set_from_x_vrt_ind[_x]:
                    egcl_label = egcl_label_from_egcl[egcl_ind]
                    _size = egcl_size_from_egcl_label[egcl_label]
                    if _size > max_size:
                        max_size = _size
                        next_label = egcl_label
                        
                prev_label = _x_ind_list[_x]
                if prev_label != next_label:
                    _x_ind_list[_x] = next_label
                    x_moved_vrt[_x] = (prev_label, next_label)
                
            # for y vertices
            if _y not in y_assigned_vrt:
                y_assigned_vrt.add(_y)
                max_size = 0 
                for egcl_ind in adj_egcl_set_from_y_vrt_ind[_y]:
                    egcl_label = egcl_label_from_egcl[egcl_ind]
                    _size = egcl_size_from_egcl_label[egcl_label]
                    if _size > max_size:
                        max_size = _size
                        next_label = egcl_label
                        
                prev_label = _y_ind_list[_y]
                if prev_label != next_label:
                    _y_ind_list[_y] = next_label
                    y_moved_vrt[_y] = (prev_label, next_label)
                
            # for z vertices
            if _z not in z_assigned_vrt:
                z_assigned_vrt.add(_z)
                max_size = 0 
                for egcl_ind in adj_egcl_set_from_z_vrt_ind[_z]:
                    egcl_label = egcl_label_from_egcl[egcl_ind]
                    _size = egcl_size_from_egcl_label[egcl_label]
                    if _size > max_size:
                        max_size = _size
                        next_label = egcl_label
                        
                prev_label = _z_ind_list[_z]
                if prev_label != next_label:
                    _z_ind_list[_z] = next_label
                    z_moved_vrt[_z] = (prev_label, next_label)
        
        moved_vrt = {'x' : x_moved_vrt, 'y' : y_moved_vrt, 'z' : z_moved_vrt}
        return (_x_ind_list, _y_ind_list, _z_ind_list), moved_vrt

class _AbstractEdgeListTripartiteModularity():
    
    def __init__(self, edges=None):
        
        if edges is not None:
            # check data type of edges
            if isinstance(edges[0], tuple):
                self.edges = edges
            else:
                self.edges = tuple(map(tuple, edges))
                
            # the number of vertices
            xyz = zip(*edges)
            self.num_x = len(set(xyz[0]))
            self.num_y = len(set(xyz[1]))
            self.num_z = len(set(xyz[2]))
            
            # parameters which is constant about the network
            self.num_eg = len(edges)
            self.adj_eg_set_from_vrt = \
                self._make_adjacency_between_vertex_edge()
            
            # parameters which will change with the assignment of communities
            self.vrt_com = None
            self.eg_set_from_com = None
            self.eg_num_from_com = None
            self.eg_set_from_corres = None
            self.eg_num_from_corres = None
            self.mod_val_from_corres = None
            self.mod_value = None
        
    def set_edges(self, edges):
        self.__init__(edges)
    
    def name(self):
        return self.__class__.__name__
    
    def get_modularity(self, vrt_com):
        """ Calculate the value of modularity.
        vrt_com should be a list(or tuple) of indexes of 
        communities of vertices as the following. 
            vrt_com[0] = [1, 1, .... ] # labels of communities in x part
            vrt_com[1] = [1, 1, .... ] # labels of communities in y part
            vrt_com[2] = [1, 1, .... ] # labels of communities in z part
        """
        self._init_for_calc(vrt_com)
        return self._calc_modularity()
    
    def _make_adjacency_between_vertex_edge(self):
        between_vrt_eg = ({}, {}, {})
        for eg in self.edges:
            between_vrt_eg[0].setdefault(eg[0], set()).add(eg)
            between_vrt_eg[1].setdefault(eg[1], set()).add(eg)
            between_vrt_eg[2].setdefault(eg[2], set()).add(eg)
            
        return between_vrt_eg
    
    def _init_for_calc(self, vrt_com):
        """ make some dictionaries for the convenience of the calculation.
        """
        
        eg_set_from_x_com, eg_set_from_y_com, eg_set_from_z_com = {}, {}, {}
        eg_set_from_corres = {}
        
        self.vrt_com = vrt_com
        x_coms, y_coms, z_coms = self.vrt_com
        for eg in self.edges:
            x_com = x_coms[eg[0]]
            y_com = y_coms[eg[1]]
            z_com = z_coms[eg[2]]
            
            # classify the correspondences
            eg_set_from_x_com.setdefault(x_com, set()).add(eg)
            eg_set_from_y_com.setdefault(y_com, set()).add(eg)
            eg_set_from_z_com.setdefault(z_com, set()).add(eg)
            corres = (x_com, y_com, z_com)
            eg_set_from_corres.setdefault(corres, set()).add(eg)
            
        self.eg_set_from_com = \
            (eg_set_from_x_com, eg_set_from_y_com, eg_set_from_z_com)
        self.eg_set_from_corres = eg_set_from_corres
        
    def _calc_modularity(self):
        self.mod_value = sum( self._calc_partial_modularity() )
        return self.mod_value
        
    def get_num_x(self):
        return self.num_x
    
    def get_num_y(self):
        return self.num_y
    
    def get_num_z(self):
        return self.num_z
    
class NeubauerTripartiteModularityOld(_AbstractEdgeListTripartiteModularity):
    
    def _calc_partial_modularity(self):
        eg_set_from_x_com = self.eg_set_from_com[0]
        eg_set_from_y_com = self.eg_set_from_com[1]
        eg_set_from_z_com = self.eg_set_from_com[2]
        mod_val_from_corres = {}
        eg_num_from_corres = {}
        eg_num_from_com = [{}, {}, {}]
        
        for corres, eg_set in self.eg_set_from_corres.iteritems():
            x_com, y_com, z_com = corres
            eg_num_in_corres = eg_num_from_corres.setdefault(
                                                 corres, len(eg_set))
            eg_num_in_x_com = eg_num_from_com[0].setdefault(
                                                x_com, 
                                                len(eg_set_from_x_com[x_com]))
            eg_num_in_y_com = eg_num_from_com[1].setdefault(
                                                y_com, 
                                                len(eg_set_from_y_com[y_com]))
            eg_num_in_z_com = eg_num_from_com[2].setdefault(
                                                z_com, 
                                                len(eg_set_from_z_com[z_com]))
            
            e_lmn = eg_num_in_corres / self.num_eg
            a_l = eg_num_in_x_com / self.num_eg
            a_m = eg_num_in_y_com / self.num_eg
            a_n = eg_num_in_z_com / self.num_eg
            _mod = e_lmn  - a_l * a_m * a_n
            alpha = e_lmn * (1 / a_l + 1 / a_m + 1 / a_n)
            mod_val_from_corres[corres] = _mod * alpha * 1/3
            
        self.mod_val_from_corres = mod_val_from_corres
        self.eg_num_from_corres = eg_num_from_corres
        self.eg_num_from_com = eg_num_from_com
        return self.mod_val_from_corres.values()
    

    def calculate_delta_value(self, moved_vrt):
        adj_eg_set_from_x_vrt = self.adj_eg_set_from_vrt[0]
        adj_eg_set_from_y_vrt = self.adj_eg_set_from_vrt[1]
        adj_eg_set_from_z_vrt = self.adj_eg_set_from_vrt[2]
        moved_vrt_x = moved_vrt['x']
        moved_vrt_y = moved_vrt['y']
        moved_vrt_z = moved_vrt['z']
        
        # affected a_xyz, moved edges
        _delta_a_x, _delta_a_y, _delta_a_z = {}, {}, {}
        moved_eg = set()
        _iteritems = moved_vrt_x.iteritems
        for x_vrt, (x_com_prev, x_com_next) in _iteritems():
            _moved_eg = adj_eg_set_from_x_vrt[x_vrt]
            moved_eg |= _moved_eg
            
            _delta_a_x.setdefault(x_com_prev, {'add':set(), 'del':set()} )['del'].update(_moved_eg)
            _delta_a_x.setdefault(x_com_next, {'add':set(), 'del':set()} )['add'].update(_moved_eg)
            
        _iteritems = moved_vrt_y.iteritems
        for y_vrt, (y_com_prev, y_com_next) in _iteritems():
            _moved_eg = adj_eg_set_from_y_vrt[y_vrt]
            moved_eg |= _moved_eg
            
            _delta_a_y.setdefault(y_com_prev, {'add':set(), 'del':set()} )['del'].update(_moved_eg)
            _delta_a_y.setdefault(y_com_next, {'add':set(), 'del':set()} )['add'].update(_moved_eg)
            
        _iteritems = moved_vrt_z.iteritems
        for z_vrt, (z_com_prev, z_com_next) in _iteritems():
            _moved_eg = adj_eg_set_from_z_vrt[z_vrt]
            moved_eg |= _moved_eg
            
            _delta_a_z.setdefault(z_com_prev, {'add':set(), 'del':set()} )['del'].update(_moved_eg)
            _delta_a_z.setdefault(z_com_next, {'add':set(), 'del':set()} )['add'].update(_moved_eg)
        
        # all communities which are affected by the move
        eg_num_from_x_com, eg_num_from_y_com, eg_num_from_z_com = \
            self.eg_num_from_com
        new_eg_num_from_x_com, new_eg_num_from_y_com, new_eg_num_from_z_com = \
            {}, {}, {} 
        _iteritems = _delta_a_x.iteritems
        for x_com, add_del_dic in _iteritems():
            new_eg_num_from_x_com.setdefault(x_com, 
                                             (eg_num_from_x_com.get(x_com, 0) + 
                                              len(add_del_dic['add']) -
                                              len(add_del_dic['del'])) )
        
        _iteritems = _delta_a_y.iteritems
        for y_com, add_del_dic in _iteritems():
            new_eg_num_from_y_com.setdefault(y_com, 
                                             (eg_num_from_y_com.get(y_com, 0) + 
                                              len(add_del_dic['add']) -
                                              len(add_del_dic['del'])) )
        
        _iteritems = _delta_a_z.iteritems
        for z_com, add_del_dic in _iteritems():
            new_eg_num_from_z_com.setdefault(z_com, 
                                             (eg_num_from_z_com.get(z_com, 0) + 
                                              len(add_del_dic['add']) -
                                              len(add_del_dic['del'])) )
        
        # all correspondences which are affected by the move
        eg_num_from_corres = self.eg_num_from_corres
        new_eg_num_from_corres = {}
        x_vrt_com, y_vrt_com, z_vrt_com = self.vrt_com
        for _x, _y, _z in moved_eg:
            if _x in moved_vrt_x:
                prev_x_ind, next_x_ind = moved_vrt_x[_x][0], moved_vrt_x[_x][1]
            else:
                prev_x_ind, next_x_ind = x_vrt_com[_x], x_vrt_com[_x]
                
            if _y in moved_vrt_y:
                prev_y_ind, next_y_ind = moved_vrt_y[_y][0], moved_vrt_y[_y][1]
            else:
                prev_y_ind, next_y_ind = y_vrt_com[_y], y_vrt_com[_y]
                
            if _z in moved_vrt_z:
                prev_z_ind, next_z_ind = moved_vrt_z[_z][0], moved_vrt_z[_z][1]
            else:
                prev_z_ind, next_z_ind = z_vrt_com[_z], z_vrt_com[_z]
                
            prev_corres = (prev_x_ind, prev_y_ind, prev_z_ind)
            _new_eg_num = new_eg_num_from_corres.get(
                                             prev_corres, 
                                             eg_num_from_corres[prev_corres])
            new_eg_num_from_corres[prev_corres] = _new_eg_num - 1
            
            next_corres = (next_x_ind, next_y_ind, next_z_ind)
            _new_eg_num = new_eg_num_from_corres.get(
                                         next_corres, 
                                         eg_num_from_corres.get(next_corres, 0))
            new_eg_num_from_corres[next_corres] = _new_eg_num + 1
            
        # revise modularity
        mod_val_from_corres = self.mod_val_from_corres
        new_mod_val_from_corres = {}
        all_eg_num = self.num_eg
        delta = 0
        _iteritems = new_eg_num_from_corres.iteritems
        for corres, eg_num_in_corres in _iteritems():
            if eg_num_in_corres == 0:
                delta -= mod_val_from_corres.get(corres, 0)
                new_mod_val_from_corres[corres] = 0
                continue
            
            x_com, y_com, z_com = corres
            eg_num_in_x_com = new_eg_num_from_x_com.get(
                                            x_com, 
                                            eg_num_from_x_com.get(x_com, 0))
            eg_num_in_y_com = new_eg_num_from_y_com.get(
                                            y_com, 
                                            eg_num_from_y_com.get(y_com, 0))
            eg_num_in_z_com = new_eg_num_from_z_com.get(
                                            z_com, 
                                            eg_num_from_z_com.get(z_com, 0))
            
            e_lmn = eg_num_in_corres / all_eg_num
            a_l = eg_num_in_x_com / all_eg_num
            a_m = eg_num_in_y_com / all_eg_num
            a_n = eg_num_in_z_com / all_eg_num
            
            _mod = e_lmn  - a_l * a_m * a_n
            alpha = e_lmn * (1 / a_l + 1 / a_m + 1 / a_n)
            mod = _mod * alpha * 1/3
            
            delta += mod
            delta -= mod_val_from_corres.get(corres, 0)
            new_mod_val_from_corres[corres] = mod
            
        # to prevent precision loss (桁落ち回避)
        if -0.00000009 < delta < 0.00000009:
            delta = 0
        
        # may be used to reflect the moves later
        reflection_info = {'delta_mod' : delta, 
                           'new_moved_eg_info_from_com' : (_delta_a_x, 
                                                           _delta_a_y, 
                                                           _delta_a_z),
                           'moved_eg_set' : moved_eg, 
                           'new_eg_num_from_com' : (new_eg_num_from_x_com, 
                                                    new_eg_num_from_y_com, 
                                                    new_eg_num_from_z_com),
                           'new_eg_num_from_corres' : new_eg_num_from_corres,
                           'new_mod_val_from_corres' : new_mod_val_from_corres}
        
        return delta, reflection_info
        
    def reflect_the_moves(self, moved_vrt, reflection_info):
        moved_vrt_x = moved_vrt['x']
        moved_vrt_y = moved_vrt['y']
        moved_vrt_z = moved_vrt['z']
        delta_mod_val = reflection_info['delta_mod']
        new_moved_eg_info_from_com = reflection_info['new_moved_eg_info_from_com']
        moved_eg_set = reflection_info['moved_eg_set']
        new_eg_num_from_com = reflection_info['new_eg_num_from_com']
        new_eg_num_from_corres = reflection_info['new_eg_num_from_corres']
        new_mod_val_from_corres = reflection_info['new_mod_val_from_corres']
        
        # reflect community information
        new_moved_eg_info_from_x_com = new_moved_eg_info_from_com[0]
        new_moved_eg_info_from_y_com = new_moved_eg_info_from_com[1]
        new_moved_eg_info_from_z_com = new_moved_eg_info_from_com[2]
        new_eg_num_from_x_com = new_eg_num_from_com[0]
        new_eg_num_from_y_com = new_eg_num_from_com[1]
        new_eg_num_from_z_com = new_eg_num_from_com[2]
        eg_set_from_x_com = self.eg_set_from_com[0]
        eg_set_from_y_com = self.eg_set_from_com[1]
        eg_set_from_z_com = self.eg_set_from_com[2]
        eg_num_from_x_com = self.eg_num_from_com[0]
        eg_num_from_y_com = self.eg_num_from_com[1]
        eg_num_from_z_com = self.eg_num_from_com[2]
        _iteritems = new_eg_num_from_x_com.iteritems
        for x_com, eg_num in _iteritems():
            if eg_num == 0:
                del eg_set_from_x_com[x_com]
                del eg_num_from_x_com[x_com]
            else:
                _eg_set = (eg_set_from_x_com.get(x_com, set()) |
                           new_moved_eg_info_from_x_com[x_com]['add'] -
                           new_moved_eg_info_from_x_com[x_com]['del'])
                eg_set_from_x_com[x_com] = _eg_set
                eg_num_from_x_com[x_com] = eg_num
        
        _iteritems = new_eg_num_from_y_com.iteritems
        for y_com, eg_num in _iteritems():
            if eg_num == 0:
                del eg_set_from_y_com[y_com]
                del eg_num_from_y_com[y_com]
            else:
                _eg_set = (eg_set_from_y_com.get(y_com, set()) |
                           new_moved_eg_info_from_y_com[y_com]['add'] -
                           new_moved_eg_info_from_y_com[y_com]['del'])
                eg_set_from_y_com[y_com] = _eg_set
                eg_num_from_y_com[y_com] = eg_num
        
        _iteritems = new_eg_num_from_z_com.iteritems
        for z_com, eg_num in _iteritems():
            if eg_num == 0:
                del eg_set_from_z_com[z_com]
                del eg_num_from_z_com[z_com]
            else:
                _eg_set = (eg_set_from_z_com.get(z_com, set()) |
                           new_moved_eg_info_from_z_com[z_com]['add'] -
                           new_moved_eg_info_from_z_com[z_com]['del'])
                eg_set_from_z_com[z_com] = _eg_set
                eg_num_from_z_com[z_com] = eg_num
                
        # reflect correspondence information
        x_vrt_com, y_vrt_com, z_vrt_com = self.vrt_com
        eg_set_from_corres = self.eg_set_from_corres
        eg_num_from_corres = self.eg_num_from_corres
        mod_val_from_corres = self.mod_val_from_corres
        for eg in moved_eg_set:
            _x, _y, _z = eg
            if _x in moved_vrt_x:
                prev_x_ind, next_x_ind = moved_vrt_x[_x][0], moved_vrt_x[_x][1]
            else:
                prev_x_ind, next_x_ind = x_vrt_com[_x], x_vrt_com[_x]
                
            if _y in moved_vrt_y:
                prev_y_ind, next_y_ind = moved_vrt_y[_y][0], moved_vrt_y[_y][1]
            else:
                prev_y_ind, next_y_ind = y_vrt_com[_y], y_vrt_com[_y]
                
            if _z in moved_vrt_z:
                prev_z_ind, next_z_ind = moved_vrt_z[_z][0], moved_vrt_z[_z][1]
            else:
                prev_z_ind, next_z_ind = z_vrt_com[_z], z_vrt_com[_z]
                
            prev_corres = (prev_x_ind, prev_y_ind, prev_z_ind)
            eg_set_from_corres[prev_corres].discard(eg)
            eg_num_from_corres[prev_corres] = new_eg_num_from_corres[prev_corres]
            mod_val_from_corres[prev_corres] = new_mod_val_from_corres[prev_corres]
            
            next_corres = (next_x_ind, next_y_ind, next_z_ind)
            eg_set_from_corres.setdefault(next_corres, set()).add(eg)
            eg_num_from_corres[next_corres] = new_eg_num_from_corres[next_corres]
            mod_val_from_corres[next_corres] = new_mod_val_from_corres[next_corres]
        
        # reflect vertex indexes
        _iteritems = moved_vrt_x.iteritems
        for x_vrt, (_, next_x_com) in _iteritems():
            x_vrt_com[x_vrt] = next_x_com
        _iteritems = moved_vrt_y.iteritems
        for y_vrt, (_, next_y_com) in _iteritems():
            y_vrt_com[y_vrt] = next_y_com
        _iteritems = moved_vrt_z.iteritems
        for z_vrt, (_, next_z_com) in _iteritems():
            z_vrt_com[z_vrt] = next_z_com
            
        # reflect modularity value
        self.mod_value += delta_mod_val
        
        

if __name__ == '__main__':
    from modularity import NeubauerModularity
    modularity = NeubauerModularity()

    import time
    from npartite.synthetic.tripartite import SimpleCaseMaker
    corres_egnum, noise, is_scale_free = 33, 0.0, False
    
    syn_maker = SimpleCaseMaker(corres_egnum, noise, is_scale_free)
    syn_network = syn_maker.make()
    
    syn_name = syn_network.name_with_config()
    eglist = syn_network.edge_list()
    modularity = NeubauerModularity()
    
    optimization = TriFastUnfoldingForEdges()
    print 'Old implementation, fue', syn_name, optimization.name()
    s = time.time()
    ans_modval, ans_ind, modval_list = optimization.start(modularity, eglist)
    e = time.time()
    print '%f [s]' % (e-s)
    print ans_modval
    print ans_ind
    print modval_list
    print ''


