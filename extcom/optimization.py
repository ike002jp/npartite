#!/usr/bin/env python
# encoding: utf-8


#--------------------------------------------------
# modules
#--------------------------------------------------
from __future__ import division
import copy

from _status import NetworkStatus
from _status import ComEgclSynchronalManeger

#--------------------------------------------------
# static variables
#--------------------------------------------------


#--------------------------------------------------
# private classes
#--------------------------------------------------


#--------------------------------------------------
# public classes
#--------------------------------------------------
class _AbstractOptimization(object):

    def name(self):
        return self.__class__.__name__

class GreedyVertexBottomUp(_AbstractOptimization):

    def __init__(self):
        pass

    def start(self, modularity, edge_list):
        status = NetworkStatus(edge_list)
        status.add_com()
        
        # get 'n' of 'n-partite network'
        partnum = status.basic.partnum()

        # initialize
        status.com.assign_unique_com_labels()
        modval = modularity.calculate(status)
        modvals = [modval]
        ans_modval = modval
        ans_labels = copy.deepcopy( status.com.com_labels() )

        # repeat merging processes
        while( status.com.count_all_comnums() > partnum ):

            # all combinations of communities
            max_merged_modval = -1
            max_merged_coms = None
            for part in xrange(partnum):
                for com1, com2 in status.com.iter_com_combinations(
                                    part, combi_num=2):
                    status.com.merge_coms_tentatively(part, com1, com2)
                    merged_modval = modularity.calculate(status)
                    if merged_modval > max_merged_modval:
                        max_merged_modval = merged_modval
                        max_merged_coms = (part, com1, com2)
                    status.com.rollback_merging_coms()

            # maximum increase of modularity
            modval = max_merged_modval
            modvals.append(modval)

            # update the network with the maximum increase of merging
            status.com.merge_coms(*max_merged_coms)

            # maximum modularity during the algorithm
            if modval > ans_modval:
                ans_modval = modval
                ans_labels = copy.deepcopy( status.com.com_labels() )

        return ans_modval, ans_labels, modvals

class GreedyVertexBottomUpSpeedy(_AbstractOptimization):

    def __init__(self):
        pass

    def start(self, modularity, edge_list):
        status = NetworkStatus(edge_list)
        status.add_com()
        
        # get 'n' of 'n-partite network'
        partnum = status.basic.partnum()

        # initialize
        status.com.assign_unique_com_labels()
        modval = modularity.calculate(status)
        modvals = [modval]
        ans_modval = modval
        ans_labels = copy.deepcopy( status.com.com_labels() )

        # repeat merging processes
        while( status.com.count_all_comnums() > partnum ):

            # all combinations of communities
            max_merged_diff = -1
            max_merged_coms = None
            max_moving_diff_info = None
            max_modval_diff_info = None
            for part in xrange(partnum):
                for com1, com2 in status.com.iter_com_combinations(
                                    part, combi_num=2):
                    moving_diff_info = status.com.diff_of_merging_coms(part, com1, com2)
                    diff_modval, modval_diff_info = modularity.calculate_diff(status, moving_diff_info)
                    if diff_modval > max_merged_diff:
                        max_merged_diff = diff_modval
                        max_merged_coms = (part, com1, com2)
                        max_moving_diff_info = moving_diff_info
                        max_modval_diff_info = modval_diff_info

            # maximum increase of modularity
            modval += max_merged_diff
            modvals.append(modval)

            # update the network with the maximum increase of merging
            status.com.update_com_with_diff_info(max_moving_diff_info)
            modularity.update_modval_with_diff_info(max_modval_diff_info)

            # maximum modularity during the algorithm
            if modval > ans_modval:
                ans_modval = modval
                ans_labels = copy.deepcopy( status.com.com_labels() )

        return ans_modval, ans_labels, modvals

class GreedyEdgeBottomUp(_AbstractOptimization):

    def __init__(self):
        pass

    def start(self, modularity, edge_list):
        return self._start(modularity, edge_list)

    def _start(self, modularity, edge_list):
        # initialize
        status = NetworkStatus(edge_list)
        status.add_com()
        status.add_egcl()
        status_maneger = ComEgclSynchronalManeger(status)

        status_maneger.assign_unique_egcl_labels()
        modval = modularity.calculate(status)
        modvals = [modval]
        ans_modval = modval
        ans_labels = copy.deepcopy( status.com.com_labels() )

        # repeat merging processes
        while( status.egcl.calculate_egclnum() > 1 ):

            # all combinations of community 
            max_merged_modval = -1
            max_merged_egcls = None
            for egcl1, egcl2 in status.egcl.iter_egcl_combination():
                status_maneger.merge_egcls_tentatively(egcl1, egcl2)
                #status.egcl.merge_egcls_tentatively_naively(egcl1, egcl2)
                merged_modval = modularity.calculate(status)
                if merged_modval > max_merged_modval:
                    max_merged_modval = merged_modval
                    max_merged_egcls = (egcl1, egcl2)
                status_maneger.rollback_merging_egcls()
                #status.egcl.rollback_tentative_egcls_merge_naively()

            # maximum increase of modularity
            modval = max_merged_modval
            modvals.append(modval)

            # update the network with the maximum increase of merging
            status_maneger.merge_egcls(*max_merged_egcls)

            # maximum modularity during the algorithm
            if modval > ans_modval:
                ans_modval = modval
                ans_labels = copy.deepcopy( status.com.com_labels() )
        return ans_modval, ans_labels, modvals

class FastUnfoldingForEdgesNaively(_AbstractOptimization):
    def __init__(self):
        pass

    def start(self, modularity, edge_list):
        return self._start(modularity, edge_list)

    def _start(self, modularity, edge_list):
        # initialize
        status = NetworkStatus(edge_list)
        status.add_com()
        status.add_egcl(is_hierarchical=True)

        updater = ComEgclSynchronalManeger(status)

        updater.assign_unique_egcl_labels()
        ans_modval = modularity.calculate(status)
        modvals = []

        while True:
            _modval, _modvals = \
                self._optimize_modularity(modularity, status, updater)

            if _modval <= ans_modval:
                break

            status.egcl.merge_egcls_hierarchically()
            
            ans_modval = _modval
            modvals.extend(_modvals)

        return ans_modval, status.com.com_labels(), modvals
        

    def _optimize_modularity(self, modularity, status, updater):
        prev_modval = -1
        modval = modularity.calculate(status)
        modvals = [modval]
        
        while modval > prev_modval:
            prev_modval = modval

            for egcl in status.egcl.egcls_randomly():
                egcl_label = status.egcl.egcl_label(egcl)

                max_moved_modval = -1
                max_moved_egcl_label = None
                for adj_egcl in status.egcl.adj_egclset(egcl):
                    adj_egcl_label = status.egcl.egcl_label(adj_egcl)
                    if egcl_label == adj_egcl_label:
                        continue

                    updater.move_egcl_tentatively(egcl, adj_egcl_label)
                    moved_modval = modularity.calculate(status)
                    if moved_modval > max_moved_modval:
                        max_moved_modval = moved_modval
                        max_moved_egcl_label = adj_egcl_label
                    updater.rollback_moving_egcl()

                if max_moved_modval > modval:
                    modval = max_moved_modval
                    modvals.append(modval)
                    updater.move_egcl(egcl, max_moved_egcl_label)

        return modval, modvals

class FastUnfoldingForEdges(FastUnfoldingForEdgesNaively):

    def _optimize_modularity(self, modularity, status, updater):
        prev_modval = -1
        modval = modularity.calculate(status)
        modvals = [modval]
        
        while modval > prev_modval:
            prev_modval = modval

            for egcl in status.egcl.egcls_randomly():
                egcl_label = status.egcl.egcl_label(egcl)

                max_diff_modval = -1
                max_moved_egcl_label = None
                max_moving_diff_info = None
                max_modval_diff_info = None
                for adj_egcl in status.egcl.adj_egclset(egcl):
                    adj_egcl_label = status.egcl.egcl_label(adj_egcl)
                    if egcl_label == adj_egcl_label:
                        continue

                    #moving_diff_info = updater.diff_of_moving_egcl(
                    #    egcl, adj_egcl_label)
                    moving_diff_info = updater.minimal_diff_of_moving_egcl(
                        egcl, adj_egcl_label)
                    diff_modval, modval_diff_info = modularity.calculate_diff(
                        status, moving_diff_info)
                    if diff_modval > max_diff_modval:
                        max_diff_modval = diff_modval
                        max_moved_egcl_label = adj_egcl_label
                        max_moving_diff_info = moving_diff_info
                        max_modval_diff_info = modval_diff_info

                if max_diff_modval > 0:
                    modval += max_diff_modval
                    modvals.append(modval)

                    max_moving_diff_info = updater.diff_of_moving_egcl(
                        egcl, max_moved_egcl_label)
                     
                    status.egcl.move_egcl_tentatively(egcl, max_moved_egcl_label)
                    status.com.update_com_with_diff_info(max_moving_diff_info)
                    modularity.update_modval_with_diff_info(max_modval_diff_info)

        return modval, modvals


def _test_optimizing(eglist, opt, mod, answer, in_detail=True):
    import time

    print '  %-50s :' % ( '%s, %s' % (opt.name(), mod.name()), ), 
    ans_modval = _adjust_value( answer['modval'] )
    ans_labels = answer['labels']
    ans_modvals = _adjust_values( answer['modvals'] )
    
    s = time.time()
    _modval, labels, _modvals = opt.start(mod, eglist)
    e = time.time()
    modval = _adjust_value(_modval)
    modvals = _adjust_values(_modvals)

    same_mod = (modval == ans_modval)
    same_labels = (labels == ans_labels)
    same_modvals = (modvals == ans_modvals)
    same = same_mod and same_labels and same_modvals
    print '%s (modval:%s, labels:%s, modvals:%s), %f [s]' % \
              (same, same_mod, same_labels, same_modvals, e-s)
    if in_detail:
        print '    ans modval  :', answer['modval']
        print '        modval  :', _modval
        print '    ans labels  :', answer['labels']
        print '        labels  :', labels
        print '    ans modvals :', answer['modvals']
        print '        modvals :', _modvals
        print ''

def _adjust_values(values):
    return map( _adjust_value, values)

def _adjust_value(value):
    return round(value, 6)

if __name__ == '__main__':
    from modularity import MurataModularity
    from modularity import NeubauerModularity
    in_detail = False
    #in_detail = True
    
    murata = MurataModularity()
    neubauer = NeubauerModularity()
    
    vertex_bottom_up = GreedyVertexBottomUp()
    vertex_bottom_up_speedy = GreedyVertexBottomUpSpeedy()
    edge_bottom_up = GreedyEdgeBottomUp()
    fue_naively = FastUnfoldingForEdgesNaively()
    fue = FastUnfoldingForEdges()

    print 'Small network test'
    eglist = [[0, 0, 0], [0, 0, 1], [1, 1, 0], [1, 1, 1],
              [2, 2, 2], [2, 2, 3], [3, 3, 2], [3, 3, 3]]
    print len(eglist), eglist

    # greedy vertex, Murata
    answer = {'modval' : 0.729166666667,
              'labels' : [[0, 1, 2, 3], [0, 1, 2, 3], [1, 1, 3, 3]],
              'modvals': [0.4375, 0.5833333333333334, 0.7291666666666666, 0.625, 0.5104166666666666, 0.3958333333333333, 0.375, 0.20833333333333334, 0.0, 0.0],
              }
    _test_optimizing(eglist, vertex_bottom_up, murata, answer, in_detail)

    # greedy vertex, Neubauer
    answer = {'modval' : 0.75,
              'labels' : [[1, 1, 3, 3], [1, 1, 3, 3], [1, 1, 3, 3]],
              'modvals': [0.4375, 0.5833333333333333, 0.7291666666666666, 0.6145833333333333, 0.7395833333333333, 0.625, 0.75, 0.4166666666666667, 0.0, 0.0],
              }
    _test_optimizing(eglist, vertex_bottom_up, neubauer, answer, in_detail)

    # greedy vertex speedy, Neubauer
    answer = {'modval' : 0.75,
              'labels' : [[1, 1, 3, 3], [1, 1, 3, 3], [1, 1, 3, 3]],
              'modvals': [0.4375, 0.5833333333333333, 0.7291666666666666, 0.6145833333333333, 0.7395833333333333, 0.625, 0.75, 0.4166666666666667, 0.0, 0.0],
              }
    _test_optimizing(eglist, vertex_bottom_up, neubauer, answer, in_detail)

    # greedy edge, Murata
    answer = {'modval' : 0.75,
              'labels' : [[3, 3, 7, 7], [3, 3, 7, 7], [3, 3, 7, 7]],
              'modvals': [0.4375, 0.59375, 0.75, 0.75, 0.75, 0.75, 0.75, 0.0],
              }
    _test_optimizing(eglist, edge_bottom_up, murata, answer, in_detail)

    # greedy edge, Neubauer
    answer = {'modval' : 0.75,
              'labels' : [[3, 3, 7, 7], [3, 3, 7, 7], [3, 3, 7, 7]],
              'modvals': [0.4375, 0.59375, 0.75, 0.75, 0.75, 0.75, 0.75, 0.0],
              }
    _test_optimizing(eglist, edge_bottom_up, neubauer, answer, in_detail)

    # FUE, Neubauer
    answer = {'modval' : 0.75,
              'labels' : [[3, 3, 7, 7], [3, 3, 7, 7], [3, 3, 7, 7]],
              'modvals': [0.4375, 0.59375, 0.75, 0.75, 0.75, 0.75, 0.75, 0.0],
              }
    _test_optimizing(eglist, fue_naively, neubauer, answer, in_detail)
    #_test_optimizing(eglist, fue_naively, neubauer, answer, True)

    # FUE, Neubauer
    answer = {'modval' : 0.75,
              'labels' : [[3, 3, 7, 7], [3, 3, 7, 7], [3, 3, 7, 7]],
              'modvals': [0.4375, 0.59375, 0.75, 0.75, 0.75, 0.75, 0.75, 0.0],
              }
    _test_optimizing(eglist, fue, neubauer, answer, in_detail)
    #_test_optimizing(eglist, fue, neubauer, answer, True)

    print ''


    print 'Large network test (Simple Case, vrt-60 eg-66 noise-0.000 scale-False)'
    eglist = [(15, 11, 13), (1, 8, 4), (5, 7, 7), (18, 10, 15), (3, 2, 1), (17, 18, 11), (8, 1, 5), (13, 13, 14), (0, 9, 0), (14, 17, 16), (6, 6, 2), (11, 12, 19), (12, 19, 10), (4, 0, 8), (19, 15, 18), (2, 4, 6), (10, 14, 17), (7, 5, 9), (9, 3, 3), (16, 16, 12), (14, 10, 14), (8, 5, 3), (16, 19, 18), (8, 2, 6), (14, 18, 14), (11, 14, 17), (8, 1, 3), (17, 13, 14), (13, 12, 17), (1, 5, 0), (12, 18, 12), (3, 5, 7), (14, 13, 11), (15, 19, 12), (13, 15, 19), (18, 11, 14), (12, 18, 13), (10, 11, 12), (17, 13, 18), (7, 1, 9), (15, 16, 16), (5, 5, 2), (10, 16, 16), (14, 15, 11), (8, 9, 3), (5, 9, 5), (19, 19, 17), (9, 3, 4), (17, 10, 18), (10, 14, 16), (10, 19, 13), (7, 2, 2), (17, 16, 14), (10, 16, 11), (10, 18, 10), (7, 7, 7), (9, 6, 4), (4, 2, 8), (17, 17, 12), (13, 11, 16), (17, 11, 12), (18, 13, 17), (10, 17, 16), (12, 12, 18), (18, 14, 17), (18, 14, 12)]
    print len(eglist), eglist

    # greedy vertex, Murata
    #answer = {'modval' : 0.547520661157,
    #          'labels' : [[19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 18, 18, 19, 18, 19, 18, 18, 18, 18, 19], [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], [15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 19, 19, 19, 19, 15, 19, 19, 19, 19]],
    #          'modvals': [0.2988841583882907, 0.3190861785903109, 0.3292242906567976, 0.3390760219272617, 0.3489428258712006, 0.3587667306675573, 0.3684132416914787, 0.3783356985836326, 0.3879300349686025, 0.40584564654812594, 0.42492649173089964, 0.4342298095741623, 0.4431180491786555, 0.4499691589911978, 0.45846319021250176, 0.4633409392362562, 0.46771665228965525, 0.4707706078229497, 0.470207121722273, 0.4868960495682257, 0.49501674226192116, 0.5004556585135097, 0.5059826919331051, 0.5105323204496757, 0.5125682908052054, 0.5113961005834283, 0.5163909526857184, 0.5143712144400848, 0.533879891662261, 0.5354288987209098, 0.5445258368812088, 0.5444330819675173, 0.5408933225737634, 0.5404318668781479, 0.5354775950505978, 0.5304734674569385, 0.5254415133891718, 0.5252931055272653, 0.5151920954262553, 0.5244212093385646, 0.5143758521857695, 0.504015128326423, 0.502456845776405, 0.4913911845730026, 0.4901297177467976, 0.478266364285648, 0.4791579708935081, 0.47395442023541196, 0.473954420235412, 0.4455053287697915, 0.40729957982024095, 0.45065322647967276, 0.510034922225005, 0.4786779642151542, 0.5475206611570248, 0.3581267217630854, 0.0, 0.0],
    #          }
    answer = {'modval' : 0.683195592287,
              'labels' : [[9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19], [15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 19, 19, 19, 19, 19, 15, 19, 19, 19, 19]],
              'modvals': [0.2988841583882907, 0.3190861785903109, 0.32903414308373, 0.3389623971579894, 0.34864601014738744, 0.3664387214662697, 0.385170576286279, 0.39523796273107603, 0.4051766517331259, 0.4149275120349502, 0.4245485154576065, 0.4339121239947688, 0.44313660016139383, 0.44984857760339875, 0.45843536373839444, 0.4633873166931021, 0.46771085510754945, 0.4724483123243453, 0.4715509085343797, 0.4703787183126027, 0.4753735704148927, 0.47335383216925914, 0.4928625093914352, 0.49441151645008413, 0.5035084546103831, 0.5034156996966915, 0.5010666815074529, 0.5047919507285898, 0.5027536615002178, 0.5057519640852973, 0.5011721902217771, 0.4986817207891588, 0.5060232722078453, 0.5074006826761647, 0.5092465054586266, 0.5042110730815964, 0.4991049150828764, 0.49903071115192327, 0.48885549711996007, 0.4983907022474517, 0.4883360696032873, 0.4781515800799547, 0.47282744803405974, 0.5045009321868826, 0.5344004786153547, 0.5240026527905315, 0.5329163999962898, 0.5228153898952796, 0.5106366697275787, 0.47476834460305545, 0.5294473662242257, 0.5982344102178813, 0.5587486434593872, 0.6138937585218577, 0.6831955922865013, 0.38567493112947665, 0.0, 0.0],
              }
    _test_optimizing(eglist, vertex_bottom_up, murata, answer, in_detail)
    #_test_optimizing(eglist, vertex_bottom_up, murata, answer, True)

    # greedy vertex, Neubauer
    answer = {'modval' : 0.475523306395,
              'labels' : [[1, 1, 2, 3, 4, 8, 9, 8, 8, 9, 19, 18, 18, 19, 19, 19, 18, 19, 18, 19], [9, 9, 9, 6, 4, 9, 6, 7, 8, 9, 18, 19, 19, 18, 18, 18, 18, 18, 18, 19], [0, 7, 9, 9, 9, 9, 9, 7, 8, 9, 19, 19, 17, 17, 19, 17, 19, 17, 19, 19]],
              'modvals': [0.2987473448905957, 0.30811141720232627, 0.315139225125451, 0.3199889157878138, 0.3260170576286279, 0.3314091326488021, 0.33669964103848404, 0.3409350622849246, 0.345064974817041, 0.34609107604975375, 0.356080780254334, 0.356519510996095, 0.3565410765135285, 0.3633233158026548, 0.367177514353823, 0.3685799686488393, 0.37088831737592076, 0.37141553095409846, 0.37408373052312455, 0.387811408059342, 0.38896971817082027, 0.3882046557666394, 0.38923026009802875, 0.38802440799685983, 0.3904468911355964, 0.39128680726476867, 0.3950203581746281, 0.39317045570489095, 0.39280733805527196, 0.3951514635277397, 0.3934078584439142, 0.38943530819037214, 0.39237146111185695, 0.41232907904505867, 0.4516643077770766, 0.44673554355079176, 0.4635261705342562, 0.46437822379677784, 0.47332443522232814, 0.4755233063950285, 0.4704310616333623, 0.47262065731469344, 0.4623526883690386, 0.45207544393201454, 0.4609636835365076, 0.4486226422698465, 0.45790740913037087, 0.4419767527038467, 0.44019585836096886, 0.40507828743072577, 0.36515889145072367, 0.418465471616808, 0.46451072978327834, 0.3993764175914837, 0.4533854451209824, 0.3856749311294766, 0.0, 0.0],
              }
    _test_optimizing(eglist, vertex_bottom_up, neubauer, answer, in_detail)

    # greedy vertex speedy, Neubauer
    answer = {'modval' : 0.475523306395,
              'labels' : [[1, 1, 2, 3, 4, 8, 9, 8, 8, 9, 19, 18, 18, 19, 19, 19, 18, 19, 18, 19], [9, 9, 9, 6, 4, 9, 6, 7, 8, 9, 18, 19, 19, 18, 18, 18, 18, 18, 18, 19], [0, 7, 9, 9, 9, 9, 9, 7, 8, 9, 19, 19, 17, 17, 19, 17, 19, 17, 19, 19]],
              'modvals': [0.2987473448905957, 0.30811141720232627, 0.315139225125451, 0.3199889157878138, 0.3260170576286279, 0.3314091326488021, 0.33669964103848404, 0.3409350622849246, 0.345064974817041, 0.34609107604975375, 0.356080780254334, 0.356519510996095, 0.3565410765135285, 0.3633233158026548, 0.367177514353823, 0.3685799686488393, 0.37088831737592076, 0.37141553095409846, 0.37408373052312455, 0.387811408059342, 0.38896971817082027, 0.3882046557666394, 0.38923026009802875, 0.38802440799685983, 0.3904468911355964, 0.39128680726476867, 0.3950203581746281, 0.39317045570489095, 0.39280733805527196, 0.3951514635277397, 0.3934078584439142, 0.38943530819037214, 0.39237146111185695, 0.41232907904505867, 0.4516643077770766, 0.44673554355079176, 0.4635261705342562, 0.46437822379677784, 0.47332443522232814, 0.4755233063950285, 0.4704310616333623, 0.47262065731469344, 0.4623526883690386, 0.45207544393201454, 0.4609636835365076, 0.4486226422698465, 0.45790740913037087, 0.4419767527038467, 0.44019585836096886, 0.40507828743072577, 0.36515889145072367, 0.418465471616808, 0.46451072978327834, 0.3993764175914837, 0.4533854451209824, 0.3856749311294766, 0.0, 0.0],
              }
    _test_optimizing(eglist, vertex_bottom_up_speedy, neubauer, answer, in_detail)

    # FUE naively, Neubauer
    answer = {'modval' : 0.75,
              'labels' : [[3, 3, 7, 7], [3, 3, 7, 7], [3, 3, 7, 7]],
              'modvals': [0.4375, 0.59375, 0.75, 0.75, 0.75, 0.75, 0.75, 0.0],
              }
    _test_optimizing(eglist, fue_naively, neubauer, answer, in_detail)
    #_test_optimizing(eglist, fue_naively, neubauer, answer, True)

    # FUE, Neubauer
    answer = {'modval' : 0.75,
              'labels' : [[3, 3, 7, 7], [3, 3, 7, 7], [3, 3, 7, 7]],
              'modvals': [0.4375, 0.59375, 0.75, 0.75, 0.75, 0.75, 0.75, 0.0],
              }
    _test_optimizing(eglist, fue, neubauer, answer, in_detail)
    #_test_optimizing(eglist, fue, neubauer, answer, True)

    # greedy edge, Murata
    modularity = MurataModularity()
    answer = {'modval' : 0.6480924952,
              'labels' : [[56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 65, 65, 63, 48, 48, 65, 65, 48, 65, 65], [56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 48, 65, 63, 48, 65, 48, 65, 65, 48, 65], [56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 16, 48, 65, 65, 48, 63, 65, 65, 48, 48]],
              'modvals': [0.29844473198467686, 0.3331153592861583, 0.36075864243908323, 0.39004484700076997, 0.4176939273358007, 0.4426786227750416, 0.46443892552707994, 0.48401253118883997, 0.5013704538497926, 0.5171503835415683, 0.5342636651176597, 0.5489073471167135, 0.5604391017614156, 0.5716462142081977, 0.5814596840767638, 0.5883884761295229, 0.5961056849486602, 0.6012767713869643, 0.6182996169222064, 0.6242927437831018, 0.6320111120386602, 0.6380146738273461, 0.6440947584198274, 0.6462652234002096, 0.6468959568133122, 0.6474895882609382, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.6480924951999333, 0.623980855385814, 0.6238660711801207, 0.590885902180668, 0.5889194980104071, 0.4639519622301991, 0.44164788379664416, 0.11486304736993445, 0.0],
              }
    _test_optimizing(eglist, edge_bottom_up, murata, answer, in_detail)
    #_test_optimizing(eglist, edge_bottom_up, murata, answer, True)

    # greedy edge, Neubauer
    answer = {'modval' : 0.694214876033,
              'labels' : [[56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64], [56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64], [56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64]],
              'modvals': [0.2987473448905955, 0.31282972715479607, 0.3248787229503484, 0.3439189878583818, 0.35955042050083386, 0.3716607568535945, 0.3821757520104629, 0.4078052771687688, 0.41460908197532836, 0.4205419181423243, 0.42462545321759493, 0.42852881187301956, 0.44742426317205203, 0.4573292125742466, 0.474449884297266, 0.4876696557246071, 0.5163991705955955, 0.5526099846900849, 0.5699659936828271, 0.5943459823145933, 0.6367729348554345, 0.6485077565076116, 0.6523287097449777, 0.6674840693435736, 0.6774146422906755, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.6942148760330579, 0.3603706022576546, 0.33128041023357824, 0.055348478581554554, 0.0],
              }
    _test_optimizing(eglist, edge_bottom_up, neubauer, answer, in_detail)
    #_test_optimizing(eglist, edge_bottom_up, neubauer, answer, True)

    print ''



