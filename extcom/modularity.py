#!/usr/bin/env python
# encoding: utf-8


#--------------------------------------------------
# modules
#--------------------------------------------------
from __future__ import division

#--------------------------------------------------
# static variables
#--------------------------------------------------


#--------------------------------------------------
# private classes
#--------------------------------------------------
class _AbstractModularity(object):
    
    def __init__(self):
        self._modval = None

    def name(self):
        return self.__class__.__name__

    def value(self):
        return self._modval
    
    def calculate(self, status):
        """ Calculate the value of modularity.
        """
        self._modval = self._calculate(status)
        return self._modval

    def _calculate(self, status):
        msg = 'modularity class should have _calculate function.'
        raise NotImplementedError(msg)

#--------------------------------------------------
# public classes
#--------------------------------------------------
class MurataModularity(_AbstractModularity):

    def _calculate(self, status):
        partnum = status.basic.get_partnum()
        total_egnum = status.basic.get_edge_num()

        modval_from_com = [{} for _ in range(partnum)]
        max_e_from_com = [{} for _ in range(partnum)]
        modval = 0
        for corres, egnum_in_corres in status.com.iter_corres_egnum():
            e = egnum_in_corres / total_egnum

            a = [None for _ in range(partnum)]
            for part, com in enumerate(corres):
                a[part] = status.com.egnum_from_com(part, com) / total_egnum

            a_multiplied = reduce(lambda x, y: x*y, a)
            partial_mod = e  - a_multiplied
            
            # whether max elmn or not for each community
            for part, com in enumerate(corres):
                max_e = max_e_from_com[part].get(com, 0)
                if e > max_e:
                    max_e_from_com[part][com] = e
                    modval -= modval_from_com[part].get(com, 0)
                    modval_from_com[part][com] = partial_mod
                    modval += partial_mod
        
        return modval / partnum
    '''
    def _calculate(self, status):
        total_eg_num = status.basic.get_edge_num()
        max_elmn_from_com = ({}, {}, {})
        modval = 0
        checked_coms = (set(), set(), set())
        for corres, eg_num_in_corres in status.com.iter_corres_egnum():
            e_lmn = eg_num_in_corres / total_eg_num
            
            # whether max elmn or not for each community
            for part, com in enumerate(corres):
                if com in checked_coms[part]:
                    continue

                max_elmn = max_elmn_from_com[part].get(com, 0)
                if e_lmn < max_elmn:
                    continue

                checked_coms[part].add(com)
                max_elmn_from_com[part][com] = e_lmn

                corres_set = status.com.corresset_from_com(part, com)
                coms_in_parts = map(set, zip(*corres_set))
                for xcom, ycom, zcom in itertools.product(*coms_in_parts):
                    _corres = (xcom, ycom, zcom)
                    corres_egnum = status.com.egnum_from_corres(_corres)
                    egnum_xcom = status.com.egnum_from_com(0, xcom)
                    egnum_ycom = status.com.egnum_from_com(1, ycom)
                    egnum_zcom = status.com.egnum_from_com(2, zcom)

                    _elmn = corres_egnum / total_eg_num
                    _al = egnum_xcom / total_eg_num
                    _am = egnum_ycom / total_eg_num
                    _an = egnum_zcom / total_eg_num
                    part_mod = _elmn - _al * _am * _an
                    modval += part_mod

        return modval / 3
    '''

class NeubauerModularity(_AbstractModularity):

    def __init__(self):
        _AbstractModularity.__init__(self)
        self._modval_from_corres = None

    def _calculate(self, status):
        partnum = status.basic.get_partnum()
        total_egnum = status.basic.get_edge_num()

        modval_from_corres = {}
        mod = 0
        for corres, egnum_in_corres in status.com.iter_corres_egnum():
            e = egnum_in_corres / total_egnum

            a = [None for _ in range(partnum)]
            a_inv = [None for _ in range(partnum)]
            for part, com in enumerate(corres):
                a[part] = status.com.egnum_from_com(part, com) / total_egnum
                a_inv[part] = 1 / a[part]

            a_multiplied = reduce(lambda x, y: x*y, a)
            alpha = e * sum(a_inv) * (1 / partnum)
            _mod = e  - a_multiplied
            partial_mod = _mod * alpha

            mod += partial_mod
            modval_from_corres[corres] = partial_mod
            
        self._modval_from_corres = modval_from_corres
        return mod

    def calculate_diff(self, status, moving_diff_info):
        #new_egset_from_com = moving_diff_info['new_egset_from_com']
        new_egnum_from_com = moving_diff_info['new_egnum_from_com']
        #new_egset_from_corres = moving_diff_info['new_egset_from_corres']
        new_egnum_from_corres = moving_diff_info['new_egnum_from_corres']

        # revise modularity
        partnum = status.basic.get_partnum()
        modval_from_corres = self._modval_from_corres
        new_modval_from_corres = {}
        all_egnum = status.basic.get_edge_num()
        delta = 0
        for corres, egnum in new_egnum_from_corres.iteritems():
            if egnum == 0:
                delta -= modval_from_corres.get(corres, 0)
                new_modval_from_corres[corres] = None
                continue

            e = egnum / all_egnum
            a = [None for _ in range(partnum)]
            a_inv = [None for _ in range(partnum)]
            for part, com in enumerate(corres):
                egnum_of_com = new_egnum_from_com[part].get(
                                com, status.com.egnum_from_com(part, com) )

                a[part] = egnum_of_com / all_egnum
                a_inv[part] = 1 / a[part]

            a_multiplied = reduce(lambda x, y: x*y, a)
            alpha = e * sum(a_inv) * (1 / partnum)
            _mod = e  - a_multiplied
            partial_mod = _mod * alpha

            delta += partial_mod
            delta -= modval_from_corres.get(corres, 0)
            new_modval_from_corres[corres] = partial_mod

        # to prevent precision loss (桁落ち回避)
        if -0.00000009 < delta < 0.00000009:
            delta = 0

        # may be used later to reflect the moves
        modval_diff_info = {'delta_mod' : delta, 
                            'new_modval_from_corres' : new_modval_from_corres}
        
        return delta, modval_diff_info
        
    def update_modval_with_diff_info(self, modval_diff_info):
        delta_modval = modval_diff_info['delta_mod']
        new_modval_from_corres = modval_diff_info['new_modval_from_corres']

        # reflect correspondence information
        modval_from_corres = self._modval_from_corres
        for corres, modval in new_modval_from_corres.iteritems():
            if modval is None:
                del modval_from_corres[corres]
            else:
                modval_from_corres[corres] = modval

        # reflect modularity value
        self._modval += delta_modval

class ThresholdModularity(_AbstractModularity):

    def _calculate(self, status):
        partnum = status.basic.get_partnum()
        total_egnum = status.basic.get_edge_num()
        mod = 0
        for corres, egnum_in_corres in status.com.iter_corres_egnum():
            e = egnum_in_corres / total_egnum

            a = [None for _ in range(partnum)]
            a_inv = [None for _ in range(partnum)]
            for part, com in enumerate(corres):
                a[part] = status.com.egnum_from_com(part, com) / total_egnum
                a_inv[part] = 1 / a[part]

            a_multiplied = reduce(lambda x, y: x*y, a)
            alpha = e * sum(a_inv) * (1 / partnum)
            _mod = e  - a_multiplied
            partial_mod = _mod * alpha

            if _mod > 0:
                mod += partial_mod
            
        return mod

class PowerModularity(_AbstractModularity):

    def __init__(self, power=2):
        _AbstractModularity.__init__(self)
        self.power = power

    def name(self):
        return self.__class__.__name__ + '(%2.3f)' % self.power

    def _calculate(self, status):
        partnum = status.basic.get_partnum()
        total_egnum = status.basic.get_edge_num()
        mod = 0
        for corres, egnum_in_corres in status.com.iter_corres_egnum():
            e = egnum_in_corres / total_egnum

            a = [None for _ in range(partnum)]
            a_inv = [None for _ in range(partnum)]
            for part, com in enumerate(corres):
                a[part] = status.com.egnum_from_com(part, com) / total_egnum
                a_inv[part] = 1 / a[part]

            a_multiplied = reduce(lambda x, y: x*y, a)
            alpha = e * sum(a_inv) * (1 / partnum)
            _mod = e  - a_multiplied
            partial_mod = _mod * alpha

            if _mod > 0:
                mod += pow(partial_mod, self.power)
            
        return mod


#--------------------------------------------------
# test
#--------------------------------------------------
def _test_modularity(status, mod_list, answer_list):
    for ind, mod in enumerate(mod_list):
        modname = mod.name()
        modval = mod.calculate(status)
        ans = answer_list[ind]
        tn = modval == ans
        print "    %-40s : %.5f (ans: %.5f), %s" % (modname, modval, ans, tn)
    print ''

def _test_delta_calculation(modularity, moved_vrts_info, answer):
    moving_diff = status.com.diff_of_moving_vrts(moved_vrts_info)
    delta, modval_diff = modularity.calculate_diff(status, moving_diff)
    modval = modularity.value() + delta
    ans = answer
    tn = modval == answer
    print "    mod: %.5f (del: %.5f), ans: %.5f, %s" % (modval, delta, ans, tn)
    status.com.update_com_with_diff_info(moving_diff)
    modularity.update_modval_with_diff_info(modval_diff)

if __name__ == "__main__":
    from _status import NetworkStatus
    edge_list = [[0, 0, 0], [0, 0, 1], [1, 1, 0], [1, 1, 1],
                 [2, 2, 2], [2, 2, 3], [3, 3, 2], [3, 3, 3]]
    status = NetworkStatus(edge_list)
    status.add_com()

    ##############################################
    # TEST: modularity calculation
    ##############################################
    print 'Modularity calcuation test'
    mod_list = [NeubauerModularity(),
                MurataModularity(),
                ThresholdModularity(),
                PowerModularity(power=2)]

    # modularity calculation, test 1
    print '  All vertices are same community'
    com_labels = [[0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0]]
    answer_list = [0.0, 0.0, 0.0, 0.0]
    status.com.set_com_labels(com_labels)
    _test_modularity(status, mod_list, answer_list)

    # modularity calculation, test 2
    print '  Idealistic community structure for this network (probably)'
    com_labels = [[0, 0, 1, 1],
                       [0, 0, 1, 1],
                       [0, 0, 1, 1]]
    answer_list = [0.75, 0.75, 0.75, 0.28125]
    status.com.set_com_labels(com_labels)
    _test_modularity(status, mod_list, answer_list)

    # modularity calculation, test 3
    print '  Asymetric community structure'
    com_labels = [[0, 1, 1, 1],
                       [0, 1, 1, 1],
                       [0, 0, 0, 0]]
    answer_list = [0.3125, 0.3125, 0.3125, 0.04931640625]
    status.com.set_com_labels(com_labels)
    _test_modularity(status, mod_list, answer_list)

    #exit()

    ##############################################
    # TEST: delta calculation
    ##############################################
    print 'Delta calcuation test'
    com_labels = [[0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0]]
    modularity = NeubauerModularity()
    status.com.set_com_labels(com_labels)
    modval = modularity.calculate(status)
    print "    %.5f (first modularity value)" % modval

    # delta calculation, test 1
    moved_vrts_info = ({2 : [0, 1],
                        3 : [0, 1]},
                       {2 : [0, 1],
                        3 : [0, 1]},
                       {2 : [0, 1],
                        3 : [0, 1]})
    answer = 0.75
    _test_delta_calculation(modularity, moved_vrts_info, answer)

    # delta calculation, test 2
    moved_vrts_info = ({1 : [0, 1]},
                       {1 : [0, 1]},
                       {2 : [1, 0],
                        3 : [1, 0]})
    answer = 0.3125
    _test_delta_calculation(modularity, moved_vrts_info, answer)

    # delta calculation, test 3
    moved_vrts_info = ({0 : [0, 1]},
                       {0 : [0, 1]},
                       {})
    answer = 0.0
    _test_delta_calculation(modularity, moved_vrts_info, answer)
  
