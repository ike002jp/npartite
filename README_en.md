# npartite

**This package contains some algorithms for community detection from n-partite networks.**  

## Abstract
This contains some algorithms for community detection from n-partite networks, 
such as bi-partite networks or tri-partite networks.
The algorithms for community detection are based on optimizing modularity.

This contains some other algorithms, not only the algorithms for community detection.
For example, the algorithm for generating synthetic n-partite networks
and evaluation measure for the results of community detection are contained.
We can do some experiments on the synthetic networks and the real networks with these algorithms.

The following image is an example of community detection from a tri-partite network.
![community detection from tri-partite networks](https://raw2.github.com/ike002jp/npartite/master/community.png  "community detection from tri-partite networks")

## Contents

* 4 types of n-partite modularities
    * Murata's modularity \[[Murata 2010](#Murata2010)\]
    * Neubauer's modularity \[[Neubauer 2010](#Neubauer2010)\]
    * Threshold modularity \[[Ikematsu 2014](#Ikematsu2014)\]
    * Power modularity \[[Ikematsu 2014](#Ikematsu2014)\]

* 3 types of modularity optimization methods
    * Clustering the vertices based on greedy bottom up manner
    * Clustering the edges based on greedy bottom up manner
    * FUE(Fast Unfolding for Edges) \[[Ikematsu 2013](#Ikematsu2013)\]

* Evaluation measure for the results of community detection
    * NMI (Normalized Mutual Information) \[[Ana 2003](#Ana2003)\]

* Algorithms for generating synthetic n-partite networks
    * 6 types of tri-partite networks have already defined \[[Ikematsu 2014](#Ikematsu2014)\]
        * SIMPLE Case
        * OVERLAPPING Case
        * CONTRADICTIVE Case
        * SIMPLE 3 Case
        * SIMPLE PLUS Case
        * SIMPLE 3 PLUS Case

##Required modules
[numpy](http://www.numpy.org/) is needed if you want to calculate the value of NMI.

##How to use
The following example makes a synthetic tri-partite network which is SIMPLE case.
After that, it detects communities on the network.
And then, it compares the detected communities and the correct community with NMI.

```python

from npartite.extcom.evaluation import calculate_nmi
from npartite.extcom.modularity import NeubauerModularity
from npartite.extcom.optimization import GreedyVertexBottomUp
from npartite.synthetic.tripartite import SimpleCaseMaker

# make synthetic network
syn_maker = SimpleCaseMaker(corres_egnum=40, noise_ratio=0)
syn_network = syn_maker.make()

# get edge list and correct community structure
edge_list = syn_network.edge_list()
correct_labels = syn_network.correct_community_labels()

# initialize modularity and optimization method
modularity = NeubauerModularity()
optimization = GreedyVertexBottomUp()

# detect communities
results = optimization.start(modularity, edge_list)

# modularity value, 
# detected community labels of the each vertex
# all modularity values through the process of the optimization 
mod_value, detected_labels, mod_values = results

# print each value
print 'modularity value: %f' % (mod_value, )
print 'detected community labels: %s' % (detected_labels, )
print 'correct community labels : %s' % (correct_labels, )
print 'all modularity values through the process: %s' % (mod_values, )


# calculate the value of NMI
nmi = calculate_nmi(detected_labels, correct_labels)
print 'NMI: %f' % (nmi, )

```

##References
* <a name="Ana2003"></a> \[Ana 2003\] Ana, L., Jain, A.: Robust data clustering. In: Proceedings of 2003 IEEE Computer Society Conference on Computer Vision and Pattern Recognition, vol. 2, pp. II–128–II–133 (2003)
* <a name="Murata2010"></a> \[Murata 2010\] Murata, T.: Modularity for Heterogeneous Networks, in Proceedings of the 21st ACM Conference on Hypertext and Hypermedia, pp. 129-134 (2010)
* <a name="Neubauer2010"></a> \[Neubauer 2010\] Neubauer, N., Obermayer, K.: Community detection in tagging-induced hypergraphs. In: Workshop on Information in Networks (2010)
* <a name="Ikematsu2013"></a> \[Ikematsu 2013\] Kyohei Ikematsu, Tsuyoshi Murata. A Fast Method for Detecting Communities from Tripartite Networks.  In Proceedings of the 5th International Conference on Social Informatics (SocInfo2013), pp.192-205, November 2013, Kyoto Japan
* <a name="Ikematsu2014"></a> \[Ikematsu 2014\] 池松恭平，村田剛志．3部モジュラリティの改善とその最適化手法．人工知能学会論文誌，Vol.29，No.2，2014，245-258


