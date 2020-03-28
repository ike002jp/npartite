# npartite
　  
**このパッケージは，n部ネットワークからのコミュニティ抽出に関連するアルゴリズムを含んでいます．**  

**This package contains some algorithms about community detection from n-partite networks.**  
([English is here](https://github.com/ike002jp/npartite/blob/master/README_en.md))

## 概要
2部ネットワークや3部ネットワークなどのn部ネットワークから，
コミュニティを抽出するためのアルゴリズムを含んでいます．
コミュニティ抽出アルゴリズムはモジュラリティ値の最適化に基づくものです．

コミュニティ抽出アルゴリズムだけではなく，
人工的なn部ネットワークの生成アルゴリズム，コミュニティ抽出結果の評価指標なども含んでいます．
各種アルゴリズムを用いて，人工ネットワークや実ネットワーク上で様々な実験を行うことが可能です．

以下の画像は，3部ネットワークからのコミュニティ抽出の例です．
![3部ネットワークからのコミュニティ抽出](https://raw2.github.com/ike002jp/npartite/master/community.png  "3部ネットワークからのコミュニティ抽出")

## パッケージの内容

* 4種類のn部モジュラリティ
    * Murataのモジュラリティ \[[Murata 2010](#Murata2010)\]
    * Neubauerのモジュラリティ \[[Neubauer 2010](#Neubauer2010)\]
    * 閾値モジュラリティ \[[Ikematsu 2014](#Ikematsu2014)\]
    * べき乗モジュラリティ \[[Ikematsu 2014](#Ikematsu2014)\]

* 3種類のモジュラリティ最適化手法
    * ノードをボトムアップにクラスタリングする手法
    * エッジををボトムアップにクラスタリングする手法
    * FUE(Fast Unfolding for Edges) \[[Ikematsu 2013](#Ikematsu2013)\]
        * FUEは2種類あります．1つはどのようなn部ネットワークにも適用可能です．
          もう1つは3部ネットワークのみにしか適用できませんが，高速です．

* コミュニティ抽出結果の評価指標
    * NMI（正規化相互情報量）\[[Ana 2003](#Ana2003)\]

* 人工的なn部ネットワークの生成アルゴリズム
    * （ほぼ）任意の正解構造を持つn部ネットワークを定義して，生成することが可能
    * 生成したn部ネットワークはエッジリストとして利用可能
    * 3部ネットワークについては，既に6種類の人工ネットワークを定義済み \[[Ikematsu 2014](#Ikematsu2014)\]
        * SIMPLE ケース
        * OVERLAPPING ケース
        * CONTRADICTIVE ケース
        * SIMPLE 3 ケース
        * SIMPLE PLUS ケース
        * SIMPLE 3 PLUS ケース
    * 2部ネットワークは，3部ネットワークにおけるSIMPLE ケースと同様の構造のみ定義済み

## 必要なモジュール
NMI値を計算するためには[numpy](http://www.numpy.org/)が必要です．

## 使い方
以下の例では，SIMPLE ケースのような正解構造を持つ人工3部ネットワークを作成し，
それに対してコミュニティ抽出を行っています．
また得られたコミュニティ抽出結果を，NMIを用いて正解構造と比較しています．

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

## 参考文献
* <a name="Ana2003"></a> \[Ana 2003\] Ana, L., Jain, A.: Robust data clustering. In: Proceedings of 2003 IEEE Computer Society Conference on Computer Vision and Pattern Recognition, vol. 2, pp. II–128–II–133 (2003)
* <a name="Murata2010"></a> \[Murata 2010\] Murata, T.: Modularity for Heterogeneous Networks, in Proceedings of the 21st ACM Conference on Hypertext and Hypermedia, pp. 129-134 (2010)
* <a name="Neubauer2010"></a> \[Neubauer 2010\] Neubauer, N., Obermayer, K.: Community detection in tagging-induced hypergraphs. In: Workshop on Information in Networks (2010)
* <a name="Ikematsu2013"></a> \[Ikematsu 2013\] Kyohei Ikematsu, Tsuyoshi Murata. A Fast Method for Detecting Communities from Tripartite Networks.  In Proceedings of the 5th International Conference on Social Informatics (SocInfo2013), pp.192-205, November 2013, Kyoto Japan
* <a name="Ikematsu2014"></a> \[Ikematsu 2014\] 池松恭平，村田剛志．3部モジュラリティの改善とその最適化手法．人工知能学会論文誌，Vol.29，No.2，2014，245-258


