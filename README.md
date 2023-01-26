# Single Source Personalized PegeRank
This repository is for "Cache-Efficient Approach for Index-Free Personalized PageRank"

## Citation
```
@article{tsuchida2023cache,
  author={Tsuchida, Kohei and Matsumoto, Naoki and Shin, Andrew and Kaneko, Kunitake},
  journal={IEEE Access}, 
  title={Cache-Efficient Approach for Index-Free Personalized PageRank}, 
  year={2023},
  volume={11},
  pages={6944-6957},
  doi={10.1109/ACCESS.2023.3237738}}
```

## Environmnets
- C++ 17
- G++ 9.4.0
- Ubuntu 20.04
- SFMT 1.5.1 (SIMD-Oriented Fast Mersenne Twister)

## Compile
```
$ make
```

## Usage


```
$ ./SSPPR -algo <algorithm> -graph <graph-path> -attribute <attribute-path> -query <query-path>
                            -query_size <must be greater than 1> [-alpha <must be between 0 and 1>] [-epsilon <must be between 0 and 1>]
```

- algorithm
  - ForwardPush
  - MonteCarlo
  - [Fora](https://dl.acm.org/doi/abs/10.1145/3097983.3098072) (KDD'17)
  - [ResAcc](https://ieeexplore.ieee.org/abstract/document/9101811) (ICDE'20)
  - [SpeedPPR](https://dl.acm.org/doi/abs/10.1145/3448016.3457298) (SIGMOD'21)
  - CachePPR (**the proposed method**)
- parameters
  - graph: graph file
  - attribute: attribute file
  - query: query file
  - query_size: this value must be greater than one 
  - alpha: a termination probability of random walk (default 0.2)
  - epsilon: an error bound (default 0.5)

### PPR computation of the DBLP dataset (Example)
```
$ ./SSPPR -algo CachePPR -graph data/dblp/graph.txt -attribute data/dblp/graph.attribute -query data/dblp/graph.query -query_size 10
```

## Graph Format
Input grpah should follow [SNAP](http://snap.stanford.edu/data/) format.
Attribute file contains the number of nodes and edges respectively.

You can see the example in ```./data/dblp/``` folder.



※ We only accept directed graphs. If necesarry, you can convert undirected graphs -> directed graphs by ```ToDirected``` algorithm.
```
./SSPPR -algo ToDirected -graph <graph-path> -attributed <attributed-path>
```

## Query Generation
Query file has the list of query nodes of PPR computation.
You can generate query nodes by ```QueryGenerate``` algorithm.
```
./SSPPR -algo QueryGenerate -graph <graph-path> -attributed <attributed-path> -query <query-path> -query_size <query_size>
```
