# SOWalker
an novel second-order graph processing system for random walk

```
./bin/test/node2vec <dataset> [weighted] [sorted] [skip] [blocksize] [nthreads] [dynamic] [cache_size] [max_iter] [walkpersource] [length] [p] [q]

- dataset:       the dataset path
- weighted:      whether the dataset is weighted
- sorted:        whether the vertex neighbors is sorted
- skip:          whether to skip the interactive preprocess query
- blocksize:     the size of each block
- nthreads:      the number of threads to walk
- dynamic:       whether the blocksize is dynamic, according to the number of walks
- cache_size:    the size(GB) of cache
- max_iter:      the maximum number of iteration for simulated annealing scheduler
- walkpersource: the number of walks for each vertex
- length:        the number of step for each walk
- p:             node2vec parameter
- q:             node2vec parameter
```

## Node2vec

```
./bin/test/node2vec /dataset/livejournal/w-soc-livejournal.txt  length 20 walkpersource 1
```

