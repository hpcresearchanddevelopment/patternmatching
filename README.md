Distribute Pattern Matching on Large Metadata Graphs

<img src="https://github.com/hpcresearchanddevelopment/patternmatching/blob/master/examples/doc/tree_0011.png" width="200" height="200">

git@github.com:hpcresearchanddevelopment/patternmatching.git

srun -N1 --ntasks-per-node=4 --distribution=block ./src/generate_rmat -s 18 -p 1 -f 1 -o /dev/shm/rmat -b /urs/graph/rmat

cd  patternmatching/build/quartz/

srun -N1 --ntasks-per-node=4 --distribution=block ./src/ run_pattern_matching_beta -i /dev/shm/rmat -b /usr/graph/rmat -p ../../examples/rmat_log2_tree_pattern/ -o ../../examples/results/
