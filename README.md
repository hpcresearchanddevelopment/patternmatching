<p>Distribute Pattern Matching on Large Metadata Graphs</p>

<p>Here, we present an example of searching a pattern in a R-MAT generated graph using our program. The code is developed on top of HavoqGT.</p>

<p>Clone (with SSH) the code from git@github.com:hpcresearchanddevelopment/patternmatching.git
git clone git@github.com:hpcresearchanddevelopment/patternmatching.git</p>

<p>You will require the latest releases of OpenMPI or MAVPICH2 and the Boost library (some Boost releases have bugs, e.g., 1.58, the code works fine with 1.57) to run HavoqGT. The code has only been tested on latest generation of Linux distributions. Once you have checked out the code, make sure you are on the master branch.</p>

<p>Go to the directory, build/quartz/
Setup CMake environment by running the following script: 
./scripts/quartz/do\_cmake.sh
(Make necessary adjustments to the script for CMake to work within your environment.)</p>

<img src="https://github.com/hpcresearchanddevelopment/patternmatching/blob/master/examples/doc/tree_0011.png" width="200" height="200">


srun -N1 --ntasks-per-node=4 --distribution=block ./src/generate_rmat -s 21 -p 1 -f 1 -o /dev/shm/rmat -b /urs/graph/rmat

cd  patternmatching/build/quartz/

srun -N1 --ntasks-per-node=4 --distribution=block ./src/ run_pattern_matching_beta -i /dev/shm/rmat -b /usr/graph/rmat -p ../../examples/rmat_log2_tree_pattern/ -o ../../examples/results/
