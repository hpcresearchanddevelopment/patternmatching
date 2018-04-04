<p>Distribute Pattern Matching on Large Metadata Graphs</p>

<p>Here, we present an example of searching a pattern in a R-MAT generated graph using our program. The code is developed on top of HavoqGT.</p>

<p>Clone (with SSH) the code from git@github.com:hpcresearchanddevelopment/patternmatching.git:
<br/>
git clone git@github.com:hpcresearchanddevelopment/patternmatching.git</p>

<p>You will require the latest releases of OpenMPI or MAVPICH2 and the Boost library (some Boost releases have bugs, e.g., 1.58, the code works fine with 1.57) to run HavoqGT. The code has only been tested on latest generation of Linux distributions. Once you have checked out the code, make sure you are on the master branch.</p>

<p>Go to the directory, build/quartz/
Setup CMake environment by running the following script: 
<br/>
./scripts/quartz/do\_cmake.sh
<br/>
(Make necessary adjustments to the script for CMake to work within your environment.)</p>

<p>The next step is to generate a graph in HavoqGT format. Go to the directory, build/quartz/ and build the R-MAT generator:
<br/>
make generate\_rmat</p>

<p>Create a directory, e.g., /usr/graph/, to store the generated graph. Assuming you are in the Slurm environment, run the following command to generate a R-MAT graph:
<br/>
srun -N1 --ntasks-per-node=4 --distribution=block ./src/generate\_rmat -s 18 -p 1 -f 1 -o /dev/shm/rmat -b /urs/graph/rmat
</p>

<p>This will create a graph with four partitions, to be run of four MPI processes. Note that this is a Scale 18 graph. (Notice the parameter for the -s flag). The mmap/binary graph file will be store in /usr/graph/</p>

<div><center><img src="https://github.com/hpcresearchanddevelopment/patternmatching/blob/master/examples/doc/tree_0011.png" width="200" height="200"></center></div>


srun -N1 --ntasks-per-node=4 --distribution=block ./src/generate_rmat -s 21 -p 1 -f 1 -o /dev/shm/rmat -b /urs/graph/rmat

cd  patternmatching/build/quartz/

srun -N1 --ntasks-per-node=4 --distribution=block ./src/ run_pattern_matching_beta -i /dev/shm/rmat -b /usr/graph/rmat -p ../../examples/rmat_log2_tree_pattern/ -o ../../examples/results/
