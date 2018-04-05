<p>Distribute Pattern Matching on Large Metadata Graphs</p>

<p>Here, we present an example of searching a pattern in a R-MAT generated graph using our program. The code is developed on top of HavoqGT.</p>

<p>Clone (with SSH) the code from git@github.com:hpcresearchanddevelopment/patternmatching.git:
<br/>
git clone git@github.com:hpcresearchanddevelopment/patternmatching.git</p>

<p>You will require the latest releases of OpenMPI or MAVPICH2 and the Boost library (some Boost releases have bugs, e.g., 1.58, the code works fine with 1.57) to run HavoqGT. The code has only been tested on latest generation of Linux distributions. Once you have checked out the code, make sure you are on the master branch.</p>

cd  patternmatching/build/quartz/

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
srun -N1 --ntasks-per-node=4 --distribution=block ./src/generate\_rmat -s 21 -p 1 -f 1 -o /dev/shm/rmat -b /urs/graph/rmat
</p>

<p>This will create a graph with four partitions, to be run of four MPI processes. Note that this is a Scale 21 graph. (Notice the parameter for the -s flag). The mmap/binary graph file will be store in /usr/graph/</p>

<div align="center"><img src="https://github.com/hpcresearchanddevelopment/patternmatching/blob/master/examples/doc/tree_0011.png" width="200" height="200"></div>

<p>We use degree information to create numeric vertex labels, computed using the formula. ceil(log\_2(d(v\_i)+1)). Here, d(v\_i) is the degree of a vertex v\_i.</p>

<p>The input pattern is available in the following dircetory: patternmatching/examples/rmat\_log2\_tree\_pattern/</p>

<p>The program requires a predefined directory structure to output results: patternmatching/examples/results/ </p>

<p>Next, use the following command to search the pattern stored in patternmatching/examples/rmat\_log2\_tree\_pattern/.
<br/>
<br/>
Note that we do not need to provide vertex labels for the Tree pattern as we will use labels based on vertex degree and the program will generate labels when no input label is provided, i.e., the -l flag is not set.
<br/>
srun -N1 --ntasks-per-node=4 --distribution=block ./src/ run\_pattern\_matching\_beta -i /dev/shm/rmat -b /usr/graph/rmat -p ../../examples/rmat\_log2\_tree\_pattern/ -o ../../examples/results/
</p>

<p>The program logs status information to the standard output so you know the current state of the execution.</p>

<p>See the instructions on the readme page regarding how to interpret the results and retrieve the pruned graph. You will find scripts (written in python) that will help you to parse the result files.
The directory /examples/results/0/all\_ranks\_active\_vertices\_count/ contains the number of active vertices at the end of each iteration. However the results are distributed, 
The last entry in the file indicate the final number of active vertiices. (Same for edges as well). 

To accumuate the results 

python ../../examples/scripts/total\_active\_count.py ../../examples/results/0/all\_ranks\_active\_vertices\_count/ > vertices\_count

The file ../../examples/results/0/result\_superstep contains the timing information

</p>

<h2>Enumeration</h2>
<p>For the Tree pattern in this example, the last step in the execution enumerates the pattern in the pruned graph. The input for enumeration is the last entry in the file examples/rmat_log2_tree_pattern/0/pattern_non_local_constraint. Here, the same NLCC code walks the full template with work aggregation turned off.</p>

