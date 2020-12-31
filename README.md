# rna_quant

squant reads a cs423 file and estimates the true count for each of the transcripts provided in the file header, based on approximations on the read alignments in the file alignment blocks.

## Usage
squant runs as the project description specifies, with two mandatory flags("--in", "--out"). The order of the mandatory flags may be switched. 

squant uses an external JAR library from apache commons (commons-math3).

The following commands should be run inside the terminal after locating to the project directory.

The following command builds the project.

javac -cp ".:./commons-math3-3.6.1/commons-math3-3.6.1.jar" squant.java

The following commands execute the standard full EM algorithm. The specific parameters should be changed to the name of actual input files (with full directory path) if they differ.

java -cp ".:./commons-math3-3.6.1/commons-math3-3.6.1.jar" squant --in /data/alignments.cs423 --out results.tsv
or equivalently
java -cp ".:./commons-math3-3.6.1/commons-math3-3.6.1.jar" squant --out results.tsv --in data/alignments.cs423

# squant
This project is about simulating the quantification of RNA transcripts as RNA sequencing would. 
We were provided with a set of RNA transcripts and their alignments to the original DNA sequence. 

For the transcripts, we were given with a name for the transcript and its length.

For alignments, we were given with a name of the transcript in question, its orientation, the position of the alignment, and the probability for alignment.

We used the Expectation Maximization algorithm to achieve this and followed a provided model with a set of computation steps.

| Method      | Description |
| ----------- | ----------- |
| fullEM()     | The full EM algorithm. Its core function follows the implementation logic provided by the hackmd.io page. It is also responsible for calling storeEffL() which calculates and stores the effective lengths of the transcripts, setUpReadsP() which evaluates and stores the p_2, p_3, and p_4 values for all reads, and initialEStep() which performs the first iteration of the EM algorithm. It uses converged() to determine if the EM algorithm is completed.      |
| converged()   | fullEM()’s utility method. Determines if the EM algorithm has converged by comparing all current p_1 values with the p_1 values from the previous iteration. It calls equalUpTo() to determine convergence between two values.        |
| equalUpTo()    | converged()’s utility method. Determines if two doubles can be considered equal by being equal up to a certain decimal place. |
| initialEStep() | fullEM()’s utility method. Carries out the first E-step of the EM algorithm. Evaluates the probability count of each transcript appearing using p_1 =1/M. |
| storeEffL()    | fullEM()'s utility method. Calculates the effective length for all transcripts by calling effectiveLength(). |
| setUpReadsP()  | fullEM()'s utility method. Sets up the probability values for all reads. |
| effectiveLength()| storeEffL()’s utility method. Computes the effective lengths for all transcripts.|
| writeFile()    | Utility method for writing the resulting estimates and lengths for all transcripts.|
| readFile()     | Utility method for reading and processing all data from the input file.|
