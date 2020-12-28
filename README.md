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
