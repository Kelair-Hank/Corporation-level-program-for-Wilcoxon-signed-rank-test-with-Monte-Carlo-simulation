# Enterprise-level-program-for-Wilcoxon-signed-rank-test-with-Monte-Carlo-simulation

Introduction<br />
This is a Enterprise level program for implementing Wilcoxon signed rank test with 1000 times resampling Monte Carlo simulation.<br />
<br />
Features:<br />
This program tests genomic sequence alignment accuracy of multiple new algorithms against an original algorithm. The accuracy is determined by mean correct rate of each subread of 3rd generation genomic sequencing technology (nanopore). The results generated by the original algorithm and each new algorithm are stored in different directories.<br />
<br />
Implementation:<br />
1. Read through all the notes in the program: "Monte Carlo Simulation with Wilcoxon signed rank test.jl"<br />
2. Change directory names to those directories each algorithm running results are located.<br />
3. Change file names to the file you record mean correct rate of each subread.<br />
4. Change algorithm names to the algorithms you are testing.
5. Run the program at the directory 1 layer above those directories each algorithm running results are located with the following code:<br />
       include("Monte Carlo Simulation with Wilcoxon signed rank test.jl")
