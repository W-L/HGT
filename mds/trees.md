The alignments from the OMA standalone run were analyzed with prottest to find the best model which was then used by RAxML to build a tree. Using the script tree_building.sh 

prottest3 -i $i -AIC

raxmlHPC-PTHREADS-SSE3 -f a -m PROTGAMMA$model -s $i -T 7 -p 12345 -n $name -# 20 -x 1234

-f a to run mode a, which means bootstrap estimation and maximum-likelihood tree estimation
-m PROTGAMMA$model parses the best model from prottest to RAxML
-s is the alignment file again
-T for number of threads
-p, -x random seeds
-n output name
-#20 number of bootstraps

The trees are plotted with distances of the branches in blue and bootstrap values in red, using ape.

