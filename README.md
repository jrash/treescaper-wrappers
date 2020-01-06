# treescaper-wrappers

Collection of scripts for analyzing sets of trees with treescaper.  These scripts were used to conduct the analyses described in this [poster](https://www.math.fsu.edu/~whuang2/pdf/JRAposterEvolution2014Final.pdf) presentation given at the Evolution conference in Raleigh, NC on June, 2014. 

See Genevieve Mount's updates to this workflow [here](https://github.com/jrash/TreeScaper_wrappers_2017).

## SeqSimRep.sh

Makes a folder for each sequence length and branch length pair specified.  Makes a folder for each rep.  Simulates sequences with conflicting signal and then does a mr bayes and garli bootstrap analysis.

Usage: SeqSimRep.sh [sequence lengths] [scale of branch lengths] [rep number from] [rep number to]

Example: ./SeqSimRep.sh '1000 5000' '.1 1' 1 5
 * Produces 4 folders - sequence length 1000 with .1 and 1 branch length scales and 5000 with .1 and 1 branch length scales
 * Each folder contains replicates 1-5

Input files:
 * Tree[1/2].phy: the guide trees used to simulate sequence alignments
 
Output files:
 * SimSeqs/Tree[1/2].seqs.nex: the sequence alignments simulated from either guide tree
 * SimSeqs/Tree.cat.seqs.nex: Both sequence alignments concatenated together

 * Bayes/Tree.cat.seqs.nex.burn.t: the set of posterior trees from the 5 runs, with burnin removed
 * Bayes/Tree.cat.seqs.nex.con.tre: the consensus tree of the above tree set
 * Bayes/Tree.con.Garli.${1}bp${2}L.pdf: pdf of the consensus tree

 * MLboot/Tree.cat.seqs.nex.boot.tre:  The bootstrap tree set from the garli bootstrap analysis
 * MLboot/Tree.cat.seqs.nex.con: the consensus tree of the above tree set
 * MLboot/Tree.con.Garli.${1}bp${2}L.pdf: pdf of the consensus tree

 * MLbest/Tree.cat.seqs.nex.best.tre: garli maximum likelihood tree

## treescaperWrapperV2.py

Runs community detection on the tree set with range of lambda values specified. Finds the plateau of communities detected and prints to ouput files. The input tree set needs to be a newick file named "all_trees.new".

Usage: treescaperWrapperV2.py [model] [interval] [network]
 * [model] can be CNM/CPM/ERNM/NNM
 * [interval] size of the lambda intervals in the grid search
 * [network] can be Covariance or Affinity

Output files:

 * ['treeset']_['type']WholeCommunity_results.out: community results over the whole range of lambda values
  
   * ['type'] can be Covariance or Affinity
   * ['treeset'] tree set name
  
 * ['treeset']_CovPlateauCommunity.out: community structure of the plateau
 * ['treeset']_comKey.out: key showing you which bipartitions are in which communities
 * AffinityCom[number].nex: a nexus file of the trees in an affinity community
 * AffinityCom[number].nex.con: consensus tree of an affinity community
 * AffinityCom[number].nex.con.pdf: pdf of consensus tree of an affinity community


## treescaperWrapperKnownPlateau.py

Useful if you have found the plateau with the automatic search function of the treescaper GUI.  If you enter the lambda values where the plateau was found for both affinity and covariance matrices, you will get all the useful output files of treescaperWrapperV2.py.

Usage: treescaperWrapperKnownPlateau.py [model] [plateau] [network] [rooted]
 * [model] can be CNM/CPM/ERNM/NNM
 * [plateau] is a numeric lambda value
 * [network] can be Covariance or Affinity
 * [rooted] is binary

Output files: same as treescaperWrapperV2.py

## AffinityCommunities.py 

Usage: AffinityCommunities.py 'Path/To/Treeset' [model] [plateau]
 * [model] can be CNM/CPM/ERNM/NNM
 * [plateau] is a numeric lambda value

Output files:

* AffinityCom[number].nex: a nexus file of the trees in an affinity community
* AffinityCom[number].nex.con: consensus tree of an affinity community
* AffinityCom[number].nex.con.pdf: pdf of consensus tree of an affinity community
* AffinityCommunitiesTreeCount.txt: a file to count the frequency and relative frequency of trees in each affinity community

# Example workflow

Below is a typical workflow to run in a bash shell. This was used to generate the output files in the repository.

```
CLVTreeScaper -trees -f all_trees.nex -ft Trees -w 0 -r 0 -o Community -t Covariance -cm CPM -lm auto -hf .95 -lf .05 > all_trees_CovAuto.out
CLVTreeScaper -trees -f all_trees.nex  -ft Trees -w 0 -r 0 -o BipartMatrix -bfm matrix 
CLVTreeScaper -trees -f all_trees.nex -ft Trees -w 0 -r 0 -o Community -t Affinity -cm CPM -lm auto -dm URF -am Exp > all_trees_AffAuto.out

treescaperWrapperKnownPlateau.py CPM .06 Covariance 1
treescaperWrapperKnownPlateau.py CPM .06 Affinity 1
```

