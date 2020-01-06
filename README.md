# treescaper-wrappers

See Genevieve Mount's updates to this workflow [here](https://github.com/jrash/TreeScaper_wrappers_2017).

Collection of scripts for analyzing sets of trees with treescaper.  These scripts were used to conduct the analyses described in this [poster](https://www.math.fsu.edu/~whuang2/pdf/JRAposterEvolution2014Final.pdf) presentation given at the Evolution conference in Raleigh, NC on June, 2014. 

## treescaperWrapperV2.py

Runs community detection on the tree set with range of lambda values specified. Finds the plateau of communities detected and prints to ouput files. The input tree set needs to be a newick file named "all_trees.new". The lambda ranges and increments may need to be changed for affinity and covariance community detection. These sections are marked in the main() function.

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

# Example worklow

Below is a typical workflow to run in a bash shell. This was used to generate the output files in the repository.

```
CLVTreeScaper -trees -f all_trees.nex -ft Trees -w 0 -r 0 -o Community -t Covariance -cm CPM -lm auto -hf .95 -lf .05 > all_trees_CovAuto.out
CLVTreeScaper -trees -f all_trees.nex  -ft Trees -w 0 -r 0 -o BipartMatrix -bfm matrix 
CLVTreeScaper -trees -f all_trees.nex -ft Trees -w 0 -r 0 -o Community -t Affinity -cm CPM -lm auto -dm URF -am Exp > all_trees_AffAuto.out

treescaperWrapperKnownPlateau.py CPM .06 Covariance 1
treescaperWrapperKnownPlateau.py CPM .06 Affinity 1
```

