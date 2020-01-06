# treescaper-wrappers

Collection of scripts for analyzing sets of trees with treescaper.  These scripts were used to conduct the analyses described in this [poster](https://www.math.fsu.edu/~whuang2/pdf/JRAposterEvolution2014Final.pdf) presentation given at the Evolution conference in Raleigh, NC on June, 2014. 

## treescaperWrapperKnownPlateau.py

Usage: treescaperWrapperKnownPlateau.py [model] [plateau] [network] [rooted]
 * [model] can be CNM/CPM/ERNM/NNM
 * [plateau] is a numeric lambda value
 * ['network'] can be Covariance or Affinity
 * [rooted] is binary

Useful if you have found the plateau with the automatic search function of the treescaper GUI.  If you enter the lambda values where the plateau was found for
both affinity and covariance matrices, you will get all the useful output of treescaperWrapperV2.py.  See treescaperWrapperV2.py for usage and output.

output files

 * ['treeset']_['type']WholeCommunity_results.out: community results over the whole range of lambda values
  
   * ['type'] can be Covariance or Affinity
   * ['treeset'] tree set name
  
 * ['treeset']_CovPlateauCommunity.out: community structure of the plateau
 * ['treeset']_comKey.out: key showing you which bipartitions are in which communities
 * AffinityCom[number].nex: a nexus file of the trees in an affinity community
 * AffinityCom[number].nex.con: consensus tree of an affinity community
 * AffinityCom[number].nex.con.pdf: pdf of consensus tree of an affinity community
