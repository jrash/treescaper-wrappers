#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Usage: treescaperWrapperKnownPlateau.py [model] [plateau] [network] [rooted]
#[model] can be CNM/CPM/ERNM/NNM

#***Check the ./CLVTreescaper settings within the script.  You might want to run ./CLVTreescaper with different options.  Make sure the
#numbering of all the indices is correct. The developers have changed between starting at 1 and starting at 0.  This includes
#the affinityCommunities script that is imported into this script.

#useful if you have found the plateau with the automatic search function of the treescaper GUI.  If you enter the lambda values where the plateau was found for
#both affinity and covariance matrices, you will get all the useful output of treescaperWrapperV2.py.  See treescaperWrapperV2.py for usage and output.

#output files

# ['type'] can be Covariance or Affinity
# ['treeset'] tree set name
# ['plateauLambda'] the lambda value were the plateau was found

#['treeset']_['type']WholeCommunity_results.out: community results over the whole range of lambda values
#['treeset']_CovPlateauCommunity.out: community structure of the plateau
#['treeset']_comKey.out: key showing you which bipartitions are in which communities
# AffinityCom[number].nex: a nexus file of the trees in an affinity community
# AffinityCom[number].nex.con: consensus tree of an affinity community
# AffinityCom[number].nex.con.pdf: pdf of consensus tree of an affinity community


###Numbering looks right now####

import re
import os
import sys
import numpy as np
import fnmatch
from collections import Counter
import filecmp
from AffinityCommunities import affinityCommunityConsensus
import dendropy
import StringIO
import heapq
from shutil import move, copyfile
from tempfile import mkstemp


def make_list(ls,convert):
	ls = ls.split("\t")
	ls.remove("\n")
	if convert == "int":
		ls = [int(i) for i in ls]
	if convert == "float":
		ls = [float(i) for i in ls]
# 	ls = ls[1:]
	return ls


def reg_ex_match(file, pattern):
#returns the first match of a reg ex search

	file.seek(0)
	for line in file:
		m = pattern.match(line)
		if m:
			return m.group(1)


def mode_function(lst):
#Returns a list of all the possible plateaus


	counterLst = Counter(lst)
	_,val = counterLst.most_common(1)[0]
	modeLS = [x for x,y in counterLst.items() if y == val]
	lstNoMode = [x for x in lst if x not in modeLS] #remove all the values that are the mode of the list
	if len(counterLst.most_common(2)) == 2: # find the next mode of the list
		counterLstNoMode = Counter(lstNoMode)
		_,val2 = counterLstNoMode.most_common(1)[0]
		if val - val2 >= 2:
			#if the largest plateau is bigger than the second largest by one increment, than its possible that the true size of the second largest
			#plateau is bigger than the largest. The second largest plateau could encompass the space between the lambda value immediately below and immediately above its lambda range.  The second largest plateau could
			#grow by 1.9999... increment and largest plateau could not grow at all. But if the largest plateau is bigger than the second largest by two or more increments, you are sure that plateau is the largest
			#****newly added feature, not used for any of my results
			return modeLS
		else:
			return [x for x,y in counterLst.items() if y == val or y == val2]
	else:
		return modeLS
	
	
def get_plateau(treeSet, treeSetTrunc, type, model, rooted):


	if type == "Covariance":
		comResults = open("%s_CovWholeCommunity_results.out" % treeSetTrunc)
	if type == "Affinity":
		comResults = open("%s_AffWholeCommunity_results.out" % treeSetTrunc)
	labelLS = comResults.readline()
	labelLS = make_list(labelLS, "int")
	print labelLS

	lambdaLS = comResults.readline()
	lambdaLS = make_list(lambdaLS,"float")
	print lambdaLS

# 	comNumLS = comResults.readline()
# 	conNumLS = make_list(comNumLS, "int")
# 	print conNumLS
# 	m = max(labelLS)
#
# 	labelLStrim= [i for i in labelLS if i != m and i != 0]
# 	print labelLStrim
#
# 	if type == "Affinity" and labelLStrim != []:
#
# 		m = max(labelLStrim)
# 		print m
#
# 		labelLStrim= [i for i in labelLStrim if i != m]
# 		print labelLStrim
#
#
# 	if labelLStrim != []:  #if using traditional search for plateau, call mode_function
# 		plateauLabel = mode_function(labelLStrim)
# 	else:  #if using automatic search and feeding the results, use the only result in the list
# 		plateauLabel = []
# 		plateauLabel.append(labelLS[1])
# 	print plateauLabel
#
#
# 	if len(plateauLabel) == 1:
#
# 		plateauIndex = labelLS.index(plateauLabel[0])

	plateauLambda = lambdaLS[0]
	print plateauLambda

	if type == "Covariance":
		os.system("/home/vestige/Documents/BrownLabJash/Programs/bin/treescaper_scripts_2017/CLVTreeScaper -trees "+\
		"-f %s -ft Trees -w 0 -r %s -o Community -t Covariance -cm %s -lm manu -lp %s -ln 1 -hf .95 -lf .05" % (treeSet, rooted, model, plateauLambda)+\
		" > %s_CovPlateauCommunity.out" %  treeSetTrunc)


		comFile = open("%s_CovPlateauCommunity.out" %  treeSetTrunc, 'r')
		pattern = re.compile('Number of communities: (\d+)')
		coms = int(reg_ex_match(comFile, pattern))
		comKey = open("%s_comKey.out" %  treeSetTrunc, 'w')

		for i in range(1, coms+1): # makes a key to decode the bipartitions in each community
			pattern = re.compile('Community '+str(i)+' includes nodes: (.+)')
			comStr = reg_ex_match(comFile, pattern)
			comLS = comStr.split(",")
			print comLS
			comLS = filter(None, comLS)
			comKey.write("Com %s:\n" % i)
			for j in comLS:
				j = int(j)+1
				pattern = re.compile("bipartition "+str(j)+" : ([0-1]+?), appear times: ([0-9]+?)$")
				print "bipartition "+str(j)+" : ([0-1]+?), appear times: ([0-9]+)"
				fh, absPath = mkstemp()
				copyfile("%s_CovPlateauCommunity.out" %  treeSetTrunc, absPath)
				comTempFile = open(absPath,'r')
				for line in comTempFile:
					m = pattern.match(line)
					if m:
						print m
						bipart = m.group(1)
#						print bipart
						freq = m.group(2)
						bipartLS = []
						for c in bipart:
							bipartLS.append(c)
#						print bipartLS
						indices = [x+1 for x, y in enumerate(bipartLS) if y == '1']
						for k in indices:
							pattern2 = re.compile("(.+) , "+str(k))
							taxon = reg_ex_match(comFile, pattern2)
							indices[indices.index(k)] = taxon
				 #close temp file
				comTempFile.close()
				os.close(fh)

				comKey.write("%s %s\n" % (indices,freq))
			comKey.write("\n")
		comKey.close()










	if type == "Affinity":
		os.system("/home/vestige/Documents/BrownLabJash/Programs/bin/treescaper_scripts_2017/CLVTreeScaper -trees "+\
		"-f %s -w 0 -r %s -o Community -t Affinity -cm %s -lm manu -dm URF -am Exp -lp %s -ln 1 " % (treeSet, rooted, model, plateauLambda)+\
		" > %s_AffPlateauCommunity.out" %  treeSetTrunc)#outputs community structure for current lambda values


	return plateauLambda




def main():
	model = sys.argv[1]
	plateau = sys.argv[2]
	network = sys.argv[3]
	rooted = sys.argv[4]
	for file in os.listdir('.'):
		for file in os.listdir('.'):
			if fnmatch.fnmatch(file, 'all_trees.nex'):
				if rooted == '1':
					tlst = dendropy.TreeList.get_from_path("all_trees.nex", "nexus", rooting='force-rooted')
				else:
					tlst = dendropy.TreeList.get_from_path("all_trees.nex", "nexus", rooting='force-unrooted')
			elif fnmatch.fnmatch(file, 'all_trees.new'):
				if rooted == '1':
					tlst = dendropy.TreeList.get_from_path("all_trees.new", "newick", rooting='force-rooted')
				else:
					tlst = dendropy.TreeList.get_from_path("all_trees.new", "newick", rooting='force-unrooted')

	tlst.write_to_path("all_trees.pre.nex",'nexus', simple=True, translate_tree_taxa=True)

	with open("all_trees.nex", "w") as fout:
		with open("all_trees.pre.nex", "r") as fin:
			for line in fin:
				fout.write(re.sub('END;\n', 'END;',line))

	os.system("rm all_trees.pre.nex")
	treeSet="all_trees.nex"
	treeSetIndex = treeSet.find(".")
	treeSetTrunc = treeSet[:treeSetIndex]
	os.system("echo 'hi'")
	
	if network == 'Covariance':

	 	os.system("/home/vestige/Documents/BrownLabJash/Programs/bin/treescaper_scripts_2017/CLVTreeScaper -trees "+\
	 	"-f %s -ft Trees -w 0 -r %s -o Community -t Covariance -cm %s -lm manu -lp %s -ln 1 -hf .95 -lf .05" % (treeSet, rooted, model, plateau)+\
	 	" > %s_CovCommunity.out" %  treeSet)#outputs community structure for current lambda values
	 	os.system("mv %s_Covariance\ Matrix_community_results.out %s_CovWholeCommunity_results.out" % (treeSetTrunc,treeSetTrunc))

	 	plateauLambda = get_plateau(treeSet, treeSetTrunc, "Covariance", model, rooted)

	if network == 'Affinity':

		os.system("/home/vestige/Documents/BrownLabJash/Programs/bin/treescaper_scripts_2017/CLVTreeScaper -trees "+\
		"-f %s -w 0 -r %s -o Community -t Affinity -cm %s -lm manu -dm URF -am Exp -lp %s -ln 0 " % (treeSet, rooted, model, plateau)+\
		" > %s_AffCommunity.out" %  treeSet)#outputs community structure for current lambda values
		os.system("mv %s_Affinity-URF_community_results.out %s_AffWholeCommunity_results.out" % (treeSetTrunc, treeSetTrunc))

		plateauLambda = get_plateau(treeSet, treeSetTrunc, "Affinity", model, rooted)

		if plateauLambda:
			affinityCommunityConsensus(treeSet, model ,plateauLambda, rooted)

 	os.system("sumtrees.py -r -o all_trees.con all_trees.nex")
 	os.system("cat ~/Documents/BrownLabJash/Programs/bin/SeqSim/FigTreeBlock.txt >> all_trees.con")
 	os.system("figtree -graphic PDF all_trees.con all_trees.pdf")

	os.system("echo 'treescaperWrapperKnownPlateau.py %s\n' >> commands.txt" % ' '.join(sys.argv[1:]))


	#os.system("cat RAxML_bestTree.G1_gene* > RAxML_allTree.G1.tre")
if __name__=='__main__':
	main()
