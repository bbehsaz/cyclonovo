#!/usr/bin/env python
import sys, string, random

def suffix(string):
	return string[1:]
def prefix(string):
	return string[:len(string)-1]

def constructDeBruijn(kmers):
	adjacencies = {}
	nodes2aminos = {}
	n = 0
	edges = []
	nodes = set()
	
	#To get the deburijn graph that actually nodes are the k-mer size ... only for plotting purposes!
	'''
	for kmer1 in kmers:
		nodes.add(kmer1)
		for kmer2 in kmers:
			if suffix(kmer1) == prefix(kmer2):
				edges.append((kmer1,kmer2))
				if kmer1 not in adjacencies:
					adjacencies[kmer1]= [kmer2]
				else:
					adjacencies[kmer1].append(kmer2)
	print len(nodes)
	print len(edges)
	'''
	for kmer in kmers:
		nodes.add(prefix(kmer))
		nodes.add(suffix(kmer))
		if prefix(kmer) in adjacencies:
			adjacencies[prefix(kmer)].append(suffix(kmer))
			edges.append((prefix(kmer),suffix(kmer)))
		else:
			adjacencies[prefix(kmer)] = [suffix(kmer)]
			edges.append((prefix(kmer),suffix(kmer)))

	return adjacencies, list(nodes), edges