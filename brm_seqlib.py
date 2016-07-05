import os
import sys
sys.path.append('/home/bradon/scripts/BRM_LIB')

def build_blosum_62():
	"""Returns a dict of the BLOSUM62 matrix. Stop codons are coded as '*'"""
	DATA = {}
	with open('/home/bradon/scripts/BRM_LIB/BLOSUM62.lst','r') as f:
		for ln in f:
			ln = ln.replace('\n','')
			info = ln.split('\t')
			if not DATA.has_key(info[0]):
				DATA[info[0]] = {}
			DATA[info[0]][info[1]] = int(info[2])
	return DATA

def calc_protein_bit(seq,BLOSUM):
	"""Calculates the max bitscore for a protein sequence"""
	if not isinstance(seq,str):
		print '\n{0}\n'.format(seq)
		raise TypeError("calc_protein_bit() requires a string input")

	bit = 0
	seq = list(seq)
	for i in seq:
		bit += BLOSUM[i][i]
	return bit


def revcomp (seq):
	"Generates reverse complement"

	if not isinstance(seq,str):
		raise TypeError("revcomp() requires a string input")
	f_lst = ("G", "C", "A", "T", "g", "c", "a", "t")
	r_lst = ("C", "G", "T", "A", "c", "g", "t", "a")
	r_seq = []
	seq = list(seq.strip())
	for i in reversed(seq):
		for n in range(len(f_lst)):
			if f_lst[n] is i:
				r_seq.append(r_lst[n])
				break
	return ''.join(r_seq)

def gc_calc (seq):
	"Calculates GC content"

	if not isinstance(seq,str):
		raise TypeError("gc_calc() requires a string input")
	seq = seq.upper()
	gc = 0.0
	at = 0.0
	gc += seq.count('G')
	gc += seq.count('C')
	at += seq.count('A')
	at += seq.count('T')
	total = gc + at
	perc = round(gc / total, 4) * 100
	return perc

def fasta2dict(f, length_only=False,full_header=False):
	"""Generates a dictionary from a fasta file: dict[name] = sequence or dict[name]=seq_length"""
	import collections as cl

	DATA = {}
	if os.path.isdir(f):
		print "ERROR: cannot read a directory"
		return
	if length_only == False:
		with open(f,'r') as fasta:
			for ln in fasta:
				ln = ln.replace('\n','')
				if '>' in ln:
					if not full_header:
						name = ln.split()[0].replace('>','')
					else:
						name = ln.replace('>','')
					DATA[name] = cl.deque()
				else:
					DATA[name].append(ln)
		for i in DATA:
			DATA[i] = ''.join(DATA[i])
	else:
		with open(f,'r') as fasta:
			for ln in fasta:
				ln = ln.replace('\n','')
				if '>' in ln:
					if not full_header:
						name = ln.split()[0].replace('>','')
					else:
						name = ln.replace('>','')
					DATA[name] =0
				else:
					DATA[name] += len(ln)
	return DATA

def break_lines(seq,num):
	"""Adds newlines to a long sequence at set intervals"""
	import collections as cl

	seqLen = len(seq)
	seqPrint = cl.deque()

	start = 0
	end = int(num)
	while end < seqLen:
		seqPrint.append("{0}\n".format(seq[start:end]))
		start += num
		end += num

	seqPrint.append("{0}\n".format(seq[start:]))
	seqPrint = ''.join(seqPrint)

def oligo_builder(total):
	"""Builds list of all possible oligos of a set length"""
	oligo = []
	maxL = 4**total
	for a in range(maxL+1):
		oligo.append("")
	length = 1
	div = maxL / 4**length
	x = div
	rep = 1
	nuc =["A", "T", "G", "C"]
	nuc_num = 0

	while length < total:
		if rep <= x:
			oligo[rep] += nuc[nuc_num]
			rep += 1
		if rep == x + 1:
			nuc_num += 1
			if nuc_num == 4:
				nuc_num = 0
			x += div
		if rep == maxL + 1:
			rep = 1
			length += 1
			div = maxL / 4**length
			x = div

	nuc_num = 0
	for a in range(maxL+1):
		if a == 0:
			pass
		else:
			oligo[a] += (nuc[nuc_num])
			nuc_num += 1
			if nuc_num == 4:
				nuc_num = 0
	del oligo[0]
	tuple(oligo)
	return oligo

def list2string(List,slength):
	"Converts a list of characters into a list of strings: List, length of strings"
	length = len(List) / 4
	string_list = []
	oligo = ""
	for a in List:
		oligo += a
		if len(oligo) == slength:
			string_list.append(oligo)
			oligo = ""
	return string_list

def exact_match(data,findSeq,wLen,wStep):
	"""Identifies locations of a substring, and returns counts in sequence windows"""
	if not isinstance(data,str):
		raise TypeError('Sequence must be str')
	if not isinstance(find,str):
		raise TypeError('Sequence to search for must be a str')
	import collections as cl

	class Win(object):
		def __init__(self,start,end):
			self.start = start
			self.end = end
			self.count = 0

	pos = cl.deque()	#Start position of each match
	targetLength = len(findSeq)

	seg = data.split(findSeq)
	currentPos = 0
	for i in seg:
		pos.append(len(i) + 1 + currentPos)
		currentPos += len(i) + 1 + currentPos + targetLength
	pos.pop()
	pos = tuple(pos)
	totalMatch = len(pos)

	WIN = {}
	count = 0
	wStart = 0
	wEnd = 0
	while wStart <= len(data):
		wStart += wStep
		wEnd = wStart + wLen
		WIN[count] = Win(wStart,wEnd)
		count += 1
	for i in WIN:
		for x in pos:
			if x >= WIN[i].start and x <= WIN[i].end:
				WIN[i].count += 1
	return pos, WIN

def perc_id(seq1,seq2,ignore_ns=True):
	"""
	Calculate percent identity of an alignment.
	Faster if ignore_ns = False, but Ns and Xs will be counted as normal nucleotides
	"""

	if not ignore_ns:
		seq1 = list(seq1)
		seq2 = list(seq2)
		diff = float([x == y for x,y in zip(seq1,seq2)].count(True))
		total = float(len(seq1))
		perc = diff / total
		return (perc,diff,total)

	else:
		seq1 = tuple(seq1.upper())
		seq2 = tuple(seq2.upper())
		diff = 0.0
		comparisons = 0.0

		for i in range(0,len(seq1)):
			if seq1[i] != 'X' and seq1[i] != '-' and seq1[i] != 'N':
				if seq2[i] != 'X' and seq2[i] != '-' and seq1[i] != 'N':
					if seq1[i] != seq2[i]:
						diff += 1
					comparisons += 1

		if comparisons == 0:
			return (0,diff,comparisons)
		perc = 1 - (diff / comparisons)
		return (perc,diff,comparisons)

def codon_table():
	"""Builds a lookup table for codon -> aa"""
	CODON_TABLE = {
		"ACC":"T","ATG":"M","ACA":"T","ACG":"T","ATC":"I","AAC":"N","ATA":"I","AGG":"R",
		"CCT":"P","CTC":"L","AGC":"S","AAG":"K","AGA":"R","CAT":"H","AAT":"N","ATT":"I",
		"CTG":"L","CTA":"L","ACT":"T","CAC":"H","AAA":"K","CAA":"Q","AGT":"S","CCA":"P",
		"CCG":"P","CCC":"P","CTT":"L","TAT":"Y","GGT":"G","TGT":"C","CGA":"R","CAG":"Q",
		"CGC":"R","GAT":"D","CGG":"R","TTT":"F","TGC":"C","GGG":"G","TGA":"*","GGA":"G",
		"TGG":"W","GGC":"G","TAC":"Y","TTC":"F","TCG":"S","TAG":"*","TTG":"L","CGT":"R",
		"GAA":"E","TCA":"S","GCA":"A","GTA":"V","GCC":"A","GTC":"V","GCG":"A","GTG":"V",
		"GAG":"E","GTT":"V","GCT":"A","TTA":"L","GAC":"D","TCC":"S","TAA":"*","TCT":"S",
		'NNN':'x'}
	return CODON_TABLE