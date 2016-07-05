from operator import itemgetter
import collections as cl
#from __future__ import print_function

#Some common color definitions
COLOR = {'blue':'#335B96','red':'#E41B17','green':'#347235','purple':'#A900B5','gold':'#FDD017',
	'black':'#2B1B17','teal':'#008080','maroon':'#800000','brown':'#966F33'}

def is_number(s):
	"""Figure out if something is convertable to a number"""
	try:
		float(s)		
		return True
	except ValueError:	
		return False


def range_overlap(aStart,aEnd,bStart,bEnd,max_overlap=0):
	"""Determines if two ranges of numbers overlap"""
	iCoords = set(range(aStart, aEnd))
	jCoords = set(range(bStart, bEnd))

	if len(iCoords.intersection(jCoords)) > max_overlap:
		return True
	else:
		return False

def clade2dict(f,invert=False):
	"""Parse fasta-type clade divisions"""
	DATA = {}
	if invert == False:
		with open(f,'r') as info:
			for ln in info:
				ln = ln.replace('\n','')
				if '>' in ln:
					clade = ln.replace('>','')
					DATA[clade] = {}
				elif ln:
					ln = ln.split('\t')
					if len(ln) > 1:
						DATA[clade][ln[0]] = ln[1]
					else:
						DATA[clade][ln[0]] = ''
	elif invert == True:
		with open(f,'r') as info:
			for ln in info:
				ln = ln.replace('\n','')
				if '>' in ln:
					clade = ln.replace('>','')
				elif ln:
					ln = ln.split('\t')
					g = ln[0]
					if not DATA.has_key(g):
						DATA[g] = {}
						DATA[g][clade] = ''
					else:
						DATA[g][clade] = ''
	return DATA

#string search with overlaps
def match_overlap(string, sub):
	count = start = 0
	while True:
		start = string.find(sub, start) + 1
		if start > 0:
			count+=1
		else:
			return count

def r_fet(a,b,c,d,name):
	"""Returns text required to run FET in R"""

	fet_input = ("""
	over <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
	pval <- fisher.test(over, alternative=\"greater\")$p.value
	odds <- fisher.test(over, alternative=\"greater\")$estimate
	fet_data_over <- rbind(fet_data_over, c("{4}",pval,odds))

	under <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
	pval <- fisher.test(under, alternative=\"less\")$p.value
	odds <- fisher.test(under, alternative=\"less\")$estimate
	fet_data_under <- rbind(fet_data_under, c("{4}",pval,odds))

	both <- matrix(c({0}, {1}, {2}, {3}), nrow=2)
	pval <- fisher.test(both, alternative=\"two.sided\")$p.value
	odds <- fisher.test(both, alternative=\"two.sided\")$estimate
	fet_data_both <- rbind(fet_data_both, c("{4}",pval,odds))\n""".format(a,b,c,d,name))

	return fet_input

def jaccard_index(a,b):
	set_1 = set(a)
	set_2 = set(b)
	n = len(set_1.intersection(set_2))
	return round(n / float(len(set_1) + len(set_2) - n) ,5)

def goodman_kruskal_gamma(m, n):
	from itertools import permutations

	mSet = set(m)
	if len(mSet) < len(m):
		print 'Every domain must only appear once for goodman_kruskal_gamma'
		return 0
	nSet = set(n)
	if len(nSet) < len(n):
		print 'Every domain must only appear once for goodman_kruskal_gamma'
		return 0

	MPAIRS = {}
	for idx, i in enumerate(m):
		for idx2, j in enumerate(m):
			if idx2 > idx:
				MPAIRS['{0}{1}'.format(i,j)] = ''

	NPAIRS = {}
	for idx, i in enumerate(n):
		for idx2, j in enumerate(n):
			if idx2 > idx:
				NPAIRS['{0}{1}'.format(i,j)] = ''
	print MPAIRS
	print NPAIRS

	all_pairs = set(MPAIRS.keys() + NPAIRS.keys())
	total_pairs = len(set(MPAIRS.keys() + NPAIRS.keys()))
	pos = 0.0
	rev = 0.0

	for i in all_pairs:
		#Both have the pair, in the same order
		if NPAIRS.has_key(i) and MPAIRS.has_key(i):
			pos += 1
		#Reversed in second protein, forward in first
		elif NPAIRS.has_key(i[::-1]) and MPAIRS.has_key(i):
			rev += 1
		#Reversed in first protein, forward in second
		elif MPAIRS.has_key(i[::-1]) and NPAIRS.has_key(i):
			rev += 1

	statistic = (1 + ((pos - rev) / total_pairs)) / 2

	if pos + rev <= 1:
		return (0,pos,rev,total_pairs)
	else:
		return (round(statistic,5),pos,rev,total_pairs)

def tbl2dict(f,remove='',invert=False,valueCol=1):
	"""Converts a tab-delimited table into a dictionary: dict[col0] = valueCol"""
	import collections as cl
	DATA = cl.OrderedDict()
	with open(f,'r') as info:
		for ln in info:
			ln = ln.replace('\n','')
			if remove:
				ln = ln.replace(remove,'')
			ln = ln.split('\t')
			if len(ln) > 1 and invert == False:
				DATA[ln[0]] = ln[valueCol]
			if len(ln) > 1 and invert == True:
				DATA[ln[valueCol]] = ln[0]
			elif len(ln) == 1:
				DATA[ln[0]] = ''
	return DATA

def tbl2dictdict(f,remove=''):
	"""Converts a tab-delimited table into a dictionary: dict[col1] = {other cols}"""
	DATA = cl.OrderedDict()
	with open(f,'r') as info:
		for ln in info:
			ln = ln.replace('\n','')
			if remove:
				ln = ln.replace(remove,'')
			ln = ln.split('\t')
			if ln[0]:
				name = ln.pop(0)
				DATA[name] = cl.OrderedDict()
				if len(ln) > 0:
					for i in ln:
						DATA[name][i] = ''
	return DATA

def tbl2dictlist(f,remove=''):
	"""Converts a tab-delimited table into a dictionary: dict[col1] = [other cols]"""
	DATA = cl.OrderedDict()
	with open(f,'r') as info:
		for ln in info:
			ln = ln.replace('\n','')
			if remove:
				ln = ln.replace(remove,'')
			ln = ln.split('\t')
			if ln[0]:
				name = ln.pop(0)
				DATA[name] = []
				if len(ln) > 0:
					for i in ln:
						DATA[name].append(i)
	return DATA

def dict2matrix(DATA,f,rowOrder=[],colOrder=[]):
	"""Prints a 2 layer dict 'DATA' as a matrix to an output file 'f'"""
	COL = {}
	for i in DATA:
		for j in DATA[i]:
			COL[j] = ''
	for i in DATA:
		for j in COL:
			if not DATA[i].has_key(j):
				DATA[i][j] = 0

	out = open(f,'w')
	out.write('GENOME')
	if rowOrder and not colOrder:
		for i in sorted(COL.items(), key=lambda x: x[1],reverse=True):
			out.write('\t{0}'.format(i[0]))
		out.write('\n')
		for g in rowOrder:
			out.write(g)
			for i in sorted(COL.items(), key=lambda x: x[1], reverse=True):
				out.write('\t{0}'.format(DATA[g][i[0]]))
			out.write('\n')
	elif not rowOrder and colOrder:
		for i in colOrder:
			out.write('\t{0}'.format(i[0]))
		out.write('\n')
		for g in DATA:
			out.write(g)
			for i in colOrder:
				out.write('\t{0}'.format(DATA[g][i[0]]))
			out.write('\n')
	elif rowOrder and colOrder:
		for i in colOrder:
			out.write('\t{0}'.format(i[0]))
		out.write('\n')
		for g in rowOrder:
			out.write(g)
			for i in colOrder:
				out.write('\t{0}'.format(DATA[g][i[0]]))
			out.write('\n')
	else:
		for i in sorted(COL.items(), key=lambda x: x[1],reverse=True):
			out.write('\t{0}'.format(i[0]))
		out.write('\n')
		for g in DATA:
			out.write(g)
			for i in sorted(COL.items(), key=lambda x: x[1], reverse=True):
				out.write('\t{0}'.format(DATA[g][i[0]]))
			out.write('\n')
	out.close()

def matrix2dict(f,sep='\t'):
	"""Converts a 2 dimensional matrix to a dictionary"""

	MAT1 = {}
	data1 = open(f,'r')
	f = data1.readlines()
	data1.close()
	header = f.pop(0).replace('\n','').split(sep)
	cols1 = header[1:]

	for ln in f:
		ln = ln.replace('\n','')
		ln = ln.split(sep)
		row_name = ln.pop(0)
		MAT1[row_name] = {}
		for idx, i in enumerate(ln):
			MAT1[row_name][cols1[idx]] = i
	return MAT1

def ortho2dict(f,invert=False,setup='list'):
	"""Converts an orthoMCL-style table into a dict:\ndict[fam] = [prots]
or
dict[fam] = {prots} (setup='dict'). Can also invert dictionary using invert=True:
dict[prot] = fam or group by genome dict[gemome] = {fams}"""

	CAT = {}
	with open(f,'r') as tbl_data:
		if invert == False:
			for ln in tbl_data:
				ln = ln.replace('\n','')
				data = ln.split('\t')
				cat = data[0].replace(':','')
				if setup == 'list':
					CAT[cat] = []
				elif setup == 'dict':
					CAT[cat] = {}

				genes = data[1].split(' ')
				for i in genes:
					if '|' not in i:
						genome = i.split('_')[0]
						i = '{0}|{1}'.format(genome,i)
					if setup == 'list':
						CAT[cat].append(i)
					elif setup == 'dict':
						CAT[cat][i] = ''

		elif invert == True:
			for ln in tbl_data:
				ln = ln.replace('\n','')
				data = ln.split('\t')
				cat = data[0].replace(':','')

				genes = data[1].split(' ')
				for i in genes:
					if '|' not in i:
						genome = i.split('_')[0]
						i = '{0}|{1}'.format(genome,i)
					CAT[i] = cat
		elif setup == 'genome':
			for ln in tbl_data:
				ln = ln.replace('\n','')
				data = ln.split('\t')
				cat = data[0].replace(':','')

				genes = data[1].split(' ')
				for i in genes:
					genome = i.split('|')[0]
					if not CAT.has_key(genome):
						CAT[genome] = {}
					CAT[genome][cat] = ''
	return CAT


def most_common_item(L):
	"""Get the most common item in a list"""
	import itertools
	import operator

	# get an iterable of (item, iterable) pairs
	SL = sorted((x, i) for i, x in enumerate(L))
	# print 'SL:', SL
	groups = itertools.groupby(SL, key=operator.itemgetter(0))
	# auxiliary function to get "quality" for an item
	def _auxfun(g):
		item, iterable = g
		count = 0
		min_index = len(L)
		for _, where in iterable:
			count += 1
			min_index = min(min_index, where)
		# print 'item %r, count %r, minind %r' % (item, count, min_index)
		return count, -min_index
	# pick the highest-count/earliest item
	return max(groups, key=_auxfun)[0]

def sortDict(d,reverse=False):
	"""Sort a dictionary by values"""
	#http://writeonly.wordpress.com/2008/08/30/sorting-dictionaries-by-value-in-python-improved/
	return sorted(d.iteritems(), key=itemgetter(1), reverse=True)


need_to_fix = '''
def total_size(o, handlers={}, verbose=False):
	from sys import getsizeof, stderr
	from itertools import chain

	try:
		from reprlib import repr
	except ImportError:
		pass

	""" Returns the approximate memory footprint an object and all of its contents.

	Automatically finds the contents of the following builtin containers and
	their subclasses:  tuple, list, deque, dict, set and frozenset.
	To search other containers, add handlers to iterate over their contents:

		handlers = {SomeContainerClass: iter,
					OtherContainerClass: OtherContainerClass.get_elements}

	"""
	dict_handler = lambda d: chain.from_iterable(d.items())
	all_handlers = {tuple: iter,
					list: iter,
					cl.deque: iter,
					dict: dict_handler,
					set: iter,
					frozenset: iter,
				   }
	all_handlers.update(handlers)	 # user handlers take precedence
	seen = set()					  # track which object id's have already been seen
	default_size = getsizeof(0)	   # estimate sizeof object without __sizeof__

	def sizeof(o):
		if id(o) in seen:	   # do not double count the same object
			return 0
		seen.add(id(o))
		s = getsizeof(o, default_size)

		if verbose:
			print(s, type(o), repr(o), file=stderr)

		for typ, handler in all_handlers.items():
			if isinstance(o, typ):
				s += sum(map(sizeof, handler(o)))
				break
		return s

	return sizeof(o)
'''


class TinyDict(cl.MutableMapping):
	"""Space efficient dictionary with fast iteration and cheap resizes."""

	import array
	import itertools

	# Placeholder constants
	FREE = -1
	DUMMY = -2

	@staticmethod
	def _gen_probes(hashvalue, mask):
		'Same sequence of probes used in the current dictionary design'
		PERTURB_SHIFT = 5
		if hashvalue < 0:
			hashvalue = -hashvalue
		i = hashvalue & mask
		yield i
		perturb = hashvalue
		while True:
			i = (5 * i + perturb + 1) & 0xFFFFFFFFFFFFFFFF
			yield i & mask
			perturb >>= PERTURB_SHIFT

	def _lookup(self, key, hashvalue):
		'Same lookup logic as currently used in real dicts'
		assert self.filled < len(self.indices)   # At least one open slot
		freeslot = None
		for i in self._gen_probes(hashvalue, len(self.indices)-1):
			index = self.indices[i]
			if index == FREE:
				return (FREE, i) if freeslot is None else (DUMMY, freeslot)
			elif index == DUMMY:
				if freeslot is None:
					freeslot = i
			elif (self.keylist[index] is key or
				  self.hashlist[index] == hashvalue
				  and self.keylist[index] == key):
					return (index, i)

	@staticmethod
	def _make_index(n):
		'New sequence of indices using the smallest possible datatype'
		if n <= 2**7: return array.array('b', [FREE]) * n	   # signed char
		if n <= 2**15: return array.array('h', [FREE]) * n	  # signed short
		if n <= 2**31: return array.array('l', [FREE]) * n	  # signed long
		return [FREE] * n									   # python integers

	def _resize(self, n):
		'''Reindex the existing hash/key/value entries.
		   Entries do not get moved, they only get new indices.
		   No calls are made to hash() or __eq__().

		'''
		n = 2 ** n.bit_length()					 # round-up to power-of-two
		self.indices = self._make_index(n)
		for index, hashvalue in enumerate(self.hashlist):
			for i in Dict._gen_probes(hashvalue, n-1):
				if self.indices[i] == FREE:
					break
			self.indices[i] = index
		self.filled = self.used

	def clear(self):
		self.indices = self._make_index(8)
		self.hashlist = []
		self.keylist = []
		self.valuelist = []
		self.used = 0
		self.filled = 0										 # used + dummies

	def __getitem__(self, key):
		hashvalue = hash(key)
		index, i = self._lookup(key, hashvalue)
		if index < 0:
			raise KeyError(key)
		return self.valuelist[index]

	def __setitem__(self, key, value):
		hashvalue = hash(key)
		index, i = self._lookup(key, hashvalue)
		if index < 0:
			self.indices[i] = self.used
			self.hashlist.append(hashvalue)
			self.keylist.append(key)
			self.valuelist.append(value)
			self.used += 1
			if index == FREE:
				self.filled += 1
				if self.filled * 3 > len(self.indices) * 2:
					self._resize(4 * len(self))
		else:
			self.valuelist[index] = value

	def __delitem__(self, key):
		hashvalue = hash(key)
		index, i = self._lookup(key, hashvalue)
		if index < 0:
			raise KeyError(key)
		self.indices[i] = DUMMY
		self.used -= 1
		# If needed, swap with the lastmost entry to avoid leaving a "hole"
		if index != self.used:
			lasthash = self.hashlist[-1]
			lastkey = self.keylist[-1]
			lastvalue = self.valuelist[-1]
			lastindex, j = self._lookup(lastkey, lasthash)
			assert lastindex >= 0 and i != j
			self.indices[j] = index
			self.hashlist[index] = lasthash
			self.keylist[index] = lastkey
			self.valuelist[index] = lastvalue
		# Remove the lastmost entry
		self.hashlist.pop()
		self.keylist.pop()
		self.valuelist.pop()

	def __init__(self, *args, **kwds):
		if not hasattr(self, 'keylist'):
			self.clear()
		self.update(*args, **kwds)

	def __len__(self):
		return self.used

	def __iter__(self):
		return iter(self.keylist)

	def iterkeys(self):
		return iter(self.keylist)

	def keys(self):
		return list(self.keylist)

	def itervalues(self):
		return iter(self.valuelist)

	def values(self):
		return list(self.valuelist)

	def iteritems(self):
		return itertools.izip(self.keylist, self.valuelist)

	def items(self):
		return zip(self.keylist, self.valuelist)

	def __contains__(self, key):
		index, i = self._lookup(key, hash(key))
		return index >= 0

	def get(self, key, default=None):
		index, i = self._lookup(key, hash(key))
		return self.valuelist[index] if index >= 0 else default

	def popitem(self):
		if not self.keylist:
			raise KeyError('popitem(): dictionary is empty')
		key = self.keylist[-1]
		value = self.valuelist[-1]
		del self[key]
		return key, value

	def __repr__(self):
		return 'Dict(%r)' % self.items()

def numdecimal(n,decimals):
	"""Returns a float with the appropriate number of decimal places"""
	return float("{0:.{1}f}".format(n,decimals))

#Splits a single lists into multiple lists of a set length
def grouper(n, iterable, padvalue='N'):
	"grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
	return itertools.izip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

def distance_on_unit_sphere(lat1, long1, lat2, long2,unit='km'):
	#Code from http://www.johndcook.com/blog/python_longitude_latitude/

	import math
	# Convert latitude and longitude to 
	# spherical coordinates in radians.
	degrees_to_radians = math.pi/180.0
		 
	# phi = 90 - latitude
	phi1 = (90.0 - lat1)*degrees_to_radians
	phi2 = (90.0 - lat2)*degrees_to_radians
		 
	# theta = longitude
	theta1 = long1*degrees_to_radians
	theta2 = long2*degrees_to_radians
		 
	# Compute spherical distance from spherical coordinates.
		 
	# For two locations in spherical coordinates 
	# (1, theta, phi) and (1, theta, phi)
	# cosine( arc length ) = 
	#	sin phi sin phi' cos(theta-theta') + cos phi cos phi'
	# distance = rho * arc length

	cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
		math.cos(phi1)*math.cos(phi2))
	arc = math.acos( cos )

	if unit == 'km':
		arc = arc * 6373
	elif unit == 'miles':
		arc = arc * 3960
	return arc
