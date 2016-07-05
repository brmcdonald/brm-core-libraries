import json
import collections as cl

class Container(object):
	"""Empty object to be populated by JSON parser functions"""
	pass

def read_json(f):
	"""Reads JSON file and generates a dictionary of objects."""
	data = open(f,'r')
	try:
		DATA = json.load(data,object_pairs_hook=cl.OrderedDict)
	except ValueError as err:
		json_error(err,data,"Failed to parse JSON file {0}".format(f))
	data.close()
	RETURN = {}

	for obj in DATA:
		objName = str(obj)
		RETURN[objName] = Container()
		for attr in DATA[obj]:
			if type(DATA[obj][attr]) == unicode:
				RETURN[objName].__setattr__(attr,str(DATA[obj][attr]))
			elif type(DATA[obj][attr]) == cl.OrderedDict:
				RETURN[objName].__setattr__(attr,_decode_dict(DATA[obj][attr]))
			elif type(DATA[obj][attr]) == list:
				RETURN[objName].__setattr__(attr,_decode_list(DATA[obj][attr]))
			else:
				RETURN[objName].__setattr__(attr,DATA[obj][attr])
	return RETURN


def add_attr(DATA,attrList,initialVal=""):
	"""Adds attributes to a dictionary of objects"""
	if type(attrList) != list:
		attrList = [attrList] 
	for i in DATA:
		d = {key: key for key in dir(DATA[i])}
		for attr in attrList:
			if not attr.startswith('_'):
				if not d.has_key(attr): 
					DATA[i].__setattr__(attr,initialVal)
	return DATA


def obj2tsv(DATA,outFile,attrList=[],nameList=[]):
	"""Generates a TSV file from a dictionary of objects.
Can print a given set of attributes in a particular order, or prints all
by default. Attributes must be lists, strings, or numbers. All objects must have the same attributes"""
	out = open(outFile,'w')

	FOUND = {}
	if attrList:
		header = [str(x) for x in attrList]
		out.write("\t".join(header) + "\n")
	if nameList:
		for i in nameList:
			FOUND[i] = 0

	if FOUND:
		for g in nameList:
			FOUND[g] = 1
			data = []
			if attrList:
				for attr in attrList:
					if not attr.startswith('_'):
						if type(getattr(DATA[g],attr)) == list:
							data.append(','.join(getattr(DATA[g],attr)))
						else:
							data.append(str(getattr(DATA[g],attr)))
			else:
				for attr in dir(DATA[g]):
					if not attr.startswith('_'):
						if type(getattr(DATA[g],attr)) == list:
							data.append(','.join(getattr(DATA[g],attr)))
						else:
							data.append(str(getattr(DATA[g],attr)))
			out.write("\t".join(data) + '\n')
		out.close()

	else:
		for g in DATA:
			FOUND[g] = 1
			data = []
			if attrList:
				for attr in attrList:
					if not attr.startswith('_'):
						if type(getattr(DATA[g],attr)) == list:
							data.append(','.join(getattr(DATA[g],attr)))
						else:
							data.append(str(getattr(DATA[g],attr)))
			else:
				for attr in dir(DATA[g]):
					if not attr.startswith('_'):
						if type(getattr(DATA[g],attr)) == list:
							data.append(','.join(getattr(DATA[g],attr)))
						else:
							data.append(str(getattr(DATA[g],attr)))
			out.write("\t".join(data) + '\n')
		out.close()

	for i in FOUND:
		if FOUND[i] == 0:
			print "Could not find data for {0}".format(i)



def export_json(DATA,outFile):
	"""Export data from a dictionary of objects to JSON format."""
	out = open(outFile,'w')
	WRITE = cl.OrderedDict()

	inputType = type(DATA)
	if	inputType == cl.OrderedDict or inputType == dict:
		for f in DATA:
			WRITE[f] = cl.OrderedDict()
			for attr in dir(DATA[f]):
				if not attr.startswith("__"):
					WRITE[f][attr] = getattr(DATA[f],attr)

	json.dump(WRITE,out,indent=4)
	out.close()
	return 1



#Decoding JSON unicode to strings. Code found in a post by Stack Overflow user Mike Brennan
#https://stackoverflow.com/questions/956867/how-to-get-string-objects-instead-of-unicode-ones-from-json-in-python
def _decode_list(data):
	rv = []
	for item in data:
		if isinstance(item, unicode):
			item = item.encode('utf-8')
		elif isinstance(item, list):
			item = _decode_list(item)
		elif isinstance(item, dict):
			item = _decode_dict(item)
		rv.append(item)
	return rv

def _decode_dict(data):
	rv = cl.OrderedDict()
	for key, value in data.iteritems():
		if isinstance(key, unicode):
			key = key.encode('utf-8')
		if isinstance(value, unicode):
			value = value.encode('utf-8')
		elif isinstance(value, list):
			value = _decode_list(value)
		elif isinstance(value,cl.OrderedDict):
			value = _decode_dict(value)
		rv[key] = value
	return rv
