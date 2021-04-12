import xml.etree.ElementTree as ET
from py2neo import Graph, Node, Relationship, NodeMatcher, Subgraph
import csv


def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

def processASCII(string):
	import re
	findstr = re.findall('_[0-9]{2}_', string)
	for i in findstr:
		tempstr = int(i.strip('_'))
		string = string.replace(i, chr(tempstr))
	return string

class ProcessXML:

	def __init__(self, model):
		root = ET.parse('Models/' + model).getroot()
		model = root[0]
		for i in model:
			if 'comp' in i.tag.lower():
				self.listCompartments = i
			elif 'spec' in i.tag.lower():
				self.listSpecies = i
			elif 'reac' in i.tag.lower():
				self.listReactions = i
			else:
				print(i,'not included in model')

	def getCompartments(self):
		self.compartments = {}
		for i in self.listCompartments:
			self.compartments[i.attrib['id'][2:].title()] = i.attrib['name']
		return self.compartments

	def getSpecies(self):
		self.species = {}
		dicID = self.readMetaIDs()
		self.listUncSpec = []
		for i in self.listSpecies:
			id = processASCII(i.attrib['id'])[2:]
			pos = id.rfind('_')
			compart = i.attrib['compartment'][2:].title()
			id = id[:pos]
			try:
				id = dicID[id][0]
				if 'ModelSeed' in id:
					id = id.replace('=ModelSeed', '')
				else:
					id = id.replace('META:', '')
					id = id.replace('=MetaCyc', '')
			except:
				self.listUncSpec.append([processASCII(i.attrib['id']), i.attrib['name']])
			id = id + '_' + compart
			if id in self.species.keys():
				print(id, '\t', i.attrib['id'], '\t', self.species[id])
			self.species[id] = []
			try:
				pos = 1000
				if '_' in i.attrib['name']:
					pos = i.attrib['name'].rfind('_')
				self.species[id].append(i.attrib['name'][:pos])
			except:
				pos = i.attrib['id'].rfind('_')
				self.species[id].append(processASCII(i.attrib['id'])[2:pos])
			try:
				self.species[id].append(i.attrib['compartment'][2:].title())
			except:
				for comp in self.compartments.keys():
					if processASCII(i.attrib['id'])[-2:] in comp:
						self.species[id].append(comp)
		return self.species

	def getReactions(self):
		self.reactions = {}
		for i in self.listReactions:
			try:
				self.reactions[i.attrib['id']] = [i.attrib['name']]
			except:
				self.reactions[i.attrib['id']] = []
			try:
				self.reactions[i.attrib['id']].append(i.attrib['reversible'])
				self.reactions[i.attrib['id']].append('|')
			except:
				self.reactions[i.attrib['id']].append('true')
				self.reactions[i.attrib['id']].append('|')
			for child in i:
				if 'listOfReactants' in child.tag:
					for reac in child:
						id = processASCII(reac.attrib['species'])[2:]
						pos = id.rfind('_')
						compart = id[pos + 1:].title()
						id = id[:pos]
						try:
							id = dicID[id][0]
							if 'ModelSeed' in id:
								id = id.replace('=ModelSeed', '')
							else:
								id = id.replace('META:', '')
								id = id.replace('=MetaCyc', '')
						except:
							pass
						self.reactions[i.attrib['id']].append(id + '_' + compart)
						try:
							self.reactions[i.attrib['id']].append(float(reac.attrib['stoichiometry']))
						except:
							self.reactions[i.attrib['id']].append(1.0)
				elif 'listOfProducts' in child.tag:
					self.reactions[i.attrib['id']].append('|')
					for prod in child:
						id = processASCII(prod.attrib['species'])[2:]
						pos = id.rfind('_')
						compart = id[pos + 1:].title()
						id = id[:pos]
						try:
							id = dicID[id][0]
							if 'ModelSeed' in id:
								id = id.replace('=ModelSeed', '')
							else:
								id = id.replace('META:', '')
								id = id.replace('=MetaCyc', '')
						except:
							pass
						self.reactions[i.attrib['id']].append(id + '_' + compart)
						try:
							self.reactions[i.attrib['id']].append(float(prod.attrib['stoichiometry']))
						except:
							self.reactions[i.attrib['id']].append(1.0)
				else:
					pass
					# print(child, 'not included in reaction')
		return self.reactions

	def writeSpecNames(self):
		listDone = []
		with open('listNames_C4GEM.txt', 'w') as outfile:
			for name in self.species.keys():
				if '_ext' in self.species[name][0] or '_acc' in self.species[name][0] or '_pmd' in self.species[name][0]:
					if species[name][0][:-4] not in listDone:
						outfile.write(species[name][0][:-4] + '\n')
						listDone.append(species[name][0][:-4])
				elif '_c' in self.species[name][0] or '_m' in self.species[name][0] or '_p' in self.species[name][0] or '_v' in \
						self.species[name][0] or '_x' in self.species[name][0]:
					if self.species[name][0][:-2] not in listDone:
						outfile.write(self.species[name][0][:-2] + '\n')
						listDone.append(self.species[name][0][:-2])
				elif '_biomass' in self.species[name][0]:
					if self.species[name][0][:-8] not in listDone:
						outfile.write(self.species[name][0][:-8] + '\n')
						listDone.append(self.species[name][0][:-8])
				elif '_leaf' in self.species[name][0] or '_accL' in self.species[name][0]:
					if self.species[name][0][:-5] not in listDone:
						outfile.write(self.species[name][0][:-5] + '\n')
						listDone.append(self.species[name][0][:-5])
				elif '_total' in self.species[name][0]:
					if self.species[name][0][:-6] not in listDone:
						outfile.write(self.species[name][0][:-6] + '\n')
						listDone.append(self.species[name][0][:-6])
				else:
					if self.species[name][0] not in listDone:
						outfile.write(self.species[name][0] + '\n')
						listDone.append(self.species[name][0])

	def readMetaIDs(self, file='res_C4GEM.txt'):
		import pandas as pd
		ids_data = pd.read_csv(file, header=None, sep='\t')
		self.dicID = {}
		for line in ids_data.index:
			temp = ids_data.iloc[line][1]
			try:
				temp = temp.strip('{}')
				self.dicID[ids_data.iloc[line][0]] = temp.split(',')
			except:
				self.dicID[ids_data.iloc[line][0]] = None
		return self.dicID

	def getReactionDict(self, file='Unique_ModelSEED_Reaction_Aliases.txt'):
		import csv
		dicReact = {}
		with open(file, 'r') as csvfile:
			spamreader = csv.reader(csvfile, delimiter='\t')
			for row in spamreader:
				if row[1] not in dicReact.keys():
					dicReact[row[1]] = row[0]
		return dicReact


if __name__ == '__main__':
	organism = "Maize" # "Sorghum" "SugarCane"
	data = ProcessXML('166488' + organism + '_C4GEM_vs1.0.txt')

	graph = Graph('bolt://palsson.di.uminho.pt:6092', auth=('neo4j', '123'))
	compartments = data.getCompartments()
	species = data.getSpecies()
	reactions = data.getReactions()
	dicID = data.readMetaIDs()
	dicReact = data.getReactionDict()
	listUncSpec = data.listUncSpec

	with open('UnconvertedSpecies' + organism + 'C4GEM.txt', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(listUncSpec)


	'''	for comp in compartments.keys():
		node = Node("Compartment", id=comp, name=compartments[comp], model=organism + 'C4GEM')
		graph.create(node)'''

	matcher = NodeMatcher(graph)
	nodeModel = Node('Model', id=organism + 'C4GEM', name=organism + 'C4GEM')

	for spec in species.keys():
		try:
			pos = spec.rfind('_')
			nodeSpec = matcher.match('Species').where('_.id="' + spec[:pos] + '"').where('_.compartment="' + species[spec][1] + '"')
			relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
			graph.create(relation)
		except:
			pos = spec.rfind('_')
			node = Node("Species", id=spec[:pos], name=species[spec][0], compartment=species[spec][1])
			graph.create(node)
			relation = Relationship(nodeModel, 'CONTAINS', node)
			graph.create(relation)

	# nodeModel = Node('Model', id=organism + 'C4GEM', name=organism + 'C4GEM')
	listUncReac = []

	for reac in reactions.keys():
		compart = 'C'
		if 'R_TC' in reac:
			compart = reac[4]
		exists = False
		id = reac[2:]
		if '_' in id:
			pos = id.rfind('_')
			compart = id[pos+1:].title()
			id = id[:pos]
		try:
			id = dicReact[id]
		except:
			try:
				id = dicReact['R_' + id]
			except:
				listUncReac.append([reac, reactions[reac][0]])
		if matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + compart + "'").first() is None:
			node = Node('Reaction', id=id, compartment=compart, name=reac)
		else:
			node = matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + compart + "'").first()
			exists = True
		reversible = True

		if not exists:
			for j in range(len(reactions[reac])):
				if reactions[reac][j] == '|':
					break
				if reactions[reac][j] == 'true' or reactions[reac][j] == 'false':
					node['reversible'] = reactions[reac][j]
				elif reactions[reac][j] != '|' and reactions[reac][j] != 'true' and reactions[reac][j] != 'false':
					node['name'] = reactions[reac][j]
			try:
				if node['reversible'].lower() == 'false': reversible = False
			except:
				pass
			graph.create(node)

			for i in range(j + 1, len(reactions[reac])):
				try:
					pos = reactions[reac][i].rfind('_')
					id = reactions[reac][i][:pos]
					compart = reactions[reac][i][pos+1:]
				except:
					id = ''
					compart = ''
				if reactions[reac][i] == '|':
					j = i
					break
				if isfloat(reactions[reac][i]):
					relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=float(reactions[reac][i]))
					graph.create(relation)
					if reversible:
						relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=float(reactions[reac][i]))
						graph.create(relation)
				elif matcher.match('Species').where('_.id="' + id + '"').where("_.compartment='" + compart + "'").first() is None:
					pass
				else:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where("_.compartment='" + compart + "'")
					relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=1.0)
					graph.create(relation)
					if reversible:
						relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=1.0)
						graph.create(relation)

			for i in range(j + 1, len(reactions[reac])):
				try:
					pos = reactions[reac][i].rfind('_')
					id = reactions[reac][i][:pos]
					compart = reactions[reac][i][pos+1:]
				except:
					id = ''
					compart = ''
				if isfloat(reactions[reac][i]):
					relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=float(reactions[reac][i]))
					graph.create(relation)
					if reversible:
						relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=float(reactions[reac][i]))
						graph.create(relation)
				elif reactions[reac][i] == 'false' or reactions[reac][i] == 'true':
					pass
				elif matcher.match('Species').where('_.id="' + id + '"').where("_.compartment='" + compart + "'").first() is None:
					pass
				else:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where("_.compartment='" + compart + "'")
					relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=1.0)
					graph.create(relation)
					if reversible:
						relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=1.0)
						graph.create(relation)
		relation = Relationship(nodeModel, 'CONTAINS', node)
		graph.create(relation)

	with open('UnconvertedReactions' + organism + 'C4GEM.txt', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(listUncReac)
