import xml.etree.ElementTree as ET
from py2neo import Graph, Node, Relationship, NodeMatcher, Subgraph
import csv


def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False


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
			self.compartments[i.attrib['id']] = [i.attrib['name'], i.attrib['spatialDimensions'], i.attrib['constant']]
		return self.compartments

	def getSpecies(self):
		self.species = {}
		for i in self.listSpecies:
			try:
				self.species[i.attrib['id']] = [i.attrib['name']]
			except:
				self.species[i.attrib['id']] = [i.attrib['id'][:-2]]
			try:
				self.species[i.attrib['id']].append(i.attrib['compartment'])
				self.species[i.attrib['id']].append(i.attrib['hasOnlySubstanceUnits'])
				self.species[i.attrib['id']].append(i.attrib['boundaryCondition'])
				self.species[i.attrib['id']].append(i.attrib['constant'])
				for note in i[0][0]:
					self.species[i.attrib['id']].append(note.text[note.text.index(':')+1:])
			except: pass
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
				for note in i[0][0]:
					try:
						self.reactions[i.attrib['id']].append(note.text[note.text.index(':')+1:])
					except: pass
				self.reactions[i.attrib['id']].append('|')
			except:
				self.reactions[i.attrib['id']].append('true')
				for note in i[0][0]:
					try:
						self.reactions[i.attrib['id']].append(note.text[note.text.index(':')+1:])
					except: pass
				self.reactions[i.attrib['id']].append('|')
			for child in i:
				if 'listOfReactants' in child.tag:
					for reac in child:
						self.reactions[i.attrib['id']].append(reac.attrib['species'])
						try:
							self.reactions[i.attrib['id']].append(float(reac.attrib['stoichiometry']))
						except:
							self.reactions[i.attrib['id']].append(1.0)
				elif 'listOfProducts' in child.tag:
					self.reactions[i.attrib['id']].append('|')
					for prod in child:
						self.reactions[i.attrib['id']].append(prod.attrib['species'])
						try:
							self.reactions[i.attrib['id']].append(float(prod.attrib['stoichiometry']))
						except:
							self.reactions[i.attrib['id']].append(1.0)
				elif 'kineticLaw' in child.tag:
					for parameter in child[1]:
						self.reactions[i.attrib['id']].insert(self.reactions[i.attrib['id']].index('|'), float(parameter.attrib['value']))
				else:
					print(child, 'not included in reaction', i.attrib['id'])
		return self.reactions

	def writeSpecNames(self):
		listDone = []
		with open('listIDs_RiceMultiomics.txt', 'w') as outfile:
			for name in self.species.keys():
				if self.species[name][0] not in listDone:
					outfile.write(self.species[name][0] + '\n')
					listDone.append(self.species[name][0])

	def readMetaIDs(self, file='res_RiceMultiomics.txt'):
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
	data = ProcessXML('RiceMultiomics.xml')

	graph = Graph('bolt://palsson.di.uminho.pt:6092', auth=('neo4j', '123'))
	compartments = data.getCompartments()
	species = data.getSpecies()
	reactions = data.getReactions()
	data.writeSpecNames()
	dicID = data.readMetaIDs()
	dicReact = data.getReactionDict()

	'''	for comp in compartments.keys():
		node = Node("Compartment", id=comp.upper(), name=compartments[comp][0], spatialDimensions=compartments[comp][1], constant=compartments[comp][2], model='RiceMulti')
		graph.create(node)'''

	matcher = NodeMatcher(graph)
	listUncSpec = []
	nodeModel = Node('Model', id='RiceMulti', name='RiceMulti')

	for spec in species.keys():
		try:
			id = dicID[species[spec][0]][0]
			if 'ModelSeed' in id:
				id = id.replace('=ModelSeed', '')
			else:
				id = id.replace('META:', '')
				id = id.replace('=MetaCyc', '')
		except:
			id = species[spec][0]
			listUncSpec.append([spec, species[spec][0]])
		try:
			nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + species[spec][1].upper() + '"')
			relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
			graph.create(relation)
		except:
			node = Node("Species", id=id, name=species[spec][0], compartment=species[spec][1].upper(), hasOnlySubstanceUnits=species[spec][2], boundaryCondition=species[spec][3], constant=species[spec][4])
			try:
				node['formula'] = species[spec][5]
				node['charge'] = float(species[spec][6])
			except:
				pass
			graph.create(node)
			relation = Relationship(nodeModel, 'CONTAINS', node)
			graph.create(relation)

	with open('UnconvertedSpeciesRiceMulti.txt', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(listUncSpec)

	# nodeModel = Node('Model', id='RiceMulti', name='RiceMulti')
	listUncReac = []

	for reac in reactions.keys():
		exists = False
		id = reac
		if 'tex' in id[-6:]:
			id = id.replace('tex', '')
			compart = 'C'
		elif 'tv' in id[-6:]:
			id = id.replace('tv', '')
			compart = 'V'
		elif 'tu' in id[-6:]:
			id = id.replace('tu', '')
			compart = 'U'
		elif 'tx' in id[-6:]:
			id = id.replace('tx', '')
			compart = 'X'
		elif 'ts' in id[-6:]:
			id = id.replace('ts', '')
			compart = 'P'
		elif 'tr' in id[-6:]:
			id = id.replace('tr', '')
			compart = 'R'
		elif 'tm' in id[-6:]:
			id = id.replace('tm', '')
			compart = 'M'
		elif 'Biomass' in id[-6:]:
			compart = 'C'
		elif 'c' in id[-6:]:
			id = id.replace('c', '')
			compart = 'C'
		elif 'm' in id[-6:]:
			id = id.replace('m', '')
			compart = 'M'
		elif 's' in id[-6:]:
			id = id.replace('s', '')
			compart = 'P'
		elif 'r' in id[-6:]:
			id = id.replace('r', '')
			compart = 'R'
		elif 'x' in id[-6:]:
			id = id.replace('x', '')
			compart = 'X'
		elif 'v' in id[-6:]:
			id = id.replace('v', '')
			compart = 'V'
		elif 'u' in id[-6:]:
			id = id.replace('u', '')
			compart = 'U'
		else:
			compart = 'C'
		try:
			id = dicReact[id]
		except:
			try:
				id = dicReact[id[2:]]
			except:
				listUncReac.append([reac, reactions[reac][0]])
		if matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + compart + "'").first() is None:
			node = Node('Reaction', id=id, compartment=compart)
		else:
			node = matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + compart + "'").first()
			exists = True
		reversible = True

		if not exists:
			node['name'] = reactions[reac][0]
			node['reversible'] = reactions[reac][1]
			node['geneAssociation'] = reactions[reac][2]
			node['subsystem'] = reactions[reac][3]
			node['ECnumber'] = reactions[reac][4]
			node['confidenceLevel'] = reactions[reac][5]
			#node['authors'] = reactions[reac][6]
			#node['lowerBound'] = reactions[reac][7]
			#node['upperBound'] = reactions[reac][8]
			node['fluxValue'] = reactions[reac][9]
			node['objectiveCoefficient'] = reactions[reac][10]

			try:
				if node['reversible'].lower() == 'false': reversible = False
			except:
				pass
			graph.create(node)

			for i in range(reactions[reac].index('|')+1, len(reactions[reac])):
				try:
					id = dicID[species[reactions[reac][i]][0]][0]
					compart = species[reactions[reac][i]][1].upper()
					if 'ModelSeed' in id:
						id = id.replace('=ModelSeed', '')
					else:
						id = id.replace('META:', '')
						id = id.replace('=MetaCyc', '')
				except:
					try:
						id = species[reactions[reac][i]][0]
						compart = species[reactions[reac][i]][1].upper()
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
				elif reactions[reac][i] == 'false' or reactions[reac][i] == 'true':
					pass
				elif matcher.match('Species').where("_.id='" + id + "'").where("_.compartment='" + compart + "'").first() is None:
					pass
				else:
					nodeSpec = matcher.match('Species').where("_.id='" + id + "'").where("_.compartment='" + compart + "'")
					relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=1.0)
					graph.create(relation)
					if reversible:
						relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=1.0)
						graph.create(relation)

			for i in range(j + 1, len(reactions[reac])):
				try:
					id = dicID[species[reactions[reac][i]][0]][0]
					compart = species[reactions[reac][i]][1].upper()
					if 'ModelSeed' in id:
						id = id.replace('=ModelSeed', '')
					else:
						id = id.replace('META:', '')
						id = id.replace('=MetaCyc', '')
				except:
					try:
						id = species[reactions[reac][i]][0]
						compart = species[reactions[reac][i]][1].upper()
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
				elif matcher.match('Species').where("_.id='" + id + "'").where("_.compartment='" + compart + "'").first() is None:
					pass
				else:
					nodeSpec = matcher.match('Species').where("_.id='" + id + "'").where("_.compartment='" + compart + "'")
					relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=1.0)
					graph.create(relation)
					if reversible:
						relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=1.0)
						graph.create(relation)
		relation = Relationship(nodeModel, 'CONTAINS', node, lowerBound=reactions[reac][7], upperBound=reactions[reac][8])
		graph.create(relation)

	with open('UnconvertedReactionsRiceMulti.txt', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(listUncReac)
