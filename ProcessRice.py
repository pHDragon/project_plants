import xml.etree.ElementTree as ET
from py2neo import Graph, Node, Relationship, NodeMatcher, Subgraph


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
				print(i)

	def getCompartments(self):
		self.compartments = {}
		for i in self.listCompartments:
			self.compartments[i.attrib['id']] = [i.attrib['name'], i.attrib['size'], i.attrib['constant'], 'Rice']
		return self.compartments

	def getSpecies(self):
		self.species = {}
		for i in self.listSpecies:
			try:
				self.species[i.attrib['id']] = [i.attrib['name']]
			except:
				self.species[i.attrib['id']] = [i.attrib['id']]
			try:
				self.species[i.attrib['id']].append(i.attrib['compartment'])
				self.species[i.attrib['id']].append(float(i.attrib['initialConcentration']))
				self.species[i.attrib['id']].append(i.attrib['hasOnlySubstanceUnits'])
				self.species[i.attrib['id']].append(i.attrib['boundaryCondition'])
				self.species[i.attrib['id']].append(i.attrib['constant'])
			except:
				pass
		return self.species

	def getReactions(self):
		self.reactions = {}
		for i in self.listReactions:
			if 'mit_' in i.attrib['id']:
				id = i.attrib['id']
				compart = 'M'
			elif 'chl_' in i.attrib['id']:
				id = i.attrib['id']
				compart = 'H'
			else:
				id = i.attrib['id']
				compart = 'C'
			if id in self.reactions.keys():
				print(id, i.attrib['id'], self.reactions[id])
			try:
				self.reactions[id] = [i.attrib['name'], compart, i.attrib['reversible'], i.attrib['fast'], '|']
			except:
				self.reactions[id] = [i.attrib['name'], compart, 'true', i.attrib['fast'], '|']
			for child in i:
				if 'listOfReactants' in child.tag:
					for reac in child:
						self.reactions[id].append(reac.attrib['species'])
						try:
							self.reactions[id].append(float(reac.attrib['stoichiometry']))
							self.reactions[id].append(reac.attrib['constant'])
						except:
							self.reactions[id].append(1.0)
							self.reactions[id].append(reac.attrib['constant'])
				elif 'listOfProducts' in child.tag:
					self.reactions[id].append('|')
					for prod in child:
						self.reactions[id].append(prod.attrib['species'])
						try:
							self.reactions[id].append(float(prod.attrib['stoichiometry']))
							self.reactions[id].append(reac.attrib['constant'])
						except:
							self.reactions[id].append(1.0)
							self.reactions[id].append(reac.attrib['constant'])
		return self.reactions

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
	data = ProcessXML('Rice.sbml')
	graph = Graph('bolt://localhost:7687', auth=('neo4j', '123'))
	compartments = data.getCompartments()
	species = data.getSpecies()
	reactions = data.getReactions()
	dicReact = data.getReactionDict()

	'''	nodeC = Node("Compartment", id='C', name='Cytosol', size=compartments['DefaultCompartment'][1], constant=compartments['DefaultCompartment'][2], model=compartments['DefaultCompartment'][3])
	graph.create(nodeC)
	nodeM = Node("Compartment", id='M', name='Mitochondrion', size=compartments['DefaultCompartment'][1], constant=compartments['DefaultCompartment'][2], model=compartments['DefaultCompartment'][3])
	graph.create(nodeM)
	nodeX = Node("Compartment", id='E', name='Extracellular', size=compartments['DefaultCompartment'][1], constant=compartments['DefaultCompartment'][2], model=compartments['DefaultCompartment'][3])
	graph.create(nodeX)
	nodeH = Node("Compartment", id='H', name='Chloroplast', size=compartments['DefaultCompartment'][1], constant=compartments['DefaultCompartment'][2], model=compartments['DefaultCompartment'][3])
	graph.create(nodeH)'''

	nodeModel = Node('Model', id='Rice', name='Rice')
	matcher = NodeMatcher(graph)

	for spec in species.keys():
		try:
			if '_mit' in species[spec][0] or 'mit_' in species[spec][0]:
				try:
					id = species[spec][0].replace('_mit', '')
				except:
					id = species[spec][0].replace('mit_', '')
				nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'M' + '"')
				relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
				graph.create(relation)
			elif 'x_' in species[spec][0] or '_x' in species[spec][0]:
				id = species[spec][0].replace('x_', '')
				nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'E' + '"')
				relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
				graph.create(relation)
			elif 'str_' in species[spec][0] or '_str' in species[spec][0]:
				id = species[spec][0].replace('_str', '')
				nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'H' + '"')
				relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
				graph.create(relation)
			else:
				nodeSpec = matcher.match('Species').where('_.id="' + species[spec][0] + '"').where('_.compartment="' + 'C' + '"')
				relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
				graph.create(relation)
		except:
			if '_mit' in species[spec][0] or 'mit_' in species[spec][0]:
				try:
					id = species[spec][0].replace('_mit', '')
				except:
					id = species[spec][0].replace('mit_', '')
				node = Node("Species", id=id, name=species[spec][0], compartment='M', hasOnlySubstanceUnits=species[spec][3], boundaryCondition=species[spec][4], constant=species[spec][5])
				graph.create(node)
				relation = Relationship(nodeModel, 'CONTAINS', node)
				graph.create(relation)
			elif 'x_' in species[spec][0] or '_x' in species[spec][0]:
				id = species[spec][0].replace('x_', '')
				node = Node("Species", id=id, name=species[spec][0], compartment='E', hasOnlySubstanceUnits=species[spec][3], boundaryCondition=species[spec][4], constant=species[spec][5])
				graph.create(node)
				relation = Relationship(nodeModel, 'CONTAINS', node)
				graph.create(relation)
			elif 'str_' in species[spec][0] or '_str' in species[spec][0]:
				id = species[spec][0].replace('_str', '')
				node = Node("Species", id=id, name=species[spec][0], compartment='H', hasOnlySubstanceUnits=species[spec][3], boundaryCondition=species[spec][4], constant=species[spec][5])
				graph.create(node)
				relation = Relationship(nodeModel, 'CONTAINS', node)
				graph.create(relation)
			else:
				node = Node("Species", id=species[spec][0], name=species[spec][0], compartment='C', hasOnlySubstanceUnits=species[spec][3], boundaryCondition=species[spec][4], constant=species[spec][5])
				graph.create(node)
				relation = Relationship(nodeModel, 'CONTAINS', node)
				graph.create(relation)

	# nodeModel = Node('Model', id='Rice', name='Rice')
	listUncReac = []

	for reac in reactions.keys():
		exists = False
		if 'mit_' in reactions[reac][0] or 'chl_' in reactions[reac][0]:
			id = reactions[reac][0][4:]
		else:
			id = reactions[reac][0]
		try:
			id = dicReact[id]
			if matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + reactions[reac][1] + "'").first() is None:
				node = Node('Reaction', id=id, compartment=reactions[reac][1])
			else:
				node = matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + reactions[reac][1] + "'").first()
				exists = True
		except:
			listUncReac.append([reac, reactions[reac][0]])
			if matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + reactions[reac][1] + "'").first() is None:
				node = Node('Reaction', id=id, compartment=reactions[reac][1])
			else:
				node = matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + reactions[reac][1] + "'").first()
				exists = True
		reversible = True

		if not exists:
			for j in range(len(reactions[reac])):
				if reactions[reac][j] == '|':
					break
				if j == 0:
					node['name'] = reactions[reac][j]
				elif j == 2:
					node['reversible'] = reactions[reac][j]
				elif j == 3:
					node['fast'] = reactions[reac][j]
			try:
				if node['reversible'].lower() == 'false': reversible = False
			except:
				pass
			graph.create(node)

			for i in range(j + 1, len(reactions[reac])):
				try:
					if '_mit' in species[reactions[reac][i]][0] or 'mit_' in species[reactions[reac][i]][0]:
						try:
							id = species[reactions[reac][i]][0].replace('_mit', '')
							compart = 'M'
						except:
							id = species[reactions[reac][i]][0].replace('mit_', '')
							compart = 'M'
					elif 'x_' in species[reactions[reac][i]][0] or '_x' in species[reactions[reac][i]][0]:
						id = species[reactions[reac][i]][0].replace('x_', '')
						compart = 'E'
					elif 'str_' in species[reactions[reac][i]][0] or '_str' in species[reactions[reac][i]][0]:
						id = species[reactions[reac][i]][0].replace('_str', '')
						compart = 'H'
					else:
						id = species[reactions[reac][i]][0]
						compart = 'C'
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
					if '_mit' in species[reactions[reac][i]][0] or 'mit_' in species[reactions[reac][i]][0]:
						try:
							id = species[reactions[reac][i]][0].replace('_mit', '')
							compart = 'M'
						except:
							id = species[reactions[reac][i]][0].replace('mit_', '')
							compart = 'M'
					elif 'x_' in species[reactions[reac][i]][0] or '_x' in species[reactions[reac][i]][0]:
						id = species[reactions[reac][i]][0].replace('x_', '')
						compart = 'E'
					elif 'str_' in species[reactions[reac][i]][0] or '_str' in species[reactions[reac][i]][0]:
						id = species[reactions[reac][i]][0].replace('_str', '')
						compart = 'H'
					else:
						id = species[reactions[reac][i]][0]
						compart = 'C'
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
		relation = Relationship(nodeModel, 'CONTAINS', node)
		graph.create(relation)

	import csv

	with open('UnconvertedReactionsRice.txt', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(listUncReac)


'''
listDone = []
with open('listIDs_Rice.txt', 'w') as outfile:
	for name in species.keys():
		if '_ext' in species[name][0] or '_acc' in species[name][0] or '_pmd' in species[name][0]:
			if species[name][0][:-4] not in listDone:
				outfile.write(species[name][0][:-4] + '\n')
				listDone.append(species[name][0][:-4])
		elif '_c' in species[name][0] or '_m' in species[name][0] or '_p' in species[name][0] or '_v' in species[name][0] or '_x' in species[name][0]:
			if species[name][0][:-2] not in listDone:
				outfile.write(species[name][0][:-2] + '\n')
				listDone.append(species[name][0][:-2])
		elif '_biomass' in species[name][0]:
			if species[name][0][:-8] not in listDone:
				outfile.write(species[name][0][:-8] + '\n')
				listDone.append(species[name][0][:-8])
		elif '_leaf' in species[name][0] or '_accL' in species[name][0]:
			if species[name][0][:-5] not in listDone:
				outfile.write(species[name][0][:-5] + '\n')
				listDone.append(species[name][0][:-5])
		elif '_total' in species[name][0]:
			if species[name][0][:-6] not in listDone:
				outfile.write(species[name][0][:-6] + '\n')
				listDone.append(species[name][0][:-6])
		else:
			if species[name][0] not in listDone:
				outfile.write(species[name][0] + '\n')
				listDone.append(species[name][0])
				'''