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
				print(i,'not included in model')

	def getCompartments(self):
		self.compartments = {}
		for i in self.listCompartments:
			self.compartments[i.attrib['id']] = [i.attrib['name'], i.attrib['size'], i.attrib['constant']]
		return self.compartments

	def getSpecies(self):
		self.species = {}
		for i in self.listSpecies:
			try:
				self.species[i.attrib['id']] = [i.attrib['name']]
			except:
				self.species[i.attrib['id']] = [i.attrib['id']]
			try:
				if '_Cyto' in i.attrib['name']:
					self.species[i.attrib['id']].append('C')
					self.species[i.attrib['id']][0] = self.species[i.attrib['id']][0][:-5]
				elif '_Mito' in i.attrib['name']:
					self.species[i.attrib['id']].append('M')
					self.species[i.attrib['id']][0] = self.species[i.attrib['id']][0][:-5]
				elif '_Plas' in i.attrib['name']:
					self.species[i.attrib['id']].append('P')
					self.species[i.attrib['id']][0] = self.species[i.attrib['id']][0][:-5]
				elif '_Pero' in i.attrib['name']:
					self.species[i.attrib['id']].append('X')
					self.species[i.attrib['id']][0] = self.species[i.attrib['id']][0][:-5]
				elif '_Vacu' in i.attrib['name']:
					self.species[i.attrib['id']].append('V')
					self.species[i.attrib['id']][0] = self.species[i.attrib['id']][0][:-5]
				elif 'x_' in i.attrib['name']:
					self.species[i.attrib['id']].append('E')
					self.species[i.attrib['id']][0] = self.species[i.attrib['id']][0][2:]
				else:
					self.species[i.attrib['id']].append('C')

				self.species[i.attrib['id']].append(float(i.attrib['initialConcentration']))
				self.species[i.attrib['id']].append(i.attrib['hasOnlySubstanceUnits'])
				self.species[i.attrib['id']].append(i.attrib['boundaryCondition'])
				self.species[i.attrib['id']].append(i.attrib['constant'])
			except: pass
		return self.species

	def getReactions(self):
		self.reactions = {}
		for i in self.listReactions:
			try:
				self.reactions[i.attrib['id']] = [i.attrib['name']]
			except:
				self.reactions[i.attrib['id']] = []
			if '_Cyto' in i.attrib['name']:
				self.reactions[i.attrib['id']][0] = self.reactions[i.attrib['id']][0].replace('_Cyto', '')
				self.reactions[i.attrib['id']].append('C')
			elif '_Mito' in i.attrib['name']:
				self.reactions[i.attrib['id']][0] = self.reactions[i.attrib['id']][0].replace('_Mito', '')
				self.reactions[i.attrib['id']].append('M')
			elif '_Plas' in i.attrib['name']:
				self.reactions[i.attrib['id']][0] = self.reactions[i.attrib['id']][0].replace('_Plas', '')
				self.reactions[i.attrib['id']].append('P')
			elif '_Pero' in i.attrib['name']:
				self.reactions[i.attrib['id']][0] = self.reactions[i.attrib['id']][0].replace('_Pero', '')
				self.reactions[i.attrib['id']].append('X')
			elif '_Vacu' in i.attrib['name']:
				self.reactions[i.attrib['id']][0] = self.reactions[i.attrib['id']][0].replace('_Vacu', '')
				self.reactions[i.attrib['id']].append('V')
			else:
				self.reactions[i.attrib['id']].append('C')
			self.reactions[i.attrib['id']].append(i.attrib['reversible'])
			self.reactions[i.attrib['id']].append(i.attrib['fast'])
			self.reactions[i.attrib['id']].append('|')

			for child in i:
				if 'listOfReactants' in child.tag:
					for reac in child:
						self.reactions[i.attrib['id']].append(reac.attrib['species'])
						try:
							self.reactions[i.attrib['id']].append(float(reac.attrib['stoichiometry']))
							self.reactions[i.attrib['id']].append(reac.attrib['constant'])
						except:
							self.reactions[i.attrib['id']].append(1.0)
							self.reactions[i.attrib['id']].append(reac.attrib['constant'])

				elif 'listOfProducts' in child.tag:
					self.reactions[i.attrib['id']].append('|')
					for prod in child:
						self.reactions[i.attrib['id']].append(prod.attrib['species'])
						try:
							self.reactions[i.attrib['id']].append(float(prod.attrib['stoichiometry']))
							self.reactions[i.attrib['id']].append(prod.attrib['constant'])

						except:
							self.reactions[i.attrib['id']].append(1.0)
							self.reactions[i.attrib['id']].append(prod.attrib['constant'])
				else:
					print(child, 'not included in reaction', i.attrib['id'])
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
	data = ProcessXML('Tomato.xml')

	graph = Graph('bolt://localhost:7687', auth=('neo4j', '123'))
	compartments = data.getCompartments()
	species = data.getSpecies()
	reactions = data.getReactions()
	dicReact = data.getReactionDict()

	'''	nodeC = Node("Compartment", id='C', name='Cytosol', constant=True, model='Tomato')
	graph.create(nodeC)
	nodeM = Node("Compartment", id='M', name='Mitochondrion', constant=True, model='Tomato')
	graph.create(nodeM)
	nodeP = Node("Compartment", id='P', name='Plastid', constant=True, model='Tomato')
	graph.create(nodeP)
	nodeX = Node("Compartment", id='X', name='Peroxisome', constant=True, model='Tomato')
	graph.create(nodeX)
	nodeV = Node("Compartment", id='V', name='Vacuole', constant=True, model='Tomato')
	graph.create(nodeV)
	nodeE = Node("Compartment", id='E', name='Extracellular', constant=True, model='Tomato')
	graph.create(nodeE)'''

	matcher = NodeMatcher(graph)
	nodeModel = Node('Model', id='Tomato', name='Tomato')

	for spec in species.keys():
		id = species[spec][0]
		compart = species[spec][1]
		try:
			nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + compart + '"')
			relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
			graph.create(relation)
		except:
			node = Node("Species", id=id, name=species[spec][0], compartment=species[spec][1], hasOnlySubstanceUnits=species[spec][3], boundaryCondition=species[spec][4], constant=species[spec][5])
			graph.create(node)
			relation = Relationship(nodeModel, 'CONTAINS', node)
			graph.create(relation)

	# nodeModel = Node('Model', id='Tomato', name='Tomato')
	listUncReac = []

	for reac in reactions.keys():
		exists = False
		try:
			id = dicReact[reactions[reac][0]]
			if matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + reactions[reac][1] + "'").first() is None:
				node = Node('Reaction', id=id, compartment=reactions[reac][1])
			else:
				node = matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + reactions[reac][1] + "'").first()
				exists = True
		except:
			id = reactions[reac][0]
			listUncReac.append([reac, reactions[reac][0]])
			if matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + reactions[reac][1] + "'").first() is None:
				node = Node('Reaction', id=id, compartment=reactions[reac][1])
			else:
				node = matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + reactions[reac][1] + "'").first()
				exists = True
		reversible = True

		if not exists:
			node['name'] = reactions[reac][0]
			node['reversible'] = reactions[reac][2]
			node['fast'] = reactions[reac][3]

			try:
				if node['reversible'].lower() == 'false': reversible = False
			except:
				pass
			graph.create(node)

			for i in range(reactions[reac].index('|')+1, len(reactions[reac])):
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
				elif matcher.match('Species').where("_.id='" + species[reactions[reac][i]][0] + "'").where("_.compartment='" + species[reactions[reac][i]][1] + "'").first() is None:
					pass
				else:
					nodeSpec = matcher.match('Species').where("_.id='" + species[reactions[reac][i]][0] + "'").where("_.compartment='" + species[reactions[reac][i]][1] + "'")
					relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=1.0)
					graph.create(relation)
					if reversible:
						relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=1.0)
						graph.create(relation)

			for i in range(j + 1, len(reactions[reac])):
				if isfloat(reactions[reac][i]):
					relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=float(reactions[reac][i]))
					graph.create(relation)
					if reversible:
						relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=float(reactions[reac][i]))
						graph.create(relation)
				elif reactions[reac][i] == 'false' or reactions[reac][i] == 'true':
					pass
				elif matcher.match('Species').where("_.id='" + species[reactions[reac][i]][0] + "'").where("_.compartment='" + species[reactions[reac][i]][1] + "'").first() is None:
					pass
				else:
					nodeSpec = matcher.match('Species').where("_.id='" + species[reactions[reac][i]][0] + "'").where("_.compartment='" + species[reactions[reac][i]][1] + "'")
					relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=1.0)
					graph.create(relation)
					if reversible:
						relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=1.0)
						graph.create(relation)
		relation = Relationship(nodeModel, 'CONTAINS', node)
		graph.create(relation)

	import csv

	with open('UnconvertedReactionsTomato.txt', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(listUncReac)