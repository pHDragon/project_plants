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
				print(i)

	def getCompartments(self):
		self.compartments = {}
		for i in self.listCompartments:
			self.compartments[i.attrib['id']] = i.attrib['name']
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
			except:
				pass
		return self.species

	def getReactions(self):
		self.reactions = {}
		for i in self.listReactions:
			try:
				self.reactions[i.attrib['name']] = [i.attrib['reversible'], i.attrib['fast'], '|']
			except:
				self.reactions[i.attrib['name']] = ['true', i.attrib['fast'], '|']
			for child in i:
				if 'listOfReactants' in child.tag:
					for reac in child:
						self.reactions[i.attrib['name']].append(self.species[reac.attrib['species']][0])
						try:
							self.reactions[i.attrib['name']].append(float(reac.attrib['stoichiometry']))
						except:
							self.reactions[i.attrib['name']].append(1.0)
				elif 'listOfProducts' in child.tag:
					self.reactions[i.attrib['name']].append('|')
					for prod in child:
						self.reactions[i.attrib['name']].append(self.species[prod.attrib['species']][0])
						try:
							self.reactions[i.attrib['name']].append(float(prod.attrib['stoichiometry']))
						except:
							self.reactions[i.attrib['name']].append(1.0)
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
	data = ProcessXML('Ara2009.xml')
	graph = Graph('bolt://palsson.di.uminho.pt:6092', auth=('neo4j', '123'))
	compartments = data.getCompartments()
	species = data.getSpecies()
	reactions = data.getReactions()
	dicReact = data.getReactionDict()

	'''	nodeC = Node("Compartment", id='C', name='Cytosol', model='Arabidopsis')
	graph.create(nodeC)
	nodeM = Node("Compartment", id='M', name='Mitochondrion', model='Arabidopsis')
	graph.create(nodeM)
	nodeE = Node("Compartment", id='E', name='Extracellular', model='Arabidopsis')
	graph.create(nodeE)'''

	matcher = NodeMatcher(graph)
	nodeModel = Node('Model', id='Arabidopsis', name='Arabidopsis')

	for spec in species.keys():
		try:
			if '_mit' in species[spec][0]:
				id = species[spec][0].replace('_mit', '')
				nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'M' + '"')
				relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
				graph.create(relation)
			elif 'x_' in species[spec][0]:
				id = species[spec][0].replace('x_', '')
				nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'E' + '"')
				relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
				graph.create(relation)
			else:
				id = species[spec][0]
				nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'C' + '"')
				relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
				graph.create(relation)
		except:
			if '_mit' in species[spec][0]:
				id = species[spec][0].replace('_mit', '')
				node = Node("Species", id=id, name=species[spec][0], compartment='M')
				graph.create(node)
				relation = Relationship(nodeModel, 'CONTAINS', node)
				graph.create(relation)
			elif 'x_' in species[spec][0]:
				id = species[spec][0].replace('x_', '')
				node = Node("Species", id=id, name=species[spec][0], compartment='E')
				graph.create(node)
				relation = Relationship(nodeModel, 'CONTAINS', node)
				graph.create(relation)
			else:
				id = species[spec][0]
				node = Node("Species", id=id, name=species[spec][0], compartment='C')
				graph.create(node)
				relation = Relationship(nodeModel, 'CONTAINS', node)
				graph.create(relation)

	# nodeModel = Node('Model', id='Arabidopsis', name='Arabidopsis')
	listUncReac = []

	for reac in reactions.keys():
		exists = False
		id = reac
		try:
			id = dicReact[id]
			if matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + 'C' + "'").first() is None:
				node = Node('Reaction', id=id, name=reac, compartment='C')
			else:
				node = matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + 'C' + "'").first()
				exists = True
		except:
			listUncReac.append([reac, ''])
			if matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + 'C' + "'").first() is None:
				node = Node('Reaction', id=id, name=reac, compartment='C')
			else:
				node = matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + 'C' + "'").first()
				exists = True
		reversible = True

		if not exists:
			for j in range(len(reactions[reac])):
				if reactions[reac][j] == '|':
					break
				if j == 0:
					node['reversible'] = reactions[reac][j]
				elif j == 1:
					node['fast'] = reactions[reac][j]
			try:
				if node['reversible'].lower() == 'false': reversible = False
			except:
				pass
			graph.create(node)

			for i in range(j + 1, len(reactions[reac])):
				try:
					id = reactions[reac][i]
					if '_mit' in id:
						id = id.replace('_mit', '')
						compart = 'M'
					elif 'x_' in id:
						id = id.replace('x_', '')
						compart = 'E'
					else:
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
					id = reactions[reac][i]
					if '_mit' in id:
						id = id.replace('_mit', '')
						compart = 'M'
					elif 'x_' in id:
						id = id.replace('x_', '')
						compart = 'E'
					else:
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

	with open('UnconvertedReactionsArabidopsis.txt', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(listUncReac)
