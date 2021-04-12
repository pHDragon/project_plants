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
			self.compartments[i.attrib['id']] = i.attrib['name']
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
				self.species[i.attrib['id']].append(i.attrib['charge'])
			except:
				pass
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

	def getSpeciesDict(self, file='Unique_ModelSEED_Compound_Aliases.txt'):
		import csv
		with open(file, 'r') as csvfile:
			spamreader = csv.reader(csvfile, delimiter='\t')
			dicKEGG = {}
			dic2Meta = {}
			for row in spamreader:
				if row[2] == 'MetaCyc' or row[2] == 'PlantCyc':
					dic2Meta[row[0]] = row[1]
				elif row[2] == 'KEGG':
					dicKEGG[row[1]] = row[0]
		dicSubs = {}
		for i in dicKEGG.keys():
			code = dicKEGG[i]
			try:
				dicSubs[i] = dic2Meta[code]
			except:
				dicSubs[i] = i
		return dicSubs


if __name__ == '__main__':
	data = ProcessXML('Zea mays iRS1563.xml')

	graph = Graph('bolt://localhost:7687', auth=('neo4j', '123'))
	compartments = data.getCompartments()
	species = data.getSpecies()
	reactions = data.getReactions()
	dicReact = data.getReactionDict()
	dicSubs = data.getSpeciesDict()

	'''	for comp in compartments.keys():
		node = Node("Compartment", id=comp.upper(), name=compartments[comp], model='Zea mays')
		graph.create(node)'''

	matcher = NodeMatcher(graph)
	listUncSpec = []
	nodeModel = Node('Model', id='Zea mays', name='Zea mays')

	for spec in species.keys():
		id = spec[:-2]
		compart = species[spec][1].upper()
		try:
			id = dicSubs[id]
		except:
			listUncSpec.append([spec, species[spec][0]])
		try:
			nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + compart + '"')
			relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
			graph.create(relation)
		except:
			node = Node("Species", id=id, name=species[spec][0], compartment=species[spec][1].upper(), charge=int(species[spec][2]))
			graph.create(node)
			relation = Relationship(nodeModel, 'CONTAINS', node)
			graph.create(relation)

	with open('UnconvertedSpeciesZeaMays.txt', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(listUncSpec)


	# nodeModel = Node('Model', id='Zea mays', name='Zea mays')
	listUncReac = []

	for reac in reactions.keys():
		exists = False
		try:
			id = dicReact[reac[:-2]]
		except:
			try:
				id = dicReact['R_' + reac[:-2]]
			except:
				id = reac[:-2]
				listUncReac.append([reac, reactions[reac][0]])

		compart = reac[-1].upper()
		if matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + compart + "'").first() is None:
			node = Node('Reaction', id=id, compartment=compart)
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
					id = dicSubs[reactions[reac][i][:-2]]
					compart = reactions[reac][i][-1].upper()
				except:
					try:
						id = reactions[reac][i][:-2]
						compart = reactions[reac][i][-1].upper()
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
					id = dicSubs[reactions[reac][i][:-2]]
					compart = reactions[reac][i][-1].upper()
				except:
					try:
						id = reactions[reac][i][:-2]
						compart = reactions[reac][i][-1].upper()
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

	with open('UnconvertedReactionsZeaMays.txt', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(listUncReac)
