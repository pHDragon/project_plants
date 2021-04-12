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
			self.compartments[i.attrib['id']] = [i.attrib['constant'], i.attrib['size'], i.attrib['spatialDimensions']]
		return self.compartments

	def getSpecies(self):
		self.species = {}
		for i in self.listSpecies:
			try:
				self.species[i.attrib['id']] = [i.attrib['name']]
			except:
				self.species[i.attrib['id']] = [i.attrib['id']]
			try:
				self.species[i.attrib['id']].append(i.attrib['boundaryCondition'])
				self.species[i.attrib['id']].append(i.attrib['compartment'])
				self.species[i.attrib['id']].append(i.attrib['constant'])
				self.species[i.attrib['id']].append(i.attrib['hasOnlySubstanceUnits'])
				self.species[i.attrib['id']].append(i.attrib['initialAmount'])
				self.species[i.attrib['id']].append(i.attrib['metaid'])
				self.species[i.attrib['id']].append(i[0][0][0][0][0][0].attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource'][-5:])
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
				self.reactions[i.attrib['id']].append(i.attrib['compartment'])
				self.reactions[i.attrib['id']].append(i.attrib['fast'])
				self.reactions[i.attrib['id']].append(i.attrib['metaid'])
				self.reactions[i.attrib['id']].append(i.attrib['reversible'])
				self.reactions[i.attrib['id']].append(i.attrib['sboTerm'])
				self.reactions[i.attrib['id']].append('|')
			except:
				self.reactions[i.attrib['id']].append('|')
				pass
			for child in i:
				if 'listOfReactants' in child.tag:
					for reac in child:
						self.reactions[i.attrib['id']].append(reac.attrib['species'])
						try:
							self.reactions[i.attrib['id']].append(reac.attrib['constant'])
							self.reactions[i.attrib['id']].append(float(reac.attrib['stoichiometry']))
						except:
							self.reactions[i.attrib['id']].append(reac.attrib['constant'])
							self.reactions[i.attrib['id']].append(1.0)

				elif 'listOfProducts' in child.tag:
					self.reactions[i.attrib['id']].append('|')
					for prod in child:
						self.reactions[i.attrib['id']].append(prod.attrib['species'])
						try:
							self.reactions[i.attrib['id']].append(prod.attrib['constant'])
							self.reactions[i.attrib['id']].append(float(prod.attrib['stoichiometry']))

						except:
							self.reactions[i.attrib['id']].append(prod.attrib['constant'])
							self.reactions[i.attrib['id']].append(1.0)

				else:
					print(child, 'not included in reaction', i.attrib['id'])
		return self.reactions


if __name__ == '__main__':
	data = ProcessXML('Arabidopsis2017.xml')

	graph = Graph('bolt://localhost:7687', auth=('neo4j', '123'))
	compartments = data.getCompartments()
	species = data.getSpecies()
	reactions = data.getReactions()

	for spec in species.keys():
		node = Node("Species", id=spec, name=species[spec][0], boundaryCondition=species[spec][1], compartment=species[spec][2], constant=species[spec][3], hasOnlySubstanceUnits=species[spec][4], initialAmount=species[spec][5], metaid=species[spec][6])
		try:
			node["CHEBI"] = species[spec][7]
		except: pass
		graph.create(node)

	matcher = NodeMatcher(graph)

	for reac in reactions.keys():
		node = Node("Reaction", id=reac)
		reversible = True

		node['name'] = reactions[reac][0]
		node['compartment'] = reactions[reac][1]
		node['fast'] = reactions[reac][2]
		node['metaid'] = reactions[reac][3]
		node['reversible'] = reactions[reac][4]
		node['sboTerm'] = reactions[reac][5]

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
			elif matcher.match('Species').where("_.id='" + reactions[reac][i] + "'").first() is None:
				pass
			else:
				nodeSpec = matcher.match('Species').where("_.id='" + reactions[reac][i] + "'")
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
			elif matcher.match('Species').where("_.id='" + reactions[reac][i] + "'").first() is None:
				pass
			else:
				nodeSpec = matcher.match('Species').where("_.id='" + reactions[reac][i] + "'")
				relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=1.0)
				graph.create(relation)
				if reversible:
					relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=1.0)
					graph.create(relation)
