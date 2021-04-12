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
			self.compartments[i.attrib['id']] = [i.attrib['name'], i.attrib['constant']]
		return self.compartments

	def getSpecies(self):
		self.species = {}
		for i in self.listSpecies:
			try:
				self.species[i.attrib['id']] = [i.attrib['name']]
			except:
				self.species[i.attrib['id']] = [i.attrib['id'][:-2]]
			try:
				if i.attrib['id'][-1] == 'c':
					self.species[i.attrib['id']].append('cytosol')
				elif i.attrib['id'][-1] == 'm':
					self.species[i.attrib['id']].append('mitochondria')
				elif i.attrib['id'][-1] == 'p':
					self.species[i.attrib['id']].append('plastid')
				elif i.attrib['id'][-1] == 'x':
					self.species[i.attrib['id']].append('peroxisome')
				elif i.attrib['id'][-1] == 'v':
					self.species[i.attrib['id']].append('vacuole')
				self.species[i.attrib['id']].append(i.attrib['boundaryCondition'])
				self.species[i.attrib['id']].append(i.attrib['constant'])
				self.species[i.attrib['id']].append(i.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}charge'])
				self.species[i.attrib['id']].append(i.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}chemicalFormula'])
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
				self.reactions[i.attrib['id']].append(i.attrib['fast'])
				if i.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}lowerFluxBound'] == 'cobra_0_bound':
					self.reactions[i.attrib['id']].append(0)
				elif i.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}lowerFluxBound'] == 'cobra_default_lb':
					self.reactions[i.attrib['id']].append(-1000)
				if i.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}upperFluxBound'] == 'cobra_0_bound':
					self.reactions[i.attrib['id']].append(0)
				elif i.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}upperFluxBound'] == 'cobra_default_ub':
					self.reactions[i.attrib['id']].append(1000)
				self.reactions[i.attrib['id']].append('|')
			except:
				self.reactions[i.attrib['id']].append('true')
				self.reactions[i.attrib['id']].append(i.attrib['fast'])
				if i.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}lowerFluxBound'] == 'cobra_0_bound':
					self.reactions[i.attrib['id']].append(0)
				elif i.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}lowerFluxBound'] == 'cobra_default_lb':
					self.reactions[i.attrib['id']].append(-1000)
				if i.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}upperFluxBound'] == 'cobra_0_bound':
					self.reactions[i.attrib['id']].append(0)
				elif i.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}upperFluxBound'] == 'cobra_default_ub':
					self.reactions[i.attrib['id']].append(1000)
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


if __name__ == '__main__':
	data = ProcessXML('GlycineNormal.xml')

	graph = Graph('bolt://localhost:7687', auth=('neo4j', '123'))
	compartments = data.getCompartments()
	species = data.getSpecies()
	reactions = data.getReactions()

	for spec in species.keys():
		node = Node("Species", id=spec, name=species[spec][0], compartment=species[spec][1], boundaryCondition=species[spec][2], constant=species[spec][3])
		try:
			node['charge'] = species[spec][4]
			node['chemicalFormula'] = species[spec][5]
		except: pass
		graph.create(node)

	matcher = NodeMatcher(graph)

	for reac in reactions.keys():
		node = Node("Reaction", id=reac)
		reversible = True

		node['name'] = reactions[reac][0]
		node['reversible'] = reactions[reac][1]
		node['fast'] = reactions[reac][2]
		node['lowerBound'] = reactions[reac][3]
		node['upperBound'] = reactions[reac][4]

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
