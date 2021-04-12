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
				if i.attrib['id'][-5] == 'c':
					self.species[i.attrib['id']].append('cytosol')
				elif i.attrib['id'][-5] == 'm':
					self.species[i.attrib['id']].append('mitochondria')
				elif i.attrib['id'][-5] == 'p':
					self.species[i.attrib['id']].append('plastid')
				elif i.attrib['id'][-5] == 'x':
					self.species[i.attrib['id']].append('peroxisome')
				elif i.attrib['id'][-5] == 'v':
					self.species[i.attrib['id']].append('vacuole')
				self.species[i.attrib['id']].append(i.attrib['id'][-3:])
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
				self.reactions[i.attrib['id']] = [i.attrib['id']]
			if '_c' in i.attrib['id']:
				self.reactions[i.attrib['id']].append('C')
			elif '_m' in i.attrib['id']:
				self.reactions[i.attrib['id']].append('M')
			elif '_p' in i.attrib['id']:
				self.reactions[i.attrib['id']].append('P')
			elif '_x' in i.attrib['id']:
				self.reactions[i.attrib['id']].append('X')
			elif '_v' in i.attrib['id']:
				self.reactions[i.attrib['id']].append('V')
			elif '_mc' in id or '_pc' in id or '_xc' in id or '_tx' in id or '_vc' in id:
				self.reactions[i.attrib['id']].append(i.attrib['id'][-6])
			else:
				self.reactions[i.attrib['id']].append('C')
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
					pass
					# print(child, 'not included in reaction', i.attrib['id'])
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
	data = ProcessXML('GlycineSeedling.xml')

	graph = Graph('bolt://localhost:7687', auth=('neo4j', '123'))
	compartments = data.getCompartments()
	species = data.getSpecies()
	reactions = data.getReactions()
	dicReact = data.getReactionDict()

	'''	nodeC = Node("Compartment", id='C', name='Cytosol', constant=True, model='Glycine max')
	graph.create(nodeC)
	nodeM = Node("Compartment", id='M', name='Mitochondrion', constant=True, model='Glycine max')
	graph.create(nodeM)
	nodeP = Node("Compartment", id='P', name='Plastid', constant=True, model='Glycine max')
	graph.create(nodeP)
	nodeX = Node("Compartment", id='X', name='Peroxisome', constant=True, model='Glycine max')
	graph.create(nodeX)
	nodeV = Node("Compartment", id='V', name='Vacuole', constant=True, model='Glycine max')
	graph.create(nodeV)'''

	tissueCOT = Node("Tissue", id='COT', name='Cotyledons', model='Glycine max')
	graph.create(tissueCOT)
	tissueHRA = Node("Tissue", id='HRA', name='Hypocotyl/root axis', model='Glycine max')
	graph.create(tissueHRA)

	nodeModel = Node('Model', id='Glycine max', name='Glycine max')
	matcher = NodeMatcher(graph)

	for spec in species.keys():
		id = spec[2:]
		if '_COT' in spec:
			id = id.replace('_COT', '')
			if '_c' in spec:
				id = id.replace('_c', '')
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'C' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='C', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', node)
					graph.create(relation)
			elif '_m' in spec:
				id = id.replace('_m', '')
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'M' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='M', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', node)
					graph.create(relation)
			elif '_p' in spec:
				id = id.replace('_p', '')
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'P' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='P', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', node)
					graph.create(relation)
			elif '_x' in spec:
				id = id.replace('_x', '')
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'X' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='X', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', node)
					graph.create(relation)
			elif '_v' in spec:
				id = id.replace('_v', '')
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'V' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='V', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', node)
					graph.create(relation)
			else:
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'C' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='C', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueCOT, 'CONTAINS', node)
					graph.create(relation)
		elif '_HRA' in spec:
			id = id.replace('_HRA', '')
			if '_c' in spec:
				id = id.replace('_c', '')
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'C' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='C', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', node)
					graph.create(relation)
			elif '_m' in spec:
				id = id.replace('_m', '')
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'M' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='M', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', node)
					graph.create(relation)
			elif '_p' in spec:
				id = id.replace('_p', '')
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'P' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='P', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', node)
					graph.create(relation)
			elif '_x' in spec:
				id = id.replace('_x', '')
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'X' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='X', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', node)
					graph.create(relation)
			elif '_v' in spec:
				id = id.replace('_v', '')
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'V' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='V', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', node)
					graph.create(relation)
			else:
				try:
					nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + 'C' + '"')
					relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', nodeSpec.first())
					graph.create(relation)
				except:
					node = Node("Species", id=id, name=species[spec][0], compartment='C', boundaryCondition=species[spec][3], constant=species[spec][4])
					graph.create(node)
					relation = Relationship(nodeModel, 'CONTAINS', node)
					graph.create(relation)
					relation = Relationship(tissueHRA, 'CONTAINS', node)
					graph.create(relation)

	# nodeModel = Node('Model', id='Glycine max', name='Glycine max')
	listUncReac = []

	for reac in reactions.keys():
		exists = False
		id = reac.replace('R_', '')
		compart = ''  # reac[-5].upper()
		tissue = id[-3:]
		id = id[:-4]
		if '_mc' in id or '_pc' in id or '_xc' in id or '_tx' in id or '_vc' in id:
			id = id[:-3]
		elif '_c' in id or '_m' in id or '_p' in id or '_x' in id or '_v' in id:
			id = id[:-2]
		try:
			id = dicReact[id]
		except:
			try:
				id = dicReact['R_' + id]
			except:
				listUncReac.append([reac, reactions[reac][0]])
		if matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + reactions[reac][1] + "'").first() is None:
			node = Node('Reaction', id=id, compartment=reactions[reac][1])
			node['name'] = reactions[reac][0]
			node['reversible'] = reactions[reac][2]
			node['fast'] = reactions[reac][3]
			graph.create(node)
		else:
			node = matcher.match('Reaction').where("_.id='" + id + "'").where("_.compartment='" + reactions[reac][1] + "'").first()
			exists = True
		if tissue == 'COT':
			relation = Relationship(tissueCOT, 'CONTAINS', node)
			graph.create(relation)
		else:
			relation = Relationship(tissueHRA, 'CONTAINS', node)
			graph.create(relation)

		reversible = True
		try:
			if node['reversible'].lower() == 'false': reversible = False
		except:
			pass

		if not exists:
			for i in range(reactions[reac].index('|')+1, len(reactions[reac])):
				try:
					id = reactions[reac][i][2:-4]
					if '_c' in id:
						id = id[:-2]
						compart = 'C'
					elif '_m' in id:
						id = id[:-2]
						compart = 'M'
					elif '_p' in id:
						id = id[:-2]
						compart = 'P'
					elif '_x' in id:
						id = id[:-2]
						compart = 'X'
					elif '_v' in id:
						id = id[:-2]
						compart = 'V'
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
					id = reactions[reac][i][2:-4]
					if '_c' in id:
						id = id[:-2]
						compart = 'C'
					elif '_m' in id:
						id = id[:-2]
						compart = 'M'
					elif '_p' in id:
						id = id[:-2]
						compart = 'P'
					elif '_x' in id:
						id = id[:-2]
						compart = 'X'
					elif '_v' in id:
						id = id[:-2]
						compart = 'V'
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

	import csv

	with open('UnconvertedReactionsSoybean.txt', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(listUncReac)
