import xml.etree.ElementTree as ET
from py2neo import Graph, Node, Relationship, NodeMatcher, Subgraph
import re


def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

def processASCII(string):
	import re
	findstr = re.findall('__[0-9]{2}__', string)
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
			elif 'flux' in i.tag.lower():
				self.fluxBounds = i
			else:
				print(i,'not included in model')

	def getCompartments(self):
		self.compartments = {}
		for i in self.listCompartments:
			self.compartments[i.attrib['id']] = [i.attrib['name'], i.attrib['spatialDimensions'], i.attrib['size'], i.attrib['constant'], 'Medicago truncatula']
		return self.compartments

	def getSpecies(self):
		self.species = {}
		for i in self.listSpecies:
			try:
				self.species[processASCII(i.attrib['id'])] = [i.attrib['name']]
			except:
				self.species[processASCII(i.attrib['id'])] = [i.attrib['id'][:-2]]
			try:
				self.species[processASCII(i.attrib['id'])].append(i.attrib['compartment'])
				self.species[processASCII(i.attrib['id'])].append(i.attrib['initialAmount'])
				self.species[processASCII(i.attrib['id'])].append(i.attrib['boundaryCondition'])
				self.species[processASCII(i.attrib['id'])].append(i.attrib['constant'])
			except: pass
			try:
				for ids in i[1][0][0][0][0]:  # vai buscar os atributos rdf:resource
					id = list(ids.attrib.values())[0]
					self.species[processASCII(i.attrib['id'])].append(id[11:])
			except: pass
		return self.species

	def getReactions(self):
		self.reactions = {}
		for i in self.listReactions:
			if processASCII(i.attrib['id'])[0] == '_':
				id = processASCII(i.attrib['id'])[1:]
			else:
				id = processASCII(i.attrib['id'])
			try:
				self.reactions[id] = [i.attrib['name']]
			except:
				self.reactions[id] = []
			try:
				self.reactions[id].append(i.attrib['compartment'])
				self.reactions[id].append(i.attrib['reversible'])
				self.reactions[id].append(i.attrib['fast'])
				self.reactions[id].append('|')
			except:
				self.reactions[id].append(i.attrib['compartment'])
				self.reactions[id].append('true')
				self.reactions[id].append(i.attrib['fast'])
				self.reactions[id].append('|')
			for child in i:
				if 'listOfReactants' in child.tag:
					for reac in child:
						if processASCII(reac.attrib['species'])[0] == '_':
							ID = processASCII(reac.attrib['species'])[1:]
						else:
							ID = processASCII(reac.attrib['species'])
						self.reactions[id].append(ID)
						try:
							self.reactions[id].append(float(reac.attrib['stoichiometry']))
							self.reactions[id].append(reac.attrib['constant'])
						except:
							self.reactions[id].append(1.0)
							self.reactions[id].append(reac.attrib['constant'])
				elif 'listOfProducts' in child.tag:
					self.reactions[id].append('|')
					for prod in child:
						if processASCII(prod.attrib['species'])[0] == '_':
							ID = processASCII(prod.attrib['species'])[1:]
						else:
							ID = processASCII(prod.attrib['species'])
						self.reactions[id].append(ID)
						try:
							self.reactions[id].append(float(prod.attrib['stoichiometry']))
							self.reactions[id].append(prod.attrib['constant'])
						except:
							self.reactions[id].append(1.0)
							self.reactions[id].append(prod.attrib['constant'])
				elif 'annotation' in child.tag:
					for ids in child[0][0][0][0]:
						ID = list(ids.attrib.values())[0][11:]
						self.reactions[id].insert(self.reactions[id].index('|'), ID)
				else:
					print(child, 'not included in reaction', i.attrib['id'])
		self.addFluxBounds()
		return self.reactions

	def addFluxBounds(self):
		fbc = '{http://www.sbml.org/sbml/level3/version1/fbc/version1}'
		for i in self.fluxBounds:
			try:
				self.reactions[processASCII(i.attrib[fbc+'reaction'])].insert(self.reactions[processASCII(i.attrib[fbc+'reaction'])].index('|'), i.attrib[fbc+'value'])
			except: pass

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
	data = ProcessXML('MedicagoTruncatula.xml')

	graph = Graph('bolt://localhost:7687', auth=('neo4j', '123'))
	compartments = data.getCompartments()
	species = data.getSpecies()
	reactions = data.getReactions()
	dicReact = data.getReactionDict()

	#for comp in compartments.keys():
		#node = Node("Compartment", id=comp, name=compartments[comp][0], spatialDimensions=compartments[comp][1], size=compartments[comp][2], constant=compartments[comp][3], model=compartments[comp][4])
		#graph.create(node)

	nodeModel = Node('Model', id='Medicago truncatula', name='Medicago truncatula')
	matcher = NodeMatcher(graph)

	for spec in species.keys():
		id = spec
		if spec[0] == '_':
			id = spec[1:]
		if '_' in id:
			id = id[:id.rindex('_')]
		compart = spec[spec.rindex('_')+1:]
		nodeSpec = matcher.match('Species').where('_.id="' + id + '"').where('_.compartment="' + compart + '"')
		if nodeSpec.first() is None:
			node = Node("Species", id=id, name=species[spec][0], compartment=species[spec][1], charge=species[spec][2], boundaryCondition=species[spec][3], constant=species[spec][4])
			for i in range(5, len(species[spec])):
				id = ''
				try:
					ids = species[spec][i].split(':')
					if id == ids[0]:
						try:
							node[ids[0]].append(ids[1][ids[1].index('%3A') + 3:])
						except:
							node[ids[0]].append(ids[1])
						id = ids[0]
					else:
						try:
							node[ids[0]] = [ids[1][ids[1].index('%3A') + 3:]]
						except:
							node[ids[0]] = [ids[1]]
						id = ids[0]
				except: pass
			graph.create(node)
			#nodeComp = matcher.match('Compartment').where('_.id="' + spec[spec.rindex('_')+1:] + '"').where('_.model="' + 'Medicago truncatula' + '"')
			relation = Relationship(nodeModel, 'CONTAINS', node)
			graph.create(relation)
		else:
			#nodeComp = matcher.match('Compartment').where('_.id="' + spec[spec.rindex('_') + 1:] + '"').where('_.model="' + 'Medicago truncatula' + '"')
			relation = Relationship(nodeModel, 'CONTAINS', nodeSpec.first())
			graph.create(relation)

	#nodeModel = Node('Model', id='Medicago truncatula', name='Medicago truncatula')

	listUncReac = []

	for reac in reactions.keys():
		try:
			pos = re.search('_[V|M|H|L|G|X|C|R|E|tx]{1}', reac[-3:]).start()
			id = reac[:pos-3]
		except:
			id = reac
		try:
			node = Node("Reaction", id=dicReact[id], compartment=reactions[reac][1])
		except:
			node = Node("Reaction", id=id, compartment=reactions[reac][1])
			listUncReac.append([reac, reactions[reac][0]])
		reversible = True

		node['name'] = reactions[reac][0]
		node['reversible'] = reactions[reac][2]
		node['fast'] = reactions[reac][3]
		try:
			ids = reactions[reac][4].split(':')
			node[ids[0]] = ids[1][ids[1].index('%3A')+3:]
			lb = float(reactions[reac][5])
			ub = float(reactions[reac][6])
			relation = Relationship(nodeModel, 'CONTAINS', node, lowerBound=lb, upperBound=ub)
			graph.create(relation)
		except:
			try:
				lb = float(reactions[reac][4])
				ub = float(reactions[reac][5])
				relation = Relationship(nodeModel, 'CONTAINS', node, lowerBound=lb, upperBound=ub)
				graph.create(relation)
			except:
				relation = Relationship(nodeModel, 'CONTAINS', node)
				graph.create(relation)

		try:
			if node['reversible'].lower() == 'false': reversible = False
		except:
			pass
		graph.create(node)

		for i in range(reactions[reac].index('|')+1, len(reactions[reac])):
			try:
				spec = reactions[reac][i][:reactions[reac][i].rindex('_')]
				compart = reactions[reac][i][reactions[reac][i].rindex('_')+1:]
			except:
				spec = ''
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
			elif matcher.match('Species').where('_.id="' + spec + '"').where('_.compartment="' + compart + '"').first() is None:
				pass
			else:
				nodeSpec = matcher.match('Species').where('_.id="' + spec + '"').where('_.compartment="' + compart + '"')
				relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=1.0)
				graph.create(relation)
				if reversible:
					relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=1.0)
					graph.create(relation)

		for i in range(j + 1, len(reactions[reac])):
			try:
				spec = reactions[reac][i][:reactions[reac][i].rindex('_')]
				compart = reactions[reac][i][reactions[reac][i].rindex('_') + 1:]
			except:
				spec = ''
				compart = ''
			if isfloat(reactions[reac][i]):
				relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=float(reactions[reac][i]))
				graph.create(relation)
				if reversible:
					relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=float(reactions[reac][i]))
					graph.create(relation)
			elif reactions[reac][i] == 'false' or reactions[reac][i] == 'true':
				pass
			elif matcher.match('Species').where('_.id="' + spec + '"').where('_.compartment="' + compart + '"').first() is None:
				pass
			else:
				nodeSpec = matcher.match('Species').where('_.id="' + spec + '"').where('_.compartment="' + compart + '"')
				relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=1.0)
				graph.create(relation)
				if reversible:
					relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=1.0)
					graph.create(relation)
	import csv

	with open('UnconvertedReactionsMedicago.txt', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(listUncReac)