from py2neo import Graph, Node, Relationship, NodeMatcher, Subgraph
from ProcessXML import ProcessXML


graph = Graph('bolt://localhost:7687', auth=('neo4j', '123'))
data = ProcessXML('AraMeta.xml')
compartments = data.getCompartments()
species = data.getSpecies()
reactions = data.getReactions()

for spec in species.keys():
	node = Node("Species", id=spec, name=species[spec][0], compartment=species[spec][1])
	graph.create(node)

for comp in compartments.keys():
	node = Node("Compartment", id=comp, name=compartments[comp])
	graph.create(node)


def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False


for reac in reactions.keys():
	matcher = NodeMatcher(graph)
	node = Node("Reaction", id=reac)
	reversible = True

	for j in range(len(reactions[reac])):
		if reactions[reac][j] == '|':
			break
		if j==0: node['name']=reactions[reac][j]
		elif j==1: node['reversible']=reactions[reac][j]
	try:
		if node['reversible'].lower() == 'false': reversible = False
	except: pass
	graph.create(node)

	for i in range(j+1, len(reactions[reac])):
		if reactions[reac][i] == '|':
			j=i
			break
		if isfloat(reactions[reac][i]):
			relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=float(reactions[reac][i]))
			graph.create(relation)
			if reversible:
				relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=float(reactions[reac][i]))
				graph.create(relation)
		elif reactions[reac][i] == 'false' or reactions[reac][i] == 'true':
			pass
		elif matcher.match('Species').where("_.id='" + reactions[reac][i]+"'").first() is None:
			pass
		else:
			nodeSpec = matcher.match('Species').where("_.id='" + reactions[reac][i]+"'")
			relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=1.0)
			graph.create(relation)
			if reversible:
				relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=1.0)
				graph.create(relation)

	for i in range(j+1, len(reactions[reac])):
		if isfloat(reactions[reac][i]):
			relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=float(reactions[reac][i]))
			graph.create(relation)
			if reversible:
				relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=float(reactions[reac][i]))
				graph.create(relation)
		elif reactions[reac][i] == 'false' or reactions[reac][i] == 'true':
			pass
		elif matcher.match('Species').where("_.id='" + reactions[reac][i]+"'").first() is None:
			pass
		else:
			nodeSpec = matcher.match('Species').where("_.id='" + reactions[reac][i] + "'")
			relation = Relationship(node, 'PRODUCES', nodeSpec.first(), stoichiometry=1.0)
			graph.create(relation)
			if reversible:
				relation = Relationship(nodeSpec.first(), 'USES', node, stoichiometry=1.0)
				graph.create(relation)


