import xml.etree.ElementTree as ET


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
			try:
				self.species[i.attrib['id']].append(i.attrib['boundaryCondition'])
			except:
				pass
		return self.species

	def getReactions(self):
		self.reactions = {}
		for i in self.listReactions:
			try:
				self.reactions[i.attrib['id']] = [i.attrib['name'], i.attrib['reversible'], '|']
				self.__tryRecur(i, [])
			except:
				try:
					self.reactions[i.attrib['id']] = [i.attrib['name'], '|']
					self.__tryRecur(i, [])
				except:
					self.reactions[i.attrib['id']] = ['|']
					self.__tryRecur(i, [])
		return self.reactions

	def __tryRecur(self, react, lista):
		for i in react:
			a = self.__tryRecur(i, lista)
			try:
				for j in a:
					self.reactions[react.attrib['id']].append(j)
				lista = []
			except:
				pass
		if react.attrib == {}:
			if 'reactant' in react.tag.lower():
				lista.append('|')
		else:
			for child in react.attrib:
				if 'urn:' not in child or '#_' not in child:
					lista.append(react.attrib[child])
		for i in lista:
			if 'urn:' in i or '#_' in i:
				lista.remove(i)
		return lista


if __name__ == '__main__':
	test = ProcessXML('Ara2009.xml')
	test.getCompartments()
	test.getReactions()
	test.getSpecies()
	print(test.species)
	with open('nomes_linhas.txt', 'w') as f:
		for i in list(test.species.keys()):
			f.write(test.species[i][0]+'\n')

