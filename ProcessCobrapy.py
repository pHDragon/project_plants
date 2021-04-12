import reframed
from reframed import load_cbmodel
from ProcessXML import ProcessXML
from os.path import join


class ProcessCobrapy:

	def __init__(self, model):
		try:
			self.model = load_cbmodel("C:\\Users\\Pedro\\OneDrive\\Documentos\\UMinho\\WholeNewProject\\PythonThings\\Models\\" + model)
			self.useReframed = True
		except:
			self.model = ProcessXML(model)
			self.model.getCompartments()
			self.model.getSpecies()
			self.model.getReactions()
			self.useReframed = False

if __name__ == '__main__':
	test = ProcessCobrapy('MedicagoTruncatula.xml')
