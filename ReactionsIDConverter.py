import re

#from Utilities import FileUtilities
import IDConverter as IDConverter
#from definitions import TOOL_CONFIG_PATH,ROOT_DIR

class ReactionsIDConverter(IDConverter):

    def __init__(self):
        super().__init__()
        self.__configs = FileUtilities.read_conf_file(TOOL_CONFIG_PATH)
        self.__home_path__ = 'C:\\Users\\Pedro\\OneDrive\\Documentos\\UMinho\\WholeNewProject\\PythonThings'
        self.construct_reactions_converter()


    def convert_modelSeedId_into_other_dbID_reaction(self,modelSeedId,database_name):
        return self.modelSeedToDb.get(modelSeedId).get(database_name)

    def convert_dbID_into_modelSeedId_reaction(self,database_name,id):
        if id in self.dbToModelidSeed.get(database_name).keys():
            return self.dbToModelSeed.get(database_name).get(id)
        else:
            return None

    def construct_reactions_converter(self):

        # try:
        #     info = urllib.request.urlopen("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/dev/Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt")
        #     url = True
        # except:
        info = open(self.__home_path__ + self.__configs.get("reaction_converter"), "r")
        extra_info = open(self.__home_path__ + self.__configs.get("reaction_converter_extra"), "r")
        url=False

        previousModelSeed = ""
        previousOtherDb = ""
        new_dic={}
        for line in info:

            if url:
                line=line.decode("utf-8")

            l = line.split("\t")


            modelSeedId = l[0]
            externalId = l[1]
            source = l[2].strip()

            if re.search("\AModelSEED", line) == None:

                if previousModelSeed != modelSeedId:
                    self.modelSeedToDb[previousModelSeed] = new_dic
                    new_dic = {}
                    new_dic[source]=[externalId]
                    previousModelSeed = modelSeedId
                    previousOtherDb = source

                elif previousOtherDb != source:
                    new_dic[source]=[externalId]
                    previousOtherDb=source

                else:
                    new_dic[source].append(externalId)

                if source in self.dbToModelSeed.keys():
                    if externalId in self.dbToModelSeed.get(source).keys():
                        self.dbToModelSeed[source][externalId].append(modelSeedId)
                    else:
                        self.dbToModelSeed[source][externalId] = [modelSeedId]

                else:
                    dic = {}
                    dic[externalId] = [modelSeedId]
                    self.dbToModelSeed[source] = dic
        info.close()
        for line in extra_info:

            if url:
                line = line.decode("utf-8")

            l = line.split("\t")

            modelSeedId = l[0]
            externalId = l[1]
            source = l[2].strip()

            if re.search("\AModelSEED", line) == None:

                if modelSeedId in self.modelSeedToDb.keys():
                    if source in self.modelSeedToDb.get(modelSeedId) and \
                        externalId not in self.modelSeedToDb[modelSeedId][source]:
                        self.modelSeedToDb[modelSeedId][source].append(externalId)
                    else:
                        self.modelSeedToDb[modelSeedId][source] = [externalId]

                if source in self.dbToModelSeed.keys():
                    if externalId in self.dbToModelSeed.get(source).keys():
                        self.dbToModelSeed[source][externalId].append(modelSeedId)
                    else:
                        self.dbToModelSeed[source][externalId] = [modelSeedId]

                else:
                    dic = {}
                    dic[externalId] = [modelSeedId]
                    self.dbToModelSeed[source] = dic

        extra_info.close()