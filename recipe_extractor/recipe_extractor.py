from materials_entity_recognition import MatRecognition
from material_parser import MaterialParser
from synthesis_project_ceder.nlp.preprocessing import TextPreprocessor

class RecipeExtractor:
    def __init__(self):
        self.__tp = TextPreprocessor('')
        self.__mer = MatRecognition()
        self.__mp = MaterialParser()

    def get_materials(self, doi, abstract, syn_paragraph, split_mixture=False):

        structures_list = []

        #Text preprocessing: cleaning up paragraphs
        #MER
        #get dopants
        #get phase
        #build abbreviations dictionaty
        #substitute abbreviations
        #clean up materials names
        #chemical names substitutions from pre-built dictionary

        #parse chemical formulas of materials
        if split_mixture:
            #split mitures into compounds

        #substitute variables
        #append materials strictures into structures list

        return structures_list