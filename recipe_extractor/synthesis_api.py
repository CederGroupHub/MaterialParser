from synthesis_api_hub import api_method, Client
from synthesis_api_hub.apiegg import APIEgg

from recipe_extractor.recipe_extractor import RecipeExtractor

__all__ = ['RecipeExtractorWorker']


class RecipeExtractorWorker(APIEgg):
    namespace = 'recipe_extractor'
    version = '2018121300'

    def __init__(self):
        self.re = RecipeExtractor()

    @api_method
    def extract(self, doi, abstract, syn_paragraph):
        data_structure, mer_materials, fails = self.re.get_materials(
            doi, abstract, syn_paragraph
        )
        return data_structure, mer_materials, fails
