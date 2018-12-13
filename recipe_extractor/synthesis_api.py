from synthesis_api_hub import api_method, Client
from synthesis_api_hub.apiegg import APIEgg
from synthesis_project_ceder.utils.environments import request_linear_algebra_single_threaded

from recipe_extractor.recipe_extractor import RecipeExtractor

__all__ = ['RecipeExtractorWorker']


class RecipeExtractorWorker(APIEgg):
    namespace = 'recipe_extractor'
    version = '2018121300'

    def __init__(self):
        request_linear_algebra_single_threaded()

        self.re = RecipeExtractor()

    @api_method
    def infer_topics(self, documents):
        if not isinstance(documents, list):
            raise ValueError('documents must be a list')

        for x in documents:
            if not isinstance(x, dict) \
                    or 'syn_paragraph' not in x \
                    or 'doi' not in x \
                    or 'abstract' not in x:
                raise ValueError('document must be a dict with keys '
                                 '"syn_paragraph", "doi" and "abstract"')

        results = []

        for doc in documents:
            structures_list, fails = self.re.get_materials(
                doc['doi'], doc['abstract'], doc['syn_paragraph']
            )
            result = {
                'structures': structures_list,
                'failed': fails,
            }

            results.append(result)

        return results, self.version
