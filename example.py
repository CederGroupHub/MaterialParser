# coding=utf-8
import json
from pprint import pprint
import random as rnd

from material_parser.core.material_parser import MaterialParserBuilder
from material_parser.core.preprocessing_tools.additives_processing import AdditivesProcessing
from material_parser.core.preprocessing_tools.chemical_name_processing import ChemicalNameProcessing
from material_parser.core.preprocessing_tools.mixture_processing import MixtureProcessing
from material_parser.core.preprocessing_tools.phase_processing import PhaseProcessing
from material_parser.core.postprocessing_tools.substitute_additives import SubstituteAdditives
from material_parser.core.string_cleanup import cleanup_name

mp = MaterialParserBuilder() \
    .addPreprocessing(AdditivesProcessing()) \
    .addPreprocessing(ChemicalNameProcessing()) \
    .addPreprocessing(PhaseProcessing()) \
    .addPreprocessing(MixtureProcessing())\
    .addPostprocessing(SubstituteAdditives())\
    .build()

dataset = json.loads(open("material_parser/test_data.json").read())
data = dataset[rnd.randint(0, len(dataset)-1)]

material = data["material"]
output = mp.parse(cleanup_name(material))
pprint(output.to_dict())

print ('Done!')




