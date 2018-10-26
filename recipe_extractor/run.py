# coding=utf-8
import recipe_extractor

from pymongo import MongoClient
from bson.objectid import ObjectId
import random
import json

from pprint import pprint

#solid_state_papers = json.loads(open('/home/olga/Desktop/SynthesisProject/MaterialsDB/v2/solid_state_papers.json').read())
#solid_state_data = [dict(doi=p['doi'], abstract='', synthesis=''.join([s+' ' for s in p['text']])) for p in solid_state_papers

rex = recipe_extractor.RecipeExtractor(verbose=False, pubchem_lookup=True)

solid_state_data = json.loads(open('/home/olga/Desktop/SynthesisProject/Materials&Operations/v1/failed_parser.json').read())

doi = '10.1149/1.1350658'


updated = []
for data in solid_state_data:
    if data['doi'] == doi:
        text = ''.join(sent+' ' for sent in data['text'])
        extracted_materials, mer, failed_materials = rex.get_materials(data['doi'], '', text)
        reactions = rex.get_reactions(extracted_materials['targets'], extracted_materials['precursors'], split_mixture=True)

        print ('==='*30)
        pprint(extracted_materials)
        pprint (mer)
        pprint (failed_materials)
        print ('\nReactions:')
        pprint (reactions)



        #id_ = random.choice(seq=range(len(solid_state_data[0:5000])))
#id_ = 1194
#id_ = 1446
#id_ = 1401
#id_ = 3717
#id_ = 777
# 10.1007/s10854-012-0826-2
# 10.1016/j.jcrysgro.2009.11.001
#test_doi = ['10.1016/j.cplett.2012.03.076', '10.1016/j.optmat.2015.01.027', '10.1016/j.cplett.2016.05.035', '10.1111/j.1551-2916.2008.02672.x', '10.1111/j.1551-2916.2007.01825.x', '10.1111/jace.13158']
#test_doi_2 = ['10.1016/j.cplett.2012.03.076', '10.1016/j.cplett.2016.05.035', '10.1016/j.optmat.2015.01.027']
#test_doi_3 = ['10.1016/j.optmat.2015.01.027']
#test_doi_3 = ['10.1007/s10832-016-0030-5']
# test_doi_3 = ['10.1149/1.1342168']
# #id_ = 0
# #while solid_state_data[id_]['doi'] != '10.1016/j.cplett.2012.03.076':
# #     id_ = id_ + 1
# #print (id_)
# for doi in test_doi_3:
#     id_ = 0
#     while solid_state_data[id_]['doi'] != doi:
#          id_ = id_ + 1
#     doi = solid_state_data[id_]['doi']
#     abstract = solid_state_data[id_]['abstract']
#     #synthesis_par = solid_state_data[id_]['synthesis']
#     #synthesis_par = 'Electrochemical characterization was performed in a Swagelok cell or standard coin cell configuration. Anodes used were Li metal (Aldrich), Li4Ti5O12, or graphite [Osaka gas; mesocarbon microbeads (MCMB) 25-28]. Li4Ti5O12 was fabricated by solid-state synthesis. Stoichiometric amounts of Li2CO3 (Aldrich) and TiO2 (Aldrich) were mixed thoroughly in an agate mortar and pestle and fired at 800°C under an air atmosphere. Liquid electrolyte-type electrodes were fabricated by mixing 85% spinel, 10% polyvinylidene fluoride binder (Aldrich), and 5% conductive carbon black in 1-methyl- 2-pyrrolidinone solvent for 5 h. The slurry was then cast on aluminum foil and dried at 120°C for 1 h in <1% relative humidity. Disks were then punched out of the foil. All negative electrodes besides lithium were made by the PLiON plastic process. 65% active material was mixed with 6.5% carbon black (MMM super P), 10% poly(vinylidene fluoride-co-hexafluoropropylene) binder (Elf Atochem, Kynar 2801), and 18% dibutyl phthalate plasticizer in acetone. The mixture was cast and dried at 22°C for 0.5 h. Afterward, disks were punched from the freestanding tape. The disks were then placed in ether to extract the plasticizer. All liquid cells used Whatman GF/D borosilicate glass fiber separators. Electrolyte was 1 M LiPF6 dissolved in ethylene carbonate:dimethyl carbonate in a 2:1 volume ratio (Mitsubishi chemical). Plastic lithium ion cells were fabricated using the PLiON process. All electrochemical cell fabrication was performed in a helium-filled glove box with a −80°C dew-point atmosphere. '
#     synthesis_par = 'Similary, LiAlxMn2-xO4-zFz samples were also fabricated using Al2O3 and LiF as precursors .  As a percentage of F is exchanged for oxygen during the fabrication process, the fluorine compositions reported here are, unless otherwise noted, considered to be nominal .'
#
#     res, fails = extractor.get_materials(doi, abstract, synthesis_par, split_mixture=True)
#     # print (doi)
#     # print (abstract)
#     # print (synthesis_par)
#     print('Final result:')
#     pprint(res)
#     #if not fails['obtained_reactions']:
#     #    print ('No reactions found!')
#     #for r in res:
#     #    print (r['correct_formula'])
#     #    print (r['precursors'])
#     print ('----------------------------')


#print (doi)
#print (abstract)
#print (synthesis_par)
#pprint(res)

print('Done')