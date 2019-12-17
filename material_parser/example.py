# coding=utf-8
from  material_parser import MaterialParser
import json
from pprint import pprint


mp = MaterialParser(pubchem_lookup=False, verbose=False)

test_set = json.loads(open('test_data.json').read())
print(len(test_set))

for item in test_set:

    material = item['material']
    correct = item['parser_output']

    list_of_materials = mp.split_materials_list(material)
    list_of_materials = list_of_materials if list_of_materials != [] else [ material ]
    structure = []
    for m in list_of_materials:
        structure.append(mp.parse_material_string(m))
    if structure != correct:
        print("Mismatch for ", material)
        print('-'*40)
        for s, c in zip(structure, correct):
            for key in c.keys():
                if s[key] != c[key]:
                    print(key)
                    print("Correct:", c[key])
                    print("Found:", s[key])
        print('-' * 40)
        break
    else:
        print(material, "->", "ok!")

print ('Done!')




