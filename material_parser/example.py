# coding=utf-8
from material_parser import MaterialParser
import pprint
from chemdataextractor.doc import Paragraph


def find_variables(var, text):
    sents = Paragraph(text).sentences
    i = 0
    sent = sents[i]
    values = []
    mode = ''
    try:
        while len(values) == 0 and i < len(sents):
            sent = sents[i]
            values, mode = mp.get_stoichiometric_values(var, sent.text)
            i = i + 1
    except:
        print('Failed to extract ' + var + ' from text:')

    return values


mp = MaterialParser()

list_of_targets = ['(1-x)Na0.5Sm0.5TiO3-xSm(Mg0.5Ti0.5)O3', 'NST-SMTx', 'Na0.5Sm0.5TiO3',
                   'Na0.5Sm0.5[Ti1-y(Al0.5Ta0.5)y]O3', 'NST', 'NSTATx']
paragraphs = [
    'The effect of calcining and sintering conditions on the microwave dielectric properties of Na0.5Sm0.5TiO3 (NST) \
    solid solutions were investigated. The results showed that, after being calcined at 1100 °C for 2 h, the proposed \
    NST system sintered at 1350 °C for 4 h exhibited a good microwave dielectric properties. \
    Additionally, the phase assemblages, crystal structures, microstructures and microwave dielectric properties of \
    the Na0.5Sm0.5[Ti1-y(Al0.5Ta0.5)y]O3 (NSTATx, 0.1 ≤ y ≤ 0.3) and (1 - x)Na0.5Sm0.5TiO3-xSm(Mg0.5Ti0.5)O3 \
    (NST-SMTx, 0.1 ≤ x ≤ 0.4) ceramics were also investigated in this work. \
    The X-ray diffractometer results revealed that a tilted orthorhombic perovskite structure in space group Pnma \
    was refined in the NSTATx ceramics, while the Ti7O13 phase, Sm2Ti2O7 phase, and some unknown phase were detected \
    gradually in NST-SMTx ceramics with increasing x value.',
    'Na0.5Sm0.5TiO3 (NST), Na0.5Sm0.5[Ti1-x(Al0.5Ta0.5)x]O3 (NSTATx) and (1 - x)Na0.5Sm0.5TiO3-xSm(Mg0.5Ti0.5)O3 \
    (NST-SMTx) specimens were synthesized by a conventional mixed-oxide route using high purity (≥99.5 %) oxide powders \
    of Sm2O3, Al2O3, Ta2O5, MgO, Na2CO3, and TiO2. \
    The starting NST powders were weighed and mixed with ZrO2 balls in ethanol for 24 h, then dried, sieved through \
    100-mesh screen, and calcined at 900-1100 °C for 2 h in alumina crucible. \
    Additionally, the synthesized powders were weighed according to stoichiometric ratio of \
    NSTATx (y = 0.1, 0.2 and 0.3) and NST-SMTx (x = 0.1, 0.2, 0.3 and 0.4), and then ground in ethanol for 24 h with \
    ZrO2 balls. Mixtures were dried, sieved and calcined at 1100 °C for 2 h. \
    All the calcined powders were ground and mixed with 5 wt% PVA poly(vinyl alcohol) as binder. \
    The granulated powders were uni-axially pressed into pellets with 11.5 mm in diameter and 1.5–6 mm in height under \
    the pressure of 180 MPa. After being fired to remove an organic binder, these pellets were sintered at \
    1300-1520 °C for different times in air.']


print(list_of_targets)
print(paragraphs)

data = {}
data['targets'] = list_of_targets
data['paragraphs'] = paragraphs

# Find doped elements:
data['dopants'] = []
data['targets_upd'] = []
for material in list_of_targets:
    dopants, new_material = mp.get_dopants(material)

    data['dopants'].extend(d for d in dopants)
    data['targets_upd'].append(new_material)

# Resolve abbreviations:
data['abbreviations'] = mp.build_abbreviations_dict(data['targets_upd'], data['paragraphs'])
data['targets_upd_2'] = []
for material in data['targets_upd']:
    if material in data['abbreviations']:
        data['targets_upd_2'].append(data['abbreviations'][material])
    else:
        data['targets_upd_2'].append(material)

print('Initial targets:', data['targets'])
print('After removing dopants:', data['targets_upd'])
print('After assining abbreviations:', data['targets_upd_2'])
print('Found dopants', data['dopants'])
print('Abbreviations dictionary:')
pprint.pprint(data['abbreviations'])

# Resolving chemical structure
print('Chemical structure:')
data['chemical_structure'] = {}
for target in data['targets_upd_2']:
    print(target)
    try:
        t_struct = mp.get_chemical_structure(target)
        data['chemical_structure'][target] = t_struct
    except:
        print('Failed!')

# Resolving variables:
for name, struct in data['chemical_structure'].items():
    for var in struct['fraction_vars']:
        values = find_variables(var, paragraphs[1])

        if len(values) != 0:
            struct['fraction_vars'][var] = values
            print("In chemical formula: "+name)
            print("Found variable: "+var)
            print("Found values: "+str(values))

pprint.pprint(data['chemical_structure'])

print('-----------------------------------------------')






