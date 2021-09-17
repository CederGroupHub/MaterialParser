# MaterialParser

This is an old repo for the text2chem parser. The new one is available at https://github.com/CederGroupHub/text2chem.

Material Parser is a package that provides various functionality to parse material entities and convert them into a
chemical structure.

The Material Parser addresses the problem of unification of materials entities widely used in scientific literature.

This tool is build to facilitate and promote text mining research in material science.

**Material Parser** functionality includes:

 * converting standard chemical terms into a chemical formula
 * parsing chemical formula into chemical composition,
 * constructing dictionary of materials abbreviations from a text snippets,
 * finding the values of stoichiometric and elements variables from text snippets,
 * splitting mixtures/composites/alloys/solid solutions into compounds

Note: Material Parser is not intended to be used for finding material entities in the text or to perform any
chemical named entities recognition (NER) in a text.
It rather processes and standardizes already extracted chemical (material) entities and extracts relevant chemical
information from them.

The parser was tested to work well for inorganic material terms, its performance on organic terms was not
evaluated thoroughly.

### Pipeline overview

Material Parser runs as a pipeline of default, pre- and post-processing modules.
Default module of the pipeline process chemical formula and converts into a chemical structure.
Preprocessing modules deal with an input material string,
separate relevant information from an input string and fill it into an output chemical structure.
Postprocessing modules augment chemical structure with any other information that can be found in a surrounding text.

The current version of Material Parser package includes the following preprocessing modules:
 * ``PhaseProcessing`` - gets information about material phase if it is presented in the string,
e.g. for "P2-Na2/3(CoxNi1/3-xMn2/3)O2", "P2" will be separated and stored in the structure
while the rest of the string will be sent down the pipeline.
 * ``AdditivesProcessing`` - gets information about any additives, such as dopants, stabilizers or mixed compounds,
e.g. i) for "(Na0.5K0.5)NbO3 + 1.5 mol% CuF2", "CuF2" will be separated as an additive
and "(Na0.5K0.5)NbO3" will be sent down the pipeline; ii) for "Ca2BO3Cl:Sm3+, Eu3+", "Sm3+" and "Eu3+" will be separated
and "Ca2BO3Cl" will be sent down the pipeline.
 * ``ChemicalNameProcessing`` - attempts to convert sequence of chemical terms into a formula,
e.g. "zinc (II) acetate dihydrate" will be converted into "Zn(CH3COO)2·2H2O" that will be sent down the pipeline
and the chemical terms will be stored in output structure as material name.
This tool can also separate chemical formula from the rest of the terms is they are combined by tokenizer,
e.g.  "ammonium molybdate (NH4)6Mo7O24⋅4H2O" will be split into "ammonium molybdate" - material name
and "(NH4)6Mo7O24⋅4H2O" - chemical formula that will be sent down the pipeline.
 * ``MixtureProcessing`` - split mixture/solid solution/alloy/composite/hydrate into constituting compounds,
e.g. "(1-x-y)BaTiO3-xBaBiO3-y(Bi0.5Na0.5)TiO3" will be split into "BaTiO3", "BaBiO3", and "(Bi0.5Na0.5)TiO3"
with the corresponding amounts "(1-x-y)", "x" and "y", respectively.
 * ``PubchemPreprocessing`` - looks up chemical terms in PubChem database, if chemical formula is not found.
This step usually slows down the overall pipeline perfomance.

The postprocessing module includes:
 * ``Substitute_Additives`` - appends dopand elements to the formula to have integer total stoichiometry
or add mixture compound to the composition,
e.g. i) "Zn1.92-2xYxLixSiO4:0.08Mn2+" becomes "Mn0.08Zn1.92-2xYxLixSiO4" and
ii) "(Na0.5K0.5)NbO3 + 1.5 mol% CuF2" becomes "(1-x)(Na0.5K0.5)NbO3-(x)CuF2".


### Installation:

The installation can be performed directly via `pip`: 
```
pip install git+https://github.com/CederGroupHub/MaterialParser.git
```

or, in two steps: 
```
git clone https://github.com/CederGroupHub/MaterialParser.git
cd MaterialParser
pip install -r requirements.txt -e .
```

### Running default setup:

By default, Material Parser gets chemical composition from a proper material formula
(i.e. no chemical terms, dopants or mixtures):

```
from material_parser.core.material_parser import MaterialParserBuilder
mp = mp = MaterialParserBuilder().build()
output = mp.parse("Li5+xLa3Ta2-xGexO12")
print(output.to_dict())

```
Output:
```
{'additives': [],
 'amounts_x': {'x': []},
 'composition': [{'amount': '1',
                  'elements': OrderedDict([('Li', 'x+5'),
                                           ('La', '3'),
                                           ('Ta', '2-x'),
                                           ('Ge', 'x'),
                                           ('O', '12')]),
                  'formula': 'Li5+xLa3Ta2-xGexO12',
                  'species': OrderedDict([('Li', 'x+5'),
                                          ('La', '3'),
                                          ('Ta', '2-x'),
                                          ('Ge', 'x'),
                                          ('O', '12')])}],
 'elements_x': {},
 'material_formula': 'Li5+xLa3Ta2-xGexO12',
 'material_name': '',
 'material_string': 'Li5+xLa3Ta2-xGexO12',
 'oxygen_deficiency': '',
 'phase': ''}
```
Other attributes can be filled using additional pipeline blockes if required.

### Adding modules to the pipeline:
```
from material_parser.core.material_parser import MaterialParserBuilder
from material_parser.core.preprocessing_tools.additives_processing import AdditivesProcessing
from material_parser.core.preprocessing_tools.chemical_name_processing import ChemicalNameProcessing
from material_parser.core.preprocessing_tools.mixture_processing import MixtureProcessing
from material_parser.core.preprocessing_tools.phase_processing import PhaseProcessing
from material_parser.core.postprocessing_tools.substitute_additives import SubstituteAdditives

mp = MaterialParserBuilder() \
    .addPreprocessing(AdditivesProcessing()) \
    .addPreprocessing(ChemicalNameProcessing()) \
    .addPreprocessing(PhaseProcessing()) \
    .addPreprocessing(MixtureProcessing())\
    .addPostprocessing(SubstituteAdditives())\
    .build()

mp.parse("Nasicon (Na1+x+yZr2-yYySixP3-xO12)")

```
Output:
```
{'additives': [],
 'amounts_x': {'x': [], 'y': []},
 'composition': [{'amount': '1',
                  'elements': OrderedDict([('Na', 'x+y+1'),
                                           ('Zr', '2-y'),
                                           ('Y', 'y'),
                                           ('Si', 'x'),
                                           ('P', '3-x'),
                                           ('O', '12')]),
                  'formula': 'Na1+x+yZr2-yYySixP3-xO12',
                  'species': OrderedDict([('Na', 'x+y+1'),
                                          ('Zr', '2-y'),
                                          ('Y', 'y'),
                                          ('Si', 'x'),
                                          ('P', '3-x'),
                                          ('O', '12')])}],
 'elements_x': {},
 'material_formula': 'Na1+x+yZr2-yYySixP3-xO12',
 'material_name': 'Nasicon',
 'material_string': 'Nasicon (Na1+x+yZr2-yYySixP3-xO12)',
 'oxygen_deficiency': '',
 'phase': ''}

```
Note: the order of the modules may affect the resulted output.


### Running customized pipeline:

Material Parser provides capabilities for creating customized pre- and post-processing modules.
This are defined by the corresponding interface:
``core/preprocessing_tools/preprocessing_abc.py`` and ``core/postprocessing_tools/postprocessing_abc.py``.
Add the class implementing the interface into a corresponding directory and import as a regular module.

### Output

mp.parse() output the ChemicalStructure with the following attributes:

 * material_string: ``string``
 * material_name: ``string``
 * material_formula: ``string``
 * additives: ``list`` or ``string``
 * phase: ``string``
 * oxygen_deficiency: ``char``
 * amounts_x: ``Variables``
 * elements_x: ``Variables``
 * composition: ``list`` of ``Compound``


### Citing

If you use Material Parser in your work, please cite the following paper:

 * Kononova et. al "Text-mined dataset of inorganic materials synthesis recipes", Scientific Data 6 (1), 1-11 (2019)
 [10.1038/s41597-019-0224-1](https://www.nature.com/articles/s41597-019-0224-1)
