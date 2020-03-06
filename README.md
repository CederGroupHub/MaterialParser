# MaterialParser

The class providing functionality to extract chemical data from a string of chemical terms/formulas/material names

This parser was created in order to address the problem of unification of materials entitites found in scieitific publications to facilitate text mining.

**Material Parser** functionality includes:

 * converting chemical terms into chemical formula
 * parsing chemical formula into composition,
 * constructing dictionary of materials abbreviations from a text snippets,
 * finding values of stoichiometric and elements variables froma a text snippets,
 * splitting mixtures/composites/alloys/solid solutions into compounds
 
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

### Initialization:
```
from material_parser import MaterialParser
mp = MaterialParser(verbose=False, pubchem_lookup=False, fails_log=False)
```

### Material parser

#### Initialization

 ```
 verbose: <bool> print supplemental information
 pubchem_lookup: <bool> look for unknown chemical names in PubChem (not implemented, slows down computations significantly)
 fails_log: <bool> outputs log of materials for which mp.parse_material failed (useful when long list of materials in processes)
 ```

#### Primary functionality

 * Main method to compile string of chemical terms/formulas into data structure
    ```
     mp.parse_material_string(material_string)
     ```

 * Method to convert chemical name into formula
    ```
    mp.string2formula(material_string)
    ```

 * Method to compile chemical formula into data structure containing composition
     ```
     mp.formula2composition(chemical_formula)
     ```

#### Auxiliary functions

 * Extracting snippets of the string recognized as doped elements, stabilizers, coatings, activators, etc.
    ```
    mp.separate_additives(material_string)
    ```

 * Spliting mixtures, alloys, composites, etc into list of constituting compounds with their fractions
    ```
    mp.split_formula_into_compounds(material_string)
    ```
 * Extracting species from material string
    ```
    mp.get_species(material_string)
    ```

#### Additional functionality

 * Constructing dictionary of acronyms based on provided list of materials strings and text
    ```
    mp.build_acronyms_dict(list_of_materials, text)
    ```

 * Looking for the values of elements variables in the text
    ```
    mp.get_elements_values(variable, text)
    ```

 * Looking for the values of the variables for stoichiometric amounts in the text
    ```
    mp.get_stoichiometric_values(variable, text)
    ```

 * Spliting in into list of chemical names material string in the format list of cations+anion
    ```
    mp.split_materials_list(material_string)
    ```

 * Substituting doped elements into original chemical formula to complete total stoiciometry to integer value
    ```
    mp.substitute_additives(list_of_additives, data_structure)
    ```

#### Citing

If you use Material Parser in your work, please cite the following paper:

 * Kononova et. al "Text-mined dataset of inorganic materials synthesis recipes", Scientific Data 6 (1), 1-11 (2019) [10.1038/s41597-019-0224-1](https://www.nature.com/articles/s41597-019-0224-1)
