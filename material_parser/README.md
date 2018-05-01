# Material Parser

Simple parser which allows for:

 * extraction of composition from chemical formulas, 
 * finding chemical names corresponding to abbreviations in text,
 * finding values of stoichiometric and elements parameters
 
### Installation:
```
git clone https://github.com/CederGroupHub/MaterialParser.git
cd MaterialParser
pip install .
```

### Initilization:
```
from material_parser import MaterialParser
mp = MaterialParser()
```

### Functions:

 * get_structure_by_formula(chemical_formula)

 * get_chemical_structure(chemical_name):

    *Input*: chemical name or chemical formula
    
    *Return*: structure
    ```
    {
    name: chemical_name, #original name of material
    formula: chemical_formula, #chemical formula used for parsing
    composition: {Element1: stoichiometry, Element2: stoichiometry, ...}, #composition without substitutions
    elements_vars: [Var1: [Value1, ...], ...], #list of elements variables and their values
    fraction_vars: [Var1: [Value1, ...], ...], #list of stoichiomentric variables and their values
    mixture: {Part_Name1:{composition: {Element1: stoichiometry, 
                                        Element2: stoichiometry, ...}
                          fraction: fraction_value}
              Part_Name2: ... } #parts of mixture with corresponding compositions and fractions
    }
    ```

 * get_abbreviations_dict(list_of_abbreviations, text)
 
    *Return*: dictionary
    ```
    {
     Abbreviation1: CorresponsingName1,
     Abbreviation2: CorresponsingName2,
     ...
    }
    ```
    
 * substitute_abbreviations(list_of_abbreviations, abbreviations_dict)
 
   *Return:* new list with substituted abbreviations
   
 * get_dopants(chemical_name)
   
   *Return:* list of dopants, chemical_name without dopant
   
 * get_stoichiometric_values(variable, sentence)
 
    *Return:* list of variable values found in the sentence
    
 * get_elements_values(variable, sentence)
 
    *Return:* list of elements for variable found in the sentence
 
 
