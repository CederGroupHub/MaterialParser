import regex as re
from collections import OrderedDict
from material_parser.core.postprocessing_tools.postprocessing_abc import PostprocessingABC
from material_parser.core.formula_processing import process_formula
from material_parser.core.utils import is_int, simplify
import material_parser.core.chemical_sets as cs


class SubstituteAdditives(PostprocessingABC):

    def process_data(self, chemical_structure, text_sentences=[]):
        additives = chemical_structure.additives

        if len(additives) == 0:
            return chemical_structure

        additive = additives[0]
        additive_data = process_formula(additive)

        if len(additive_data["elements"]) > 1:
            """
            mixture or solid solution, add additive as a compound to composition
            e.g. (Na0.5K0.5)NbO3 + 1.5 mol% CuF2 -> (1-x)(Na0.5K0.5)NbO3-xCuF2
            """
            chemical_structure = self.__substitute_compound(additive_data, chemical_structure)

        elif all(c.elements for c in chemical_structure.composition):
            """
            doped element with stoichiometry
            e.g. LiSr1-xPO4:Eux
            """
            formula, composition = self.__substitute_element(additive,
                                                             chemical_structure.material_formula,
                                                             chemical_structure.composition)

            if formula != chemical_structure.material_formula:
                chemical_structure.additives = []
                chemical_structure.material_formula = formula
                chemical_structure.composition = composition
            else:
                chemical_structure.additives = additives
        else:
            # print('-->', "Additive is something else")
            chemical_structure.additives = additives

        return chemical_structure

    @staticmethod
    def __substitute_compound(data, chem_structure):
        for compound in chem_structure.composition:
            compound.amount = compound.amount + "-x"
        chem_structure.add_compound({"formula": data.get("formula", ""),
                                     "amount": "x",
                                     "elements": data.get("elements", OrderedDict()),
                                     "species": data.get("species", OrderedDict())})
        chem_structure.additives = []
        return chem_structure


    @staticmethod
    def __substitute_element(additive, formula, composition):

        new_composition = []
        new_formula = formula

        if additive[-1] == "+":
            additive = additive.rstrip("+0987654321")

        """
        find any stoichiometric coefficient next to the additive and split the list of additives
        """
        r = r"^[x0-9\.]+|[x0-9\.]+$"
        coeff = re.findall(r, additive)
        element = [s for s in re.split(r, additive) if s != ""][0]

        """
        if not meaningfull result, skip substitution
        """
        if coeff == [] or element not in cs.list_of_elements:
            return new_formula, composition

        """
        consider every compound in composition
        substitute element if it makes total stoichiometry to sum up to integer
        """
        for compound in composition:
            expr = "".join(["(" + v + ")+" for e, v in compound.elements.items()]).rstrip("+")

            coeff = coeff[0] if not re.match("^[0]+[1-9]", coeff[0]) else "0." + coeff[0][1:]
            expr = expr + "+(" + coeff + ")"

            if is_int(simplify(expr)):
                new_name = element + coeff + compound.formula
                new_elements = OrderedDict([(el, v) for el, v in compound.elements.items()] + [(element, coeff)])
                new_elements.move_to_end(element, last=False)
                new_species = OrderedDict([(el, v) for el, v in compound.species.items()] + [(element, coeff)])
                new_species.move_to_end(element, last=False)

                new_composition.append({"formula": new_name,
                                        "amount": compound.amount,
                                        "elements": new_elements,
                                        "species": new_species})
                new_formula = new_formula.replace(compound.formula, new_name)
            else:
                new_composition.append({"formula": compound.formula,
                                        "amount": compound.amount,
                                        "elements": compound.elements,
                                        "species": compound.species})

        return new_formula, new_composition
