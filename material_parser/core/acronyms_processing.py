import regex as re
import material_parser.core.regex_parser as rp


def __is_abbreviation_like(word):
    if all(c.isupper() for c in re.sub(r"[0-9x\-\(\)\.]", "", word)) \
            and len(re.findall(rp.re_capitals_no_O, word)) > 1:
        return True

    return False


def build_acronyms_dict(materials_list, paragraph):
    """
    constructing dictionary of acronyms appeared in material list
    :param paragraph: <list> of <str> list of sentences to look for abbreviations names
    :param materials_list: <list> of <str> list of materials entities
    :return: <dict> abbreviation: corresponding string
    """

    acronyms_dict = {t: "" for t in materials_list if __is_abbreviation_like(t.replace(" ", "")) and t != ""}
    full_names = [m for m in set(materials_list) if m not in acronyms_dict]

    """
    first find acronyms among materials in the list
    """
    for acronym in acronyms_dict.keys():
        for m_name in full_names:
            if sorted(re.findall(rp.re_capitals_no_O, acronym)) == sorted(re.findall(rp.re_capitals_no_O, m_name)):
                acronyms_dict[acronym] = m_name

    """
    for the rest of acronyms search the text
    """
    # TODO: need more robust implementation
    for acronym, m_name in acronyms_dict.items():
        sentences = " ".join([s for s in paragraph if acronym in s]).split(acronym)
        i = 0
        while not acronyms_dict[acronym] and i < len(sentences):
            sent = sentences[i]
            for tok in sent.split(" "):
                if sorted(re.findall(rp.re_capitals_no_O, tok)) == sorted(re.findall(rp.re_capitals_no_O, acronym)):
                    acronyms_dict[acronym] = tok
            i += 1

    """
    checking the case when acronyms consists of several words: e.g. KNN-BT
    """
    for acronym in acronyms_dict.keys():
        parts = re.split("-", acronym)
        if all(p in acronyms_dict for p in parts) and acronyms_dict[acronym] == "" and len(parts) > 1:
            name = "".join("(" + acronyms_dict[p] + ")" + "-" for p in parts).rstrip("-")
            acronyms_dict[acronym] = name

    return {acronym: name for acronym, name in acronyms_dict.items() if name}