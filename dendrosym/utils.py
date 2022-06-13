import re


def unique(in_list):
    unique_list = []

    for x in in_list:
        if x not in unique_list:
            unique_list.append(x)

    unique_list.sort()

    return unique_list


def find_derivtype_in_text(in_text, search_prefix="grad_"):

    re_search = r"\b" + search_prefix + r"\w+"

    found_instances = re.findall(re_search, in_text)

    return unique(found_instances)
