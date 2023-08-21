import re
import hashlib
import tomlkit as toml


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


def is_integer_num(n):
    if isinstance(n, int):
        return True
    if isinstance(n, float):
        return n.is_integer()
    return False


def calculate_sha256(filename):
    the_hash = hashlib.sha256()

    # establish the byte array  to but things into
    bytes_ = bytearray(128 * 1024)

    mem_view = memoryview(bytes_)

    with open(filename, "rb", buffering=0) as f:
        for n in iter(lambda: f.readinto(mem_view), 0):
            the_hash.update(mem_view[:n])

    return the_hash.hexdigest()


def get_toml_data(filename):
    """Function that can read a TOML file and create a table

    The main purpose of this function is to make it easier to
    load a TOML file in Python. Instead of having to write the
    `open` block of code, you can just pass through a filename
    and it will read the entire file.

    Parameters
    ----------
    filename : str
        The input file name for the TOML data to load.

    Returns
    -------
    tomlkit.Table
        A "wrapped" dict object that contains the key and
        item pairs found inside the TOML file.
    """
    with open(filename, "r") as in_file:
        toml_data = toml.parse(in_file.read())

    return toml_data

def get_toml_string(in_dict):
    """
    
    This function literally exists to not import TOML elsewhere
    """

    return toml.dumps(in_dict)


def write_inline_toml_dicts(in_dict, parent_name=""):
    """Takes an input dictionary and writes an inline toml dict

    Parameters
    ----------
    in_dict : dict
        The dictionary to write
    """

    # NOTE: the dictionary should essentially be a named list and should *not* have dicts

    out_str = "[" + parent_name + "]\n"

    for key, items in in_dict.items():
        out_str += key + " =  { "
        for ii, (subkey, subitems) in enumerate(items.items()):

            if isinstance(subitems, dict):
                raise NotImplementedError("Cannot have inline dict for this part")

            out_str += str(subkey) + " = "

            if isinstance(subitems, str):
                out_str += f'"{subitems}"'
            elif isinstance(subitems, list):
                out_str += "[" + ",".join([f"{x}" for x in subitems])
                out_str += "]"
            elif isinstance(subitems, bool):
                out_str += "true" if subitems else "false"
            else:
                out_str += f"{subitems}"

            out_str += " }\n" if ii == len(items) - 1 else ", "

    return out_str

