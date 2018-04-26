def replace_unknown_characters(item_list):
    new_items = []
    for i in item_list:
        new_items.append(i.replace(",","COMMA").replace(" ","SPACE"))
    return new_items
