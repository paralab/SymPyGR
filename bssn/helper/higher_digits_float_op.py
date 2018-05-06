from decimal import *

data = """ 
-0.015708398471300358267832564251875737681984901428222656250000 -- ADD
0.015708398471300358267832564251875737681984901428222656250000 -- SUB
"""

vars = dict()
for i in data.split('\n'):
    if i.strip()=="":
        continue
    else:
        i = i.split("--")
        var_name = i[1].strip()
        value = i[0].strip()
        vars[var_name] = Decimal(value)

# add = vars["DENDRO_28"] + vars["DENDRO_30"] + vars["DENDRO_32"]
# sub = vars["DENDRO_34"] + vars["DENDRO_36"]

print(vars["ADD"]+vars["SUB"])