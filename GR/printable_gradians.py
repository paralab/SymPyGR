import re
def change_deriv_names(str):
    c_str=str
    derivs=['agrad','grad','kograd']
    for deriv in derivs:
        key=deriv+'\(\d, \w+\[pp\]\)'
        slist=re.findall(key,c_str)
        for s in slist:
            #print(s)
            w1=s.split('(')
            w2=w1[1].split(')')[0].split(',')
            #print(w1[0]+'_'+w2[0].strip()+'_'+w2[1].strip()+';')
            rep=w1[0]
            for v in w2:
                rep=rep+'_'+v.strip()
            #rep=rep+';'
            c_str=c_str.replace(s,rep)

    derivs2=['grad2']
    for deriv in derivs2:
        key=deriv+'\(\d, \d, \w+\[pp\]\)'
        slist=re.findall(key,c_str)
        for s in slist:
            #print(s)
            w1=s.split('(')
            w2=w1[1].split(')')[0].split(',')
            #print(w1[0]+'_'+w2[0].strip()+'_'+w2[1].strip()+';')
            rep=w1[0]
            for v in w2:
                rep=rep+'_'+v.strip()
            #rep=rep+';'
            c_str=c_str.replace(s,rep)
    return c_str


bssn = open("bssn.cpp", "r")
code = bssn.read()
code = change_deriv_names(code)
bssn.close()
bssn = open("bssn.cpp", "w")
code = bssn.write(code)

bssn = open("bssnOriginal.cpp", "r")
code = bssn.read()
code = change_deriv_names(code)
bssn.close()
bssn = open("bssnOriginal.cpp", "w")
code = bssn.write(code)


