import re

x=0;
fw=open('cuda_bssneqns.h','w+');

varInArray=['alpha','chi','K','gt0','gt1','gt2','gt3','gt4','gt5','beta0','beta1',
'beta2','At0','At1','At2','At3','At4','At5','Gt0','Gt1','Gt2','B0','B1','B2'];
varOutArray=['a_rhs','chi_rhs','K_rhs','gt_rhs00','gt_rhs01','gt_rhs02','gt_rhs11','gt_rhs12','gt_rhs22','b_rhs0','b_rhs1','b_rhs2'
,'At_rhs00','At_rhs01','At_rhs02','At_rhs11','At_rhs12','At_rhs22','Gt_rhs0','Gt_rhs1','Gt_rhs2','B_rhs0','B_rhs1','B_rhs2'];

def replaceAndWrite(line):
	for i in range(len(varInArray)):
		matches=list(re.finditer('[/\+\*\-( ]' +varInArray[i]+'\[pp\]',line));
		x=0;
		while(len(matches)>0):
			line=line[0:matches[0].start(0)+1]+'dev_var_in['+varInArray[i]+'Int+pp]'+line[matches[0].end(0):];
			matches=list(re.finditer('[/\+\*\-( ]' +varInArray[i]+'\[pp\]',line)) ;
		matches=list(re.finditer(varOutArray[i]+'\[pp\]',line)) ;
		while(len(list(matches))>0):
			line=line[0:matches[0].start(0)]+'dev_var_out['+varInArray[i]+'Int+pp]'+line[matches[0].end(0):];
			matches=list(re.finditer(varOutArray[i]+'\[pp\]',line)) ;
	fw.write(line);

with open('bssneqs.cpp') as f:
	lines = f.readlines()
	for line in lines:
		replaceAndWrite(line);
		x+=1
		print x
