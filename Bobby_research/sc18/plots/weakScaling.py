#author: Milinda Fernando
#school of computing, university of utah.
import pandas as pd
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


#ifile='titan/raw/bssn_r1_im.prof'
#ofile_unique='titan/r1/bssn_r1_im1.prof'
#ofile_ws='titan/r1/bssn_r1_im_ws1.prof'

ifile='../dat/titan/raw/bssn_r10_im_g1000.prof'
ofile_unique='../dat/titan/r10/bssn_r10_im_g1000.prof'
ofile_ws='../dat/titan/r10/bssn_r10_im_ws_g1000.prof'

#ifile='titan/raw/bssn_r10_im_g100.prof'
#ofile_unique='titan/r10/bssn_r10_im_g100.prof'
#ofile_ws='titan/r10/bssn_r10_im_ws_g100.prof'

grainsz=60000 #### needed grain size to extract. 
pm=5000     ##### range it does search
grainsz_min=grainsz-pm
grainsz_max=grainsz+pm
df=pd.read_csv(ifile, sep='\t', lineterminator='\n',encoding="utf-8-sig")

df=df.loc[:, ~df.columns.str.contains('^Unnamed')]
df.columns = df.columns.str.replace('\s+', '')
#ds=pd.DataFrame(ds)
#weak_ds=ds.loc[((int(ds['dof_zip'])/float(ds['act_npes1']))<=grainsz_max) and  ((int(ds['dof_zip'])/float(ds['act_npes1'])) >=grainsz_min)]
#print(df)
df = df.drop((df[df.rkstep_mean == 0].index))
df=df.drop_duplicates('step',keep='last',inplace=False)
#print(df)
df.to_csv(ofile_unique, sep="\t", quoting=csv.QUOTE_NONE,index=False)
#print(df.columns)
#weak_df=df.loc[(((df.numOcts)/(df.act_npes))>=grainsz_min) & (((df.numOcts)/(df.act_npes)) <= grainsz_max)]
weak_df=df.loc[(((df.dof_zip)/(df.act_npes))>=grainsz_min) & (((df.dof_zip)/(df.act_npes)) <= grainsz_max)]
weak_df=weak_df.drop_duplicates('act_npes',keep='first',inplace=False)
weak_df.to_csv(ofile_ws, sep="\t", quoting=csv.QUOTE_NONE,index=False)

print(weak_df)

#x=weak_df['act_npes']
#y=weak_df['rkstep_mean']/10.0
#x_ticks=list(map(str, weak_df['act_npes'].tolist()))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
#plt.plot(x, y)
#plt.axes().xaxis.set_ticks(x)
#plt.axes().xaxis.set_ticklabels(x_ticks)
#plt.show()
#print(weak_df)
#print(df.columns)
#print(df['act_npes'])
