import string
from numpy import *

def plotTikz(data, extraData, filename):
  """Plots tabular data in data as stacked bars and a table. 
     data must be a python 2D list, where columns form the 
     x-axis and the rows form the y-axis.
     
     The tikz code is written out to the file specified by
     filename."""

  # open output file
  '''
  clrs = ['blue', 'red', 'green', 'orange', 'yellow','lightsalmon']
  rowLabels = ['communication (s)', 'unzip (s)', 'rhs (s)', 'derivatives (s)','wavelets (s)']
  exLabels = ['Total time', 'Total dofs (zipped)']
  vals = ['vala', 'valb', 'valc', 'vald', 'vale', 'valf', 'valg', 'valh', 'vali', 'valj', 'valk']'''
  #strong scaling labels
  clrs = ['blue', 'red', 'green', 'orange', 'yellow']
  rowLabels = ['communication (s)', 'unzip (s)', 'rhs (s)', 'derivatives (s)',]
  exLabels = ['Total time (s)']
  vals = ['vala', 'valb', 'valc', 'vald', 'vale', 'valf', 'valg', 'valh', 'vali', 'valj', 'valk']

  out = open(filename, "w")

  str = '\\documentclass{standalone}\n\\usepackage{tikz}\n\n\\begin{document}\\pagestyle{empty}\n\n'
  out.write(str);

  # begin a tikz figure
  str = '\\begin{tikzpicture}\n\n'
  out.write(str)

  # draw the plot area, axes and grid.
  #str = '\\draw[very thin,color=gray] (-0.1,-0.1) grid (13,4.5);\n' + '\\draw[->] (-0.2,0) -- (13.2,0) node[right] {$p$};\n' + '\\draw[->] (0,-0.2) -- (0,4.6) node[above] {$time (\mu s)$};\n\n'
  str = '\\draw[very thin,color=gray] (-0.1,-0.1) grid (5,4.5);\n' + '\\draw[->] (-0.2,0) -- (5.2,0) node[right] {$p$};\n' + '\\draw[->] (0,-0.2) -- (0,4.6) node[above] {$time (s)$};\n\n'
  out.write(str)


  # Label the proc counts ...
  str = '\\foreach \\pos/\\label in {0.5/4K,1.5/8K,2.5/16K,3.5/32K,4.5/65K}\n' + '\\draw (\\pos,0) -- (\\pos,-0.1) (\\pos cm,-2.5ex) node [anchor=base,fill=white,inner sep=1pt]  {\\scriptsize \\label};\n\n'
  # str = '\\foreach \\pos/\\label in {0.5/32,1.5/64,2.5/128,3.5/256,4.5/512,5.5/1024,6.5/2K,7.5/4K,8.5/8K,9.5/16K,10.5/32K,11.5/65K,12.5/131K}\n' + '\\draw (\\pos,0) -- (\\pos,-0.1) (\\pos cm,-2.5ex) node [anchor=base,fill=white,inner sep=1pt]  {\\scriptsize \\label};\n\n'
  out.write(str)

  # the table ...
  sz = shape(data)
  esz = len(extraData);
  str = '\\draw[very thin,color=gray, xstep=1, ystep=0.5] (0, ' + repr(-0.5*(sz[0]+esz+1)) + ') grid (5,-0.5);\n';
  out.write(str)

  ## compute the max range/scaling of the plot ...
  # mxHt = sum(data[:,-1]);
  mxHt = sum(data[:,0]);
  print(mxHt);
  ySc = 4.4/mxHt;

  # print mxHt, ySc

  # Labels
  for i in range(sz[0]):
    str = '\\draw[fill=' + clrs[i] + '!50, fill opacity=0.5] (-3,'+repr(-0.5*(i+2)) + ') rectangle +(3,0.5);\n'
    str = str + '\\draw (-1.5, ' + repr(-0.5*(i+1.5))  + ') node {\\scriptsize ' + rowLabels[i] + '};\n\n'
    out.write(str);
  for i in range(esz):
    str = '\\draw (-3,'+repr(-0.5*(sz[0]+i+2)) + ') rectangle +(3,0.5);\n'
    str = str + '\\draw (-1.5, ' + repr(-0.5*(sz[0]+i+1.5))  + ') node {\\scriptsize ' + exLabels[i] + '};\n\n'
    out.write(str);

  for i in range(5):
    val = i/ySc
    str = '\\draw (-2.5ex, ' + repr(i) + ' cm) node [anchor=base,fill=white,inner sep=1pt] {\\scriptsize ' + '%d'%val + '};\n\n'
    out.write(str);

  
  ## The actual data, into the plot and the table.
  
  # define the working variables ...
  str = '\\newdimen\\mypos\n' + '\\newdimen\\myoff\n\n'
  out.write(str)

  str = '\\foreach \\pos/'
  for i in range(sz[0]):
    str = str + '\\' + vals[i] + '/'

  str = str[:-1] + ' in { ' 
  
  for i in range(sz[1]):
    str = str + repr(i) + '/'
    for j in range(sz[0]):
      str = str + '%3.2f'%data[j][i] + '/'
    str = str[:-1] + ', '
  str = str[:-2] + '} { \n'
  out.write(str)


  str = '\\mypos=\\pos cm\n'
  str = str + '\\advance \\mypos by 0.5 cm\n'
  out.write(str)

  # Enter the text in the table ...
  for i in range(sz[0]):
    str = '\\draw (\\mypos, ' + repr(-0.5*(i+1.5))  + ' cm) node {\\scriptsize \\' + vals[i] + '};\n'
    out.write(str);



  str = '\n\\myoff=0 cm\n\n'
  str = str + '\\advance \\mypos by -0.25 cm\n'
  out.write(str)

  for i in range(sz[0]):
    str = '\\draw[fill=' + clrs[sz[0]-i-1] + '!50, fill opacity=0.5, yscale='+ repr(ySc) + '] (\\mypos,\\myoff) rectangle +(0.5,\\' + vals[sz[0]-i-1] + ');\n'
    str = str + '\\advance \\myoff by \\' +  vals[sz[0]-i-1] + ' cm\n'
    out.write(str)
 
  str = '\n}\n\n'
  out.write(str)

  for j in range(sz[1]):
    for i in range(esz):
      
      str = '\\draw (' + repr(j+0.5) + ', ' + repr(-0.5*(sz[0]+i+1.5))  + ') node {\\scriptsize ' + extraData[i][j] + '};\n'
      out.write(str);


# close the tikz figure
  str = '\\end{tikzpicture}\n'
  out.write(str)
  
  str = '\\end{document}\n'
  out.write(str)

  # close the file
  out.close()
