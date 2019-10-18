import numpy as np
from matplotlib import pyplot as pp
import os
import sys

##to run as follow: python pythplotEXT.py dataname param1 param2 param3 param4

cutoff=0.5

#extract name of data
namedir=sys.argv[1]
namedata=namedir+'.txt'

testdir=not os.path.isdir(namedir)
testfile=os.path.isfile(namedata)
if testdir:
	if testfile:
		os.system('python dirmaker.py')
		os.system('echo ''Directory did not exist and has been created.'' ')
	else:
		os.system('echo ''No data available with this name.'' ')
		exit()

#intepretation of input
var=sys.argv[2:7] #var[i] is : (ordinate axis), :: (abscissa axis) or [number]
abuff=sys.argv[2:7] #copy to avoid overwriting var

loc1=abuff.index(':')
#abuff[loc1:loc1+1]=[] #removed if 2 different symbols
loc2=abuff.index('::')#+1
llist=[0,1,2,3,4]
[pos1,pos2,pos3]=[x for x in llist if ((x!=loc1) and (x!=loc2))]
var1=int(var[pos1])
var2=int(var[pos2])
var3=int(var[pos3])

#collect size of data
#testtxt=not os.path.isfile(namedir+'/param_JJ_2.txt')
#if testtxt: os.system('cp '+namedir+'/param_JJ_2 '+namedir+'/param_JJ_2.txt')
ff=open(namedir+'/param_JJ_2.txt',"r")
fflines=ff.readlines()
ff.close()

nstepphi=fflines[9].split('\t')[1] #reminder: in python, remove 1 to the index of line...
nstepB=fflines[10].split('\t')[2]
nstepmuSC=fflines[11].split('\t')[1]
nstepmu=fflines[12].split('\t')[1]
nstepSCm=fflines[13].split('\t')[1]

nstepphi=int(nstepphi[1:-1])
nstepB=int(nstepB[1:-1])
nstepmuSC=int(nstepmuSC[1:-1])
nstepmu=int(nstepmu[1:-1])
nstepSCm=int(nstepSCm[1:-1])

#fill the array in the order of the data sheet #ORDER (of writing, better if like data file)
nstep=[nstepmu,nstepSCm,nstepmuSC,nstepB,nstepphi]

#controls the validity of input
if ((var1<1) or (var2<1) or (var3<1)):
	print("Indices start at 1.")
	exit()
if ((var1>nstep[pos1]) or (var2>nstep[pos2]) or (var3>nstep[pos3])):
	print("Maximal values for indices are:",nstep)
	exit()

##extract relevant data for the specific diagram as stated in the input
#collect the indices
lengths=[nstep[1]*nstep[2]*nstep[3]*nstep[4],nstep[2]*nstep[3]*nstep[4],nstep[3]*nstep[4],nstep[4],1]
start=1+(var1-1)*lengths[pos1]+(var2-1)*lengths[pos2]+(var3-1)*lengths[pos3] #position in terms of lines
indexsublist=[-1+start+jj*lengths[loc2]+kk*lengths[loc1] for kk in range(nstep[loc1]) for jj in range(nstep[loc2]) ]#if jj>10] #-1 for indexing & for-loops as if written with indents, ie the first changes more slowly###############################################################

#collect all the values of gap
f=open(namedir+'/'+namedata,"r")
flines=f.readlines()
fresult=[]
for x in flines:
    fresult.append(x.split('\t')[5])
f.close()
fresult[0:1]=[] #removes the title
nfrresult=np.array([float(fresult[i]) for i in indexsublist])

#focus on small gaps
nfresult=np.where(nfrresult<cutoff,nfrresult,cutoff)

#values for diagram title
varval1=flines[indexsublist[0]+1].split('\t')[pos1]
varval2=flines[indexsublist[0]+1].split('\t')[pos2]
varval3=flines[indexsublist[0]+1].split('\t')[pos3]
#values for diagram axis
minx=round(float(flines[indexsublist[0]+1].split('\t')[loc2]),2)
miny=round(float(flines[indexsublist[0]+1].split('\t')[loc1]),2)
maxx=round(float(flines[indexsublist[-1]+1].split('\t')[loc2]),2)
maxy=round(float(flines[indexsublist[-1]+1].split('\t')[loc1]),2)
if minx==maxx:
	maxx=maxx+0.1
	minx=minx-0.1
else:
    interv=(maxx-minx)/(nstep[loc2]-1)
    maxx=maxx+interv/2.
    minx=minx-interv/2.
if miny==maxy:
	maxy=maxy+0.1
	miny=miny-0.1
else:
    interv=(maxy-miny)/(nstep[loc1]-1)
    maxy=maxy+interv/2.
    miny=miny-interv/2.

#reshape to plot
nnfresult=nfresult.reshape(nstep[loc1],nstep[loc2])
nnnfresult=np.ndarray.tolist(nnfresult)
nnnfresult.reverse() #to flip the plot (phi increasing from bottom to top)

#name of file
varnames=['mu','SCm','muSC','B','phi'] #ORDER
interp='none' #interpolation = none or bilinear

sm=''
if interp=='bilinear': sm='smooth'
plotname=namedir+varnames[loc1]+varnames[loc2]+var[pos1]+'_'+var[pos2]+'_'+var[pos3]+sm+'.png'

pp.figure(figsize=(9.7,6))
pp.title('Spectral gap by scattering theory for '+varnames[pos1]+'='+varval1+' , '+varnames[pos2]+'='+varval2+' and '+varnames[pos3]+'='+varval3)
pp.xlabel(varnames[loc2])
pp.ylabel(varnames[loc1])


pp.imshow(nnnfresult,cmap='hot',interpolation=interp,extent=[minx,maxx,miny,maxy],aspect='auto') 
pp.colorbar()
pp.savefig(plotname)

os.system('mv '+plotname+' '+namedir)














