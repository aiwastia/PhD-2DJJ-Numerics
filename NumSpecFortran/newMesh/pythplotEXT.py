import numpy as np
from matplotlib import pyplot as pp
import os
import sys

#to run as follow: python pythplotEXT.py dataname param1 param2 param3 param4

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
var=sys.argv[2:6] #var[i] is : or [number]
abuff=sys.argv[2:6] #copy to avoid overwriting var

loc1=abuff.index(':')
abuff[loc1:loc1+1]=[]
loc2=abuff.index(':')+1
llist=[0,1,2,3]
[pos1,pos2]=[x for x in llist if ((x!=loc1) and (x!=loc2))]
var1=int(var[pos1])
var2=int(var[pos2])

#collect size of data
testtxt=not os.path.isfile(namedir+'/param_JJ_2.txt')
if testtxt: os.system('cp '+namedir+'/param_JJ_2 '+namedir+'/param_JJ_2.txt')
ff=open(namedir+'/param_JJ_2.txt',"r")
fflines=ff.readlines()
ff.close()

nstepphi=fflines[8].split('\t')[1] #reminder: in python, remove 1 to the index of line...
nstepB=fflines[9].split('\t')[2]
nstepmuSC=fflines[10].split('\t')[1]
nstepmu=fflines[11].split('\t')[1]

nstepphi=int(nstepphi[1:-1])
nstepB=int(nstepB[1:-1])
nstepmuSC=int(nstepmuSC[1:-1])
nstepmu=int(nstepmu[1:-1])

nstep=[nstepphi,nstepB,nstepmuSC,nstepmu]

#controls the validity of input
if ((var1<1) or (var2<1)):
	print("Indices start at 1.")
	exit()
if ((var1>nstep[pos1]) or (var2>nstep[pos2])):
	print("Maximal values for indices are:",nstep)
	exit()

##extract relevant data for the specific diagram as stated in the input
#collect the indices
lengths=[nstep[1]*nstep[2]*nstep[3],nstep[2]*nstep[3],nstep[3],1]
start=1+(var1-1)*lengths[pos1]+(var2-1)*lengths[pos2] #position in terms of lines
indexsublist=[-1+start+jj*lengths[loc2]+kk*lengths[loc1] for kk in range(nstep[loc1]) for jj in range(nstep[loc2]) ]#if jj>10] #-1 for indexing & for-loops as if written with indents, ie the first changes more slowly###############################################################

#collect all the values of gap + values of constants
f=open(namedir+'/'+namedata,"r")
flines=f.readlines()
fresult=[]
for x in flines:
    fresult.append(x.split('\t')[4])
f.close()
fresult[0:1]=[] #removes the title
nfresult=[float(fresult[i]) for i in indexsublist]

varval1=flines[indexsublist[0]+1].split('\t')[pos1]
varval2=flines[indexsublist[0]+1].split('\t')[pos2]
minx=round(float(flines[indexsublist[0]+1].split('\t')[loc2]),2)
miny=round(float(flines[indexsublist[0]+1].split('\t')[loc1]),2)
maxx=round(float(flines[indexsublist[-1]+1].split('\t')[loc2]),2)
maxy=round(float(flines[indexsublist[-1]+1].split('\t')[loc1]),2)
if minx==maxx:
	maxx=maxx+0.1
	minx=minx-0.1
if miny==maxy:
	maxy=maxy+0.1
	miny=miny-0.1
print(minx,maxx,miny,maxy)

#reshape to plot
nnfresult=np.asarray(nfresult).reshape(nstep[loc1],nstep[loc2])##########################################################################
nnnfresult=np.ndarray.tolist(nnfresult)
nnnfresult.reverse() #to flip the plot (phi increasing from bottom to top)

#name of file
varnames=['Phi','B','muSC','mu']
interp='none' #interpolation = none or bilinear

sm=''
if interp=='bilinear': sm='smooth'
plotname=namedir+varnames[loc1]+varnames[loc2]+var[pos1]+'_'+var[pos2]+sm+'.png'

pp.figure(figsize=(9.7,6))
pp.title('Spectral gap by scattering theory for '+varnames[pos1]+'='+varval1+' and '+varnames[pos2]+'='+varval2)
pp.xlabel(varnames[loc2])
pp.ylabel(varnames[loc1])


pp.imshow(nnnfresult,cmap='hot',interpolation=interp,extent=[minx,maxx,miny,maxy],aspect='auto') 
pp.colorbar()
pp.savefig(plotname)

os.system('mv '+plotname+' '+namedir)














