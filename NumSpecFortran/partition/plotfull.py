import os
import sys

#extract name of data
namedir=sys.argv[1]
namedata=namedir+'.txt'

testdir=not os.path.isdir(namedir)
testfile=os.path.isfile(namedata)
if testdir:
	os.system('echo ''No data available with this name.'' ')
	exit()

var=sys.argv[2:7] #var[i] is : (ordinate axis), :: (abscissa axis) or [number]
abuff=sys.argv[2:7] #copy to avoid overwriting var

loc1=abuff.index(':')
#abuff[loc1:loc1+1]=[] #removed if 2 different symbols
loc2=abuff.index('::')#+1
llist=[0,1,2,3,4]
[pos1,pos2,pos3]=[x for x in llist if ((x!=loc1) and (x!=loc2))]
len1=int(var[pos1])
len2=int(var[pos2])
len3=int(var[pos3])

for l in range(len1):
	for ll in range(len2):
		for lll in range(len3):
			cplot=[0]*5
			cplot[loc1]=': '
			cplot[loc2]=':: '
			cplot[pos1]=str(l+1)+' '
			cplot[pos2]=str(ll+1)+' '
			cplot[pos3]=str(lll+1)+' '
			codeplot=''
			for kk in range(5):
				codeplot=codeplot+cplot[kk]
			os.system("echo '"+namedir+' '+codeplot+"'")
			os.system('python pythplotEXT.py '+namedir+' '+codeplot)
