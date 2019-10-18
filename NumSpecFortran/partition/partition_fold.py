import numpy as np
import os
import sys

#to be used as : python partition_fold.py [name of the last subfolder of the series]

#extract information from terminal input
fullinfo=sys.argv[1]
detinfo=fullinfo.split('_')
dataname=detinfo[0]
len1=int(detinfo[1])
len2=int(detinfo[2])
len3=int(detinfo[3])
len4=int(detinfo[4][4:])

#title line of data
gettitle=open(fullinfo+'/'+fullinfo+'.txt',"r")
ff=gettitle.readlines()
gettitle.close()
linetitle=ff[0]

#fetch nstepphi in updated infile_JJ_2 of run subfolder
getnphi=open(fullinfo+'/infile_JJ_2',"r")
Wnphiline=getnphi.readlines()
getnphi.close()
nphiline=Wnphiline[9]

buffinfile2=open('infile_JJ_2_'+dataname,"r")
bufflines2=buffinfile2.readlines()
buffinfile2.close()
buffinfile22=open('infile_JJ_2_'+dataname,"w")
for li in range(len(bufflines2)):
	if li==9:
		buffinfile22.write(nphiline)
	else:
		buffinfile22.write(bufflines2[li])
buffinfile22.close()

#begin writing in the final data file
fulldata=open(dataname+'.txt',"w")
fulldata.write(linetitle)

#save time information from make.out files
tottime=0.
maxtime=0.

#loops over the subfolders to extract and sum up data
for l in range(len1):
	for ll in range(len2):
		for lll in range(len3):
			for llll in range(len4):
				foldname=dataname+'_'+str(l+1)+'_'+str(ll+1)+'_'+str(lll+1)+'_part'+str(llll+1)
				testfile=os.path.isfile(foldname+'/'+foldname+'.txt')

				if testfile: #avoid problems when no data has been computed
					partdata=open(foldname+'/'+foldname+'.txt',"r")
					datablock=partdata.readlines()
					partdata.close()
					datablock[0:1]=[] #removes title
					for r in range(len(datablock)):
						fulldata.write(datablock[r])

					timedata=open(foldname+'/'+'make.out',"r")
					timeblock=timedata.readlines()
					timedata.close()
					Wtime=timeblock[-1]
					time=float(Wtime.split('=')[1])
					tottime+=time
					if time>maxtime:
						maxtime=time
fulldata.close()

#create and fill in folder with the name of the data (provided by infile_JJ_2)
os.system('mkdir '+dataname) #doesn't exist already because of check in partition_gen.py
os.system('mv infile_JJ_2_'+dataname+' '+dataname+'/param_JJ_2.txt')
os.system('mv infile_JJ_bands_1_'+dataname+' '+dataname+'/param_JJ_1.txt')
os.system('mv '+dataname+'.txt '+dataname)

#save data information
with open('TimeData.txt',"w") as timedata:
	timedata.write('Total cumulated time for the run: '+str(tottime)+' min\n')
	timedata.write('Effective time (maximum time among the parallelized jobs): '+str(maxtime)+' min')
os.system('mv TimeData.txt '+dataname)

#remove all the subfolders
for l in range(len1):
	for ll in range(len2):
		for lll in range(len3):
			for llll in range(len4):
				foldname=dataname+'_'+str(l+1)+'_'+str(ll+1)+'_'+str(lll+1)+'_part'+str(llll+1)

				os.system('rm -r '+foldname)











