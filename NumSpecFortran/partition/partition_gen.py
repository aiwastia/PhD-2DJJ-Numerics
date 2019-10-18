import numpy as np
import math
import os
import sys

#EXTRACT information from files to update
ff=open('infile_JJ_2',"r")
fflines=ff.readlines()
ff.close()

ff2=open('ClusterRun.sh',"r")
ff2lines=ff2.readlines()
ff2.close()

Wdataname=fflines[1].split('\t')[1]
dataname=Wdataname[2:-6]

Wphilist=fflines[3].split('\t')[1]
WWphilist=Wphilist[2:-2]

WBlist=fflines[4].split('\t')[2]
WWBlist=WBlist[2:-2]

WmuSClist=fflines[5].split('\t')[1]
WWmuSClist=WmuSClist[2:-2]

Wmulist=fflines[6].split('\t')[2]
WWmulist=Wmulist[2:-2]

WSCmlist=fflines[7].split('\t')[1]
WWSCmlist=WSCmlist[2:-2]

Wpartition=fflines[15].split('\t')[1]
partition=int(Wpartition[1:-1])





##########################################################
# function that cuts apart ONLY a ':' sequence

def seq_gen(list_crypt,index):
	if len(list_crypt.split(':'))>1 :
		start=float(list_crypt.split(':')[0])
		end=float(list_crypt.split(':')[2])
		mid=list_crypt.split(':')[1]
		if len(mid.split('.'))==1 :
			step=0.
			nstep=int(mid)
		else :
			step=float(mid)
			nstep=0

		nstepnew=nstep
		stepnew=step
		if start==end :
				nstepnew=1
				stepnew=0.
		elif start>end :
			print('#Positive parameters required, and MAX >= MIN, in parameter')
			exit()
		else :
			if nstep>1 :
				stepnew=(end-start)/(nstep-1)
			elif nstep==1 :
				stepnew=0.
			elif ((step!=0.)and(step<=(end-start))) :
				nstepnew=int(round((end-start)/step))+1
			else :
				print("#Need a number of steps or a smaller step in parameter ")
				exit()

		Pnstep=nstepnew/partition #integer division (int floor)
		Pstart=start+index*Pnstep*stepnew
		Pend=min(Pstart+(Pnstep-1)*stepnew,end)
		newpartition=int(math.ceil(nstepnew/float(Pnstep))) #avoid integer division

		seq_ind=str(Pstart)+':'+str(stepnew)+':'+str(Pend)
	else :
		seq_ind=list_crypt

	return seq_ind,Pnstep,newpartition
###############################################################

#ISOLATE the discrete values of series for mu, SCm, muSC
listmu=WWmulist.split(',') #each element is a string
listSCm=WWSCmlist.split(',')
listmuSC=WWmuSClist.split(',')
Nmu=len(listmu)
NSCm=len(listSCm)
NmuSC=len(listmuSC)

#UPDATE 'partition' to 'newpartition'
seqB,NB,newpartition=seq_gen(WWBlist,0)

#LOOP that generates all the subfolders and run the bash files
for p in range(Nmu):
	mu=listmu[p]
	for pp in range(NSCm):
		SCm=listSCm[pp]
		for ppp in range(NmuSC):
			muSC=listmuSC[ppp]
			for i in range(newpartition):
				seqB,NB,P=seq_gen(WWBlist,i)

				#TEST availability of subfolders' names
				foldname=dataname+str(p+1)+'_'+str(pp+1)+'_'+str(ppp+1)+'_part'+str(i+1)
				testdir=os.path.isdir(foldname)
				while testdir:
					ow=str(input('Directory already exists. Do you want to erase previous data?(n*/y, with single quotation marks) '))
					if ow=='y':
						os.system('rm -i -r '+foldname)
					else:
						foldname=input('Provide a new name for the data (wo/ ext): ')
					testdir=os.path.isdir(foldname)

				#CREATE and FILL subfolders
				os.system('mkdir '+foldname)
				os.system('cp infile_JJ_bands_1 '+foldname)
				os.system('cp JJ_gap '+foldname)
				
					#with new infile_JJ_2
				with open('infile_JJ_2_tmp', 'w') as infile:
					infile.write(fflines[0])
					infile.write("dataname\t='"+foldname+".txt'\n")
					infile.write(fflines[2])
					infile.write(fflines[3])
					infile.write("ABlist\t\t='"+seqB+"'\n")
					infile.write("AmuSClist\t='"+muSC+"'\n")
					infile.write("Amulist\t\t='"+mu+"'\n")
					infile.write("ASCmlist\t='"+SCm+"'\n")
					infile.write(fflines[8])
					infile.write('Anstepphi\t='+str(100)+'\n') #random choice
					infile.write('AnstepB\t\t='+str(NB)+'\n')
					infile.write('AnstepmuSC\t='+str(NmuSC)+'\n')
					infile.write('Anstepmu\t='+str(Nmu)+'\n')
					infile.write('AnstepSCm\t='+str(NSCm)+'\n')
					infile.write(fflines[16])
				os.system('cp infile_JJ_2_tmp '+foldname+'/infile_JJ_2')
				
					#with new ClusterRun.sh
				fullpath=os.getcwd()
				with open('ClusterRun_tmp', 'w') as Run:
					Run.write(ff2lines[0])
					Run.write('#SBATCH -D '+fullpath+'/'+foldname+'\n')
					Run.write(ff2lines[2])
					Run.write('#SBATCH --job-name='+'P'+str(p+1)+str(pp+1)+str(ppp+1)+str(i+1)+'\n')
					Run.write(ff2lines[4])
					Run.write(ff2lines[5])
					Run.write(ff2lines[6])
					Run.write(ff2lines[7])
					Run.write(ff2lines[8])
					Run.write(ff2lines[9])
					Run.write(ff2lines[10])
					Run.write(ff2lines[11])
				os.system('cp ClusterRun_tmp '+foldname+'/ClusterRun.sh')

				#RUN bash script
				#os.system('sbatch '+foldname+'/ClusterRun.sh')


#CLEAN a little bit
os.system('rm infile_JJ_2_tmp')
os.system('rm ClusterRun_tmp')






				

