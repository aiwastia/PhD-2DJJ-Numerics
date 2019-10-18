import os

#collect name
ff=open('infile_JJ_2',"r")
fflines=ff.readlines()
ff.close()

ffresult=fflines[1].split('\t')[1]
namedata=ffresult[2:-2]

#create directory to store the files
namedir=namedata[0:len(namedata)-4]

testdir=os.path.isdir(namedir)
while testdir:
	ow=str(input('Directory already exists. Do you want to erase previous data?(n*/y, with single quotation marks) '))
	if ow=='y':
		os.system('rm -i -r '+namedir)
	else:
		namedir=input('Provide a new name for the data (wo/ ext): ')
		oldname=namedata
		namedata=namedir+'.txt'
		os.system('mv '+oldname+' '+namedata)
	testdir=os.path.isdir(namedir)

namedata2=namedir+'.dat'
os.system('mkdir '+namedir)
os.system('cp infile_JJ_bands_1 '+namedir+'/param_JJ_1')
os.system('cp infile_JJ_2 '+namedir+'/param_JJ_2')
os.system('mv '+namedata+' '+namedir)
os.system('cp test1_1_1_1.dat '+namedir+'/'+namedata2)


