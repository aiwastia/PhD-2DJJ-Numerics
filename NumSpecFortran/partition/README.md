Folder to run the spectra. Try not to push the data on the remote server of GitHub but to keep them locally while the generating files are tracked. Technically they should be a copy of gapFocus' files.

When copying the files in a new folder for a new batch:
* ClusterRun.sh < change address to the NEW FOLDER ; change job's name; adapt time
* dirmaker.py
* infile_JJ_2 < change data name ; parameters
* infile_JJ_bands_1 < parameters
* JJ_gap < check that it's compiled from the latest version
* Makefile
* pythplotEXT.py