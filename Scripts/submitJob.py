#!/usr/bin/env python

import argparse
import os.path
import sys
import tarfile
import re
import subprocess

parser = argparse.ArgumentParser(description='Submit beam sys jobs.')

parser.add_argument("-d", "--debug",action='store_true',
                    help="Will not delete submission files in the end. Useful for debugging and will only print the submission command on screen.")
parser.add_argument("-f","--fhicl-file", dest="fhiclfile",
                    required=True,
                    help="Configuration fhicl file.")
parser.add_argument("-g","--group",
                    required=True,
                    help="Group used to submit jobs.")
parser.add_argument("-i", "--input-path", dest="inputpath",
                    required=True,
                    help="Path to input BooNEG4Beam files.")
parser.add_argument("-j", "--jobidoffset", type=int, dest="jobid",
                    default=0,
                    help="Id for running job. If submitting multiple jobs each one is offset by its process number (0..n)")
parser.add_argument("-n", type=int,
                    choices=range(1,10000),
                    metavar="[1-10000]",
                    required=True,
                    help="Number of jobs to submit.")
parser.add_argument("-o", "--output-path", dest="outputpath",
                    required=True,
                    help="pnfs path where to copy final output.")

args = parser.parse_args()

#check for lib and bin
flist=[ "bin/rwgh", "lib/libBNBSysBase.so", "lib/libBNBSysCalc.so", args.fhiclfile ]
for ff in  flist:
    if (not os.path.isfile(ff)):
        print "%s does not exist."%ff
        sys.exit(1)
if (not os.path.isdir("beamData")):
    print "Missing beamData directory. (with input root files)"
    sys.exit(1)


#now create jobfiles_*.tar that is shipped with the job
#tar -cf jobfiles.tar --transform '!^[^/]*/!!' file1 file2 file3
if not os.path.exists(args.outputpath):
    os.makedirs(args.outputpath)
tarfilename="%s/jobfiles_%i.tar"%(args.outputpath,os.getpid())
outtar = tarfile.open(tarfilename, mode='w')
for f in flist:
    outtar.add(f, arcname=os.path.basename(f))
outtar.add("beamData")
outtar.close()

ofstr='''
#!/bin/bash
echo "Running $0 on "$HOSTNAME
echo "Cluster: " ${CLUSTER}
echo "Process: " ${PROCESS}

export JOBID=$((PROCESS+%(jobidoffset)s))
cd ${_CONDOR_SCRATCH_DIR}
mkdir ${JOBID}
cd ${JOBID}
cp $INPUT_TAR_FILE .
tar -vxf `basename ${INPUT_TAR_FILE}`
INPUTFILE=`printf "%%s/april07_baseline_%%04d.root" %(inputpath)s ${JOBID}`
OUTPUTFILE=`printf "hist_april07_baseline_%%04d.root" ${JOBID}`

source setup.sh
ifdh cp $INPUTFILE ./beammc_tmp.root
./rwgh -f %(fhiclfile)s -i beammc_tmp.root -o ${OUTPUTFILE}s
rm rwgh
rm setup.sh
rm libBNBSysBase.so  
rm libBNBSysCalc.so
rm `basename ${INPUT_TAR_FILE}`
rm beammc_tmp.root
cd ..
ifdh mkdir %(outputdir)s
ifdh cp -r ${JOBID} %(outputdir)s/${JOBID}
'''%{'jobidoffset':str(args.jobid),'inputpath':args.inputpath,'outputdir':args.outputpath,'fhiclfile':os.path.basename(args.fhiclfile)}

runjobfname="runjob_%i.sh"%os.getpid()
of=open(runjobfname,'w')
of.write(ofstr)
of.close()

#Create submit command
cmd="jobsub_submit --group=%s -N %i --tar_file_name=dropbox://%s file://%s"%(args.group,args.n,os.path.abspath(tarfilename),os.path.abspath(runjobfname))

if (not args.debug):
    print "Running submit cmd:"
    print cmd
    os.system(cmd)
else:
    print "Would have ran:"
    print cmd

#Delete temp files unless debugging
if (not args.debug):
    os.remove(runjobfname)
    os.remove(tarfilename)
