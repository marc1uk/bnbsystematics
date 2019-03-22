h1. Booster Neutrino Beam Systematics

The BNB systematics framework is based on the work developed by MiniBooNE as described here https://arxiv.org/abs/0806.1449. It also includes all of the subsequent updates. In particular, for charged pions splines are used to propagate the HARP error matrix, and the K+ production errors were updated to incorporate the SciBooNE measurement (https://arxiv.org/abs/1110.0417). 
The code was initially ported from MiniBooNE framework by MicroBooNE into LArSoft EventWeight module. Here the same MicroBooNE code was extracted with all LArSoft/ART dependencies removed and made into a stand-alone module.


h1. Downloading the code

The code is maintained in a git "repository":https://cdcvs.fnal.gov/redmine/projects/booster-neutrino-beamline/repository/systematics associated with this redmine project.  You may clone this repository anonymously and without authentication but in order to push commits back you must be authenticated.  

Authenticated clone (i.e., allows to modify the software and upload your modification on this site, via git push):

<pre>
$ git clone ssh://p-systematics@cdcvs.fnal.gov/cvs/projects/systematics
</pre>

Anonymous clone (no push):
<pre>
$ git clone http://cdcvs.fnal.gov/projects/systematics
</pre>

h1. Setup

To setup the code cd into top directory. Note that you should run source from top directory to get the LD_LIBRARY_PATH and FW_SEARCH_PATH env variables correctly set.
Source the setup script:
<pre>
$ source Setup/setup.sh
</pre>

h1. Build

From top directory cd into build and run:
<pre>
$ cmake ..
$ make install -j3
</pre>

h1. Run

To run the code cd into top directory and run the rwgh.
<pre>
$ bin/rwgh -h
Options:
  -h [ --help ]                    Print help message
  -i [ --input ] arg               Path pattern for input files. Put it in 
                                   quotes or escape *.
  -f [ --fhicl ] arg               Configuration fhicl file.
  -o [ --output ] arg (=hist.root) Output file name.
</pre>

For example to generate histograms using one beam ntuple: 
<pre>
$ bin/rwgh -f fcl/eventweight_microboone.fcl -i april07_baseline_0001.root -o output_0001.root
</pre>

h1. Running on grid

You can use the provided submitJob.py script.
<pre>
$ Scripts/submitJob.py -h
usage: submitJob.py [-h] [-d] -f FHICLFILE -g GROUP -i INPUTPATH [-j JOBID] -n
                    [1-10000] -o OUTPUTPATH

Submit beam sys jobs.

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           Will not delete submission files in the end. Useful
                        for debugging and will only print the submission
                        command on screen.
  -f FHICLFILE, --fhicl-file FHICLFILE
                        Configuration fhicl file.
  -g GROUP, --group GROUP
                        Group used to submit jobs.
  -i INPUTPATH, --input-path INPUTPATH
                        Path to input BooNEG4Beam files.
  -j JOBID, --jobidoffset JOBID
                        Id for running job. If submitting multiple jobs each
                        one is offset by its process number (0..n)
  -n [1-10000]          Number of jobs to submit.
  -o OUTPUTPATH, --output-path OUTPUTPATH
                        pnfs path where to copy final output.
</pre>

Example how to submit 10 jobs:
<pre>
$ Scripts/submitJob.py -i /pnfs/uboone/persistent/uboonebeam/bnb_mc/ -o /pnfs/uboone/scratch/users/zarko/test_bnb_sys -f fcl/eventweight_microboone.fcl -n 10 -g uboone
</pre>

Once the jobs are done running add all output histograms into one root file:
<pre>
$ hadd merged_hist.root /pnfs/uboone/scratch/users/zarko/bnb_sys/*/*.root 
</pre>

h1. Analyze

Simple script that demonstrates building the error matrix is provided in Scripts directory. For example to plot the numu flux and errors in root run:
<pre>
.x Scripts/analyze.c("merged_hist.root","numu")
</pre>

Note that you need to modify the script for other neutrino species to include the right set of systematics.

h1.  Notes

Several systematics reweight using flux histograms (horn current, skin depth, nucleon xsec).
For these systematics neutrino spectrum was calculated at particular location for input histograms and this location should be matched in the dk2nu nuray passed to the WeightCalc.
Right now it is hard-coded to use 1st location (0th is usually used for random direction).
