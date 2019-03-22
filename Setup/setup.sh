source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setups
setup root v6_12_06a -q e15:prof
setup clhep v2_3_4_6 -q e15:prof
setup dk2nu v01_05_01b -q e15:prof
setup fhiclcpp v4_06_08 -q e15:prof
setup cmake v3_3_2
setup jobsub_client v1_2_6
setup ifdhc v2_1_0
setup boost v1_66_0a -q e15:prof

export LD_LIBRARY_PATH=`pwd`/lib:${LD_LIBRARY_PATH}
export FW_SEARCH_PATH=${FW_SEARCH_PATH:+${FW_SEARCH_PATH}:}`pwd`

