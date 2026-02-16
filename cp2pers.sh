# Submit multiple wcsim parallel job batches
# James Minock (code based on other people's code as usual)

export SCRIPT_PATH=/exp/annie/app/users/jminock/grid

#jobsub_submit -G annie --dag file://${SCRIPT_PATH}/wcsim.dagnabbit

QUEUE=long

#OUTPUT_FOLDER=/pnfs/annie/scratch/users/jminock/ntuple
OUTPUT_FOLDER=/pnfs/annie/persistent/users/jminock/v1_3_3_world_stv_ntuples
OUTPUT_DATA=/exp/annie/data/users/jminock/temp_add_branches
#OUTPUT_TA=/pnfs/annie/persistent/analysis/v1.3.0/MC/tank_tilt_shift/

if [ ! -d ${OUTPUT_FOLDER} ]
then
    echo "${OUTPUT_FOLDER} does not exist"
    mkdir -p $OUTPUT_FOLDER
fi

redoruns=(1614 1629 1641 1650 1661 1679 1690 1695 1724 1805 1816 1822 1829 1831 1834 1837 1839 1842 1856 1885 1887 1888 1894 1903 1904 1905 1911 1916 1933 1934)

#there are 20000 events per GENIE and ANNIEDirt file
#only 1000 events per wcsim can be run on the grid in a reasonable amount of time
#loop through batches of 20 parallel jobs that use the same inputs files
for RUN in {3500..3999} #skip 4661
do
    if [[ -f ${OUTPUT_FOLDER}/PhaseIITree_0.${RUN}.0.root ]]
    then
        rm -rf ${OUTPUT_FOLDER}/PhaseIITree_0.${RUN}.0.root
    fi
    if [[ ! -f ${OUTPUT_DATA}/PhaseIITree_0.${RUN}.0.root ]]
    then
        continue
    fi
    ifdh cp -D ${OUTPUT_DATA}/PhaseIITree_0.${RUN}.0.root ${OUTPUT_FOLDER}/
    rm -rf ${OUTPUT_DATA}/PhaseIITree_0.${RUN}.0.root
done


