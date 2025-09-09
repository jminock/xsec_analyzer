#!/bin/bash

# Number of expected command-line arguments
num_expected=2

STV_PATH=/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples

echo "Combining to make 40k"
hadd ${STV_PATH}/PhaseIITree_40k_stv_ntuple.root ${STV_PATH}/PhaseIITree_0.0.*.root ${STV_PATH}/PhaseIITree_0.1.*.root

echo "Combining to make 400k"
hadd ${STV_PATH}/PhaseIITree_400k_stv_ntuple.root ${STV_PATH}/PhaseIITree_0.?.*.root ${STV_PATH}/PhaseIITree_0.1?.*.root

echo "Combining to make 1mil pt1"
hadd ${STV_PATH}/PhaseIITree_1mil_stv_ntuple_1.root ${STV_PATH}/PhaseIITree_0.?.*.root ${STV_PATH}/PhaseIITree_0.1?.*.root ${STV_PATH}/PhaseIITree_0.2?.*.root ${STV_PATH}/PhaseIITree_0.3?.*.root ${STV_PATH}/PhaseIITree_0.4?.*.root 

echo "Combining to make 1mil pt2"
hadd ${STV_PATH}/PhaseIITree_1mil_stv_ntuple_2.root ${STV_PATH}/PhaseIITree_0.5?.*.root ${STV_PATH}/PhaseIITree_0.6?.*.root ${STV_PATH}/PhaseIITree_0.7?.*.root ${STV_PATH}/PhaseIITree_0.8?.*.root ${STV_PATH}/PhaseIITree_0.9?.*.root 

echo "Combining to make 1mil pt3"
hadd ${STV_PATH}/PhaseIITree_1mil_stv_ntuple_3.root ${STV_PATH}/PhaseIITree_0.10?.*.root ${STV_PATH}/PhaseIITree_0.11?.*.root ${STV_PATH}/PhaseIITree_0.12?.*.root ${STV_PATH}/PhaseIITree_0.13?.*.root ${STV_PATH}/PhaseIITree_0.14?.*.root 

echo "Combining to make 1mil pt4"
hadd ${STV_PATH}/PhaseIITree_1mil_stv_ntuple_4.root ${STV_PATH}/PhaseIITree_0.15?.*.root ${STV_PATH}/PhaseIITree_0.16?.*.root ${STV_PATH}/PhaseIITree_0.17?.*.root ${STV_PATH}/PhaseIITree_0.18?.*.root ${STV_PATH}/PhaseIITree_0.19?.*.root 

echo "Done! Don't forget the NuWro samples and counting POT!"

