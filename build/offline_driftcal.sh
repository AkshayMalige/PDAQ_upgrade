#!/bin/bash

fn=$1

./pdaq_drift_cal /home/akshay/PDAQ/PDAQ_upgrade/build/cosy_track/$fn /home/akshay/PDAQ/PDAQ_upgrade/build/cosy_dtcal/$(basename $fn )_DT_CAL.root 1000000000
