#!/bin/bash

fn=$1

./pdaq_resolution /home/akshay/PDAQ/PDAQ_upgrade/build/cosy_dtcal/$fn /home/akshay/PDAQ/PDAQ_upgrade/build/cosy_res/$(basename $fn )_DT_CAL.root 10000
