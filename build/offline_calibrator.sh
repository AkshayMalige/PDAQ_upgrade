#!/bin/bash

fn=$1

./pdaq_calibrater /storage2/Data/cosy_raw/$fn /storage2/Data/cosy_cal/$(basename $fn )_CAL.root 1000000000

