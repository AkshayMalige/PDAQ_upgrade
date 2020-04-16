#!/bin/bash

fn=$1

./pdaq_cluster_finder_cosy /storage2/Data/cosy_cal/$fn /storage2/Data/cosy_track/$(basename $fn )_TRACK.root 1000000


