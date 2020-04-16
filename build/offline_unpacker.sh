#!/bin/bash

fn=$1

./pdaq_unpacker_HADES /storage2/Data/cosy_hld/$fn /storage2/Data/cosy_raw/$(basename $fn )_raw.root 1000000000
