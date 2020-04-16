#!/bin/bash

for f in $(cat $1); do ./offline_resolution.sh $f; done

