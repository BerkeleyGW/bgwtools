#!/bin/bash
awk '{if (NR>1){print $2,$3,$4,2,$1}}' $1 | column -t > ${1}.dat
