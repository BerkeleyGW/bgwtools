#!/bin/bash
awk '{print $2,$3,$4,$1,$5}' $1 | column -t > ${1}_LDA.dat
awk '{print $2,$3,$4,$1,$6}' $1 | column -t > ${1}_GW.dat
awk '{print $2,$3,$4,$1,$7}' $1 | column -t > ${1}_diff.dat
