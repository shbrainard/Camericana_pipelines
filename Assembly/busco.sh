#!/bin/bash
INPUT=$1

#*** Busco ***
BUSCO="$UNDERLYINGPATH/anaconda/envs/busco4.1.4_py3/bin/busco"
$BUSCO --cpu 40 -i $INPUT -l "auto" -o "${BASE%.fa}" -m "genome" &> "${BASE%.fa}.busco.log" ;