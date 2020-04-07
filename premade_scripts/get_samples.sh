#! /usr/bin/env bash


# The first argument of this command is the name of the final VCF containing all the individuals
VCF_FILE=$1

# The second argument of this command, $2 is the motif in the sample name that tell if it is resistant or not for the treatement you want to look at
RES=$2

# This command read the samples names in $VCF_FILE pass them (using |) to grep that print only the sample that contains the motif $RES
# The names returned are then directed (with >) to be written in the file named  ${RES}.resistant.sampleList.txt (with ${RES} replaced in the real name by the second argument)
bcftools query -l $VCF_FILE | grep $RES > ${RES}.resistant.sampleList.txt
# The only difference with the previous one is -v is grep that do a reverse grep, taking the lines that do not have the motif $RES
bcftools query -l $VCF_FILE | grep -v $RES > ${RES}.wild.sampleList.txt
