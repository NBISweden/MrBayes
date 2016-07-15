#!/bin/bash
all_nex_files=$(ls $1*.nex)
for i in $all_nex_files
do
sed 's/mcmc /mcmcp /; s/mcmc;/ /;' $i > ./Modified/$i
done
exit
