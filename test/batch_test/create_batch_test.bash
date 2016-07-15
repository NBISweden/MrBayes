#!/bin/bash
j=0
all_nex_files=$(ls $1*.nex)
(
echo "#NEXUS 

begin mrbayes;
	set swapseed=1 seed=1 nowarn=yes autoclose=yes;
"	
for i in $all_nex_files
do

grep_res=$(grep -i "data;" $i)
if [ "$grep_res" != "" ];
then

    let j++
    echo "	[Data set # $j]
	execute $i;"
    grep_res=$(grep -i "lset" $i)
    if [ "$grep_res" = "" ]; 
    then
	echo "	lset rates=invgamma nst=6;"
    fi

    echo "	mcmc ng=1000 checkfr=100 file=crap;
	mcmc ng=2000 append=yes file=crap;
	sumt;
	sump;
"
fi
done
echo "
end;
"
) > batch_test.batch
exit
