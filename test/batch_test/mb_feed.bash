#!/bin/bash
echo $(date) > output.log
echo $(date) > short.log
echo "$1" | grep .*.nex 
if [ "$?" -eq "0" ];
then
all_nex_files="$1"
else
all_nex_files=$(ls $1*.nex)
fi
for i in $all_nex_files
do
grep_res=$(grep -i "data;" $i)
if [ "$grep_res" != "" ];
then
(
    echo "execute $i;"
    grep_res=$(grep -i "lset" $i)
    if [ "$grep_res" = "" ]; 
    then
	echo "	lset rates=invgamma nst=6;"
    fi

    echo "	mcmc ng=200 append=no checkfr=100 file=crap;
	mcmc ng=300 append=yes file=crap;
	sumt;
	sump;
	quit;
"
) | ./mb >> output.log 2>&1
# at string abouve both error and stdout redirected to file

#check if ./mb exit with error
if [ "$?" -ne "0" ]; then
  echo "Error in file $i"
  echo "Error in file $i" >> short.log
else
  echo "$i checked"
  echo "$i checked" >> short.log
fi

else
  echo "Skip $i"
  echo "Skip $i" >> short.log
fi
done
exit
