#!/bin/bash
if [ "$1" = "clean" ];
then
echo "Removing files from ./WithoutMbBlock/ and ./WithMbBlock/"
rm ./WithoutMbBlock/*
rm ./WithMbBlock/*
exit
fi

all_nex_files=$(ls $1*.nex)
for i in $all_nex_files
do
grep_res=$(grep -i "data;" $i)
if [ "$grep_res" != "" ];
then
    grep_res=$(grep -i "mrbayes;" $i)
    if [ "$grep_res" = "" ]; 
    then
	cp $i ./WithoutMbBlock/$i
    else
	cp $i ./WithMbBlock/$i
    fi
fi
done
exit
