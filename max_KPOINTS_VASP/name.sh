#!/bin/bash
for((i=1;i<=230;i++))
do
echo $i
mv KPOINTS_$i KPOINTS_$i.txt 
done
