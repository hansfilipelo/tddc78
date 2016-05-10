#!/bin/bash

for i in {1..500}
do 
	if [[ $(./laplsolv |awk '{print $6}' | tail -n 1) == 0.49886012620009207 ]]
	then 
		echo "Correct"
	else
		echo "incorrect"
	fi
done
