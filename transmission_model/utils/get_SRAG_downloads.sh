#!/bin/bash

RELEASE_DATE="08-02-2021"

for i in $RELEASE_DATE; do
	#Base 2020 download
	if [ ! -f "INFLUD-"$i".csv" ]; then
		wget "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2020/INFLUD-"$i".csv"
	fi
	#Base 2021
	if [ ${i: -4} -gt 2020 ] && [ ! -f "INFLUD21-"$i".csv" ]; then
		wget "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2021/INFLUD21-"$i".csv"
	fi
done
