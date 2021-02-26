#!/bin/bash

#Download the releases from the following dates (wanring -- large files)
#RELEASE_DATE="07-07-2020 14-07-2020 21-07-2020 29-07-2020 03-08-2020 10-08-2020 17-08-2020 24-08-2020 31-08-2020 07-09-2020 14-09-2020 21-09-2020 28-09-2020 05-10-2020 12-10-2020 19-10-2020 26-10-2020" 
#RELEASE_DATE="02-11-2020 10-11-2020 16-11-2020 23-11-2020 30-11-2020 07-12-2020 14-12-2020 21-12-2020 28-12-2020 04-01-2021 11-01-2021"
RELEASE_DATE="04-01-2021 11-01-2021 18-01-2021"

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
