#!/bin/bash

# prompt for password
#read -s -p "Password: " password

# download all my outdata
#wget -r -nH --cut-dirs=5 -nc ftp://trevorh2:{$password}@cc-login.campuscluster.illinois.edu//assimilation-cfr/outdata/pointwise

scp -r trevorh2@cc-login.campuscluster.illinois.edu:/home/trevorh2/assimilation-cfr/outdata/pointwise ..




