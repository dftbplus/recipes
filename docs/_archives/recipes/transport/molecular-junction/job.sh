#!/bin/bash 

BIN=$2

if [[ $1 == *"c"* ]]; then

 sed 's/ContactId = "Source"/ContactId = "Source"/' dftb_in.hsd_contacts > dftb_in.hsd

 $BIN 
 
 sed 's/ContactId = "Source"/ContactId = "Drain"/' dftb_in.hsd_contacts > dftb_in.hsd

 $BIN 

fi

if [[ $1 == *"t"* ]]; then
 
 EF_source=`grep "Fermi level" shiftcont_source.dat | awk '{print $6}'`
 EF_drain=`grep "Fermi level" shiftcont_drain.dat | awk '{print $6}'`

 sed -e 's/@EF_s/'${EF_source}'/' -e 's/@EF_d/'${EF_drain}'/' dftb_in.hsd_transport > dftb_in.hsd
 
 $BIN 

fi


