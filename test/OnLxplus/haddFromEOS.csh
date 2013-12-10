#!/bin/tcsh
cmsenv
if ( $3 == "" ) then
	echo "Please input FolderCardm, EOS_Path and output folder"
	echo " EX: source haddFromEOS.csh [FolderCard] [eos/cms/store/user/jtsai/...] [output]"
	exit
endif
if ( !( -e $3 ) ) then
	echo "ERROR: Please check the output folder exist or not"
	exit
endif
set eosPath=`echo $HOME/$2`
if ( !( -e $eosPath ) ) then
	echo "ERROR: Please check eos if mount in home director"
	echo "		EX: eosmount ~/eos"
	exit
endif
set eosPath_=`echo $eosPath | sed 's/\//\\\//g'`
set fs=`cat $1`
foreach f($fs)
	set roots=`eos ls -l $2/$f | awk '{print $9}' | sed "s/^/$eosPath_\/$f\//g"`
	#echo $roots
	set output=`echo $3/$f'.root'`
	echo "hadd $f to $output..."
	hadd -f $output $roots	
end

