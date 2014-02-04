#/bin/tcsh
if ( $1 == "" ) then
	echo "ERROR: Please input work folder name."
	echo "Ex: ./checkEOSJobs.csh [name]"
	exit
endif
if ( ! ( -e $1 ) ) then
	echo "ERROR: Here is no work folder name $1 "
	exit
endif

cmsenv
set nowPath=`pwd`
cd $1
	rm -f tmp_.log
	set sampleName=`cat datasetList.txt | grep -v '#' | awk '{print $1}' | sed 's/^\///g' | sed 's/\//__/g'`
	foreach sample($sampleName)
		touch tmp_.log
		set i=0
		set notDone=0
		echo "==================================================================================="
		echo "$sample"
		set jobNum=`ls -l $sample/input | grep '.sh' | wc -l`
		set killedJobs=`grep Killed $sample/output/*.log | sed 's/.*job_\(.*\)\.sh.*/\1/g'`
		set killedNum=`echo $killedJobs | wc -w`
		set doneJobs=`/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls eos/cms/store/user/jtsai/bpTobH/backgroundEstimationSkim/$1/$sample | grep root | sed 's/bprimeTobH_\(.*\)\.root/\1/g'` 
		set doneNum=`echo $doneJobs | wc -w`	
		set realDoneNum=`echo $doneNum"-"$killedNum | bc`	
		echo "Status: $realDoneNum/$jobNum"
		if ( $realDoneNum == 0 ) then
			echo "Nothing output..."	
		else if ( $realDoneNum == $jobNum ) then
			echo "Done!"	
		else
			set done=0
			foreach job($doneJobs)	
				if ( $i == $job ) then
					set done=1	
				endif	
			end
			if ( $done == 0 ) then
				echo $i >> tmp_.log
				@ notDone++ 
			endif	
			@ i++
		endif
		if ( $notDone != 0 ) then
			set notDone=`cat tmp_.log`	
			echo "Not done: "$notDone 
		endif
		if ( $killedNum != 0 ) then
			echo "Killed jobs: "$killedJobs 
		endif
		rm -f tmp_.log
	end
	rm -f tmp_.log
cd -

