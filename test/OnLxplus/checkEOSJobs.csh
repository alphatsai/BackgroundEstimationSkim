#/bin/tcsh
echo "#############################################"
echo "###                                       ###"
echo "###  ./checkEOSJob.csh [name] <reSubmit>  ###"
echo "###                                       ###"
echo "#############################################"
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
set start=`/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls eos/cms/store/user/jtsai/bpTobH/backgroundEstimationSkim | grep $1`
if ( $start == "" ) then
	echo "Nothing output..."
	exit
endif

cd $1
	set nowPath=`pwd`
	rm -f tmp_.log
	set sampleName=`cat datasetList.txt | grep -v '#' | awk '{print $1}' | sed 's/^\///g' | sed 's/\//__/g'`
	foreach sample($sampleName)
		touch tmp_.log
		set i=0
		set notDone=0
		echo "==================================================================================="
		echo "$sample"
		set killedJobs=`grep Killed $sample/output/*.log | sed 's/.*job_\(.*\)\.sh.*/\1/g'`
		set killedNum=`echo $killedJobs | wc -w `
		set jobNum=`ls -l $sample/input | grep '.sh' | wc -l`
		set doneJobs=`/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls eos/cms/store/user/jtsai/bpTobH/backgroundEstimationSkim/$1/$sample | grep root | sed 's/bprimeTobH_\(.*\)\.root/\1/g'` 
		set doneNum=`echo $doneJobs | wc -w`	
		set realdoneNum=`echo $doneNum'-'$killedNum | bc`	
		echo "Status(root): $doneNum/$jobNum"
		echo "Status(real): $realdoneNum/$jobNum"
		if ( $doneNum == 0 ) then
			echo "Nothing output..."	
		else if ( $doneNum == $jobNum ) then
			echo "Done!"	
		else
			while ( $i < $jobNum )
				set done=0
				#echo $doneJobs
				foreach job($doneJobs)	
					if ( $i == $job ) then
						set done=1	
					endif	
				end
				if ( $done == 0 ) then
					#echo $i
					echo $i >> tmp_.log
					@ notDone++ 
				endif	
				@ i++
			end
		endif
		if ( $notDone != 0 ) then
			set notDonelist=`cat tmp_.log`	
			echo "No root Jobs: "$notDonelist 
		endif
		if ( $killedNum != 0 ) then
			echo "Killed Jobs: "$killedJobs 
		endif
		rm -f tmp_.log

		if ( $2 == 'reSubmit' && $notDone != 0 ) then
			foreach nn($notDonelist)
				mv $nowPath/$sample/output/job_$nn.log $nowPath/$sample
				echo resubmit job_$nn.sh...
				bsub -q 1nh -o $nowPath/$sample/output/job_$nn.log source $nowPath/$sample/input/job_$nn.sh
			end	
		endif
		if ( $2 == 'reSubmit' && $killedNum != 0 ) then
			foreach kn($killedJobs)
				mv $nowPath/$sample/output/job_$kn.log $nowPath/$sample
				echo resubmit job_$kn.sh...
				bsub -q 1nh -o $nowPath/$sample/output/job_$kn.log source $nowPath/$sample/input/job_$kn.sh
			end	
		endif
	end
	rm -f tmp_.log
cd -

