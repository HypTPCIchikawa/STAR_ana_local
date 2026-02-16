#!/bin/csh
foreach i (`cat run.list`)
	set job=`condor_q onoritaka | tail | grep query | cut -d' ' -f4`
	while( $job >= 3000 )
		echo "waiting jobs finished...."
		echo `condor_q onoritaka | tail | grep query`
		sleep 15m
		set job=`condor_q onoritaka | tail | grep query | cut -d' ' -f4`
	end
	set RunNumber = ${i}
	star-submit-template -template JobNoritaka.xml -entities run=${RunNumber}
	echo "RunNumber ${RunNumber} has been submitted."
end

