#/bin/bash

if [ $# -lt 2 ]; then
	echo $0 pamFoldDir prefix
	exit 1
fi

pamFoldDir=$1
prefix=$2
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

echo pamFoldDor=$pamFoldDir prefix=$prefix


for i in $pamFoldDir/${prefix}*.bin; do
	$SCRIPTPATH/JACKIE.sortToBed $pamFoldDir/A.ref.txt $i ${i/.bin/}.bed 1 0
done
