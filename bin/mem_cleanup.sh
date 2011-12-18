#!/bin/bash
# This script cleans up the shared memory segments that might have remained
# on the machines if ACESIII crashed in the middle of its operation.
# You may see some messages which say "ipcrm: invalid key", which is normal.

echo -e "\n------------Memory cleanup script running-------------\n"

#Check if hosts.out exists
if [[ ! -r hosts.out ]]; then
	echo 'File hosts.out doesnt exist, so exiting.'
	echo -e "\n------Memory cleanup script finished executing--------\n"
	exit 0
fi

exec 3< hosts.out
read MEMID <&3
echo "Memid is $MEMID"
while read HOSTNAME <&3
do 
     echo "Working on host $HOSTNAME"
     ssh $HOSTNAME "ipcrm -M '$MEMID'"
done
exec 3>&-

echo -e "\n------Memory cleanup script finished executing--------\n"
