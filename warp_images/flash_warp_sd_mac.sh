#!/bin/bash
# Flashing utility for WARP bin files. Tested on Mac OS X Mavericks.
#
# To try not to nuke my computer, this checks to make sure that the size
# of the target SD card is actually 2.0 GB. This is a MAC-ONLY implementation
# because it uses diskutil.
#
# Sometimes the OS doesn't mount the SD card very fast and you will run this script,
# select the image files, and it will say that the drive cannot be found. Try again.
#
# (c) ryan@guerra.rocks 2015
# http://www.apache.org/licenses/LICENSE-2.0

# !!!!!!!!!!!!!!!!  EDIT THIS PATH WITH CARE !!!!!!!!!!!!!!!!!!!
# !!!!!! THIS DRIVE WILL BE THE ONE THAT IS OVERWRITTEN !!!!!!!!
disk_name="/dev/disk2"
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

current_dir=`pwd`
bin_dir="warp_configs"
sd_size=`diskutil info $disk_name | grep "Total Size:" | awk '{print $3}'`

echo "sd_flasher.sh utility v2.0"
echo "This will totally nuke everything on ${disk_name}. Are you sure? [Y/N]"
read RESPONSE

# If the user says no, then don't do anything...
if [[ "${RESPONSE}" == "y" || "${RESPONSE}" == "Y" ]]; then
  FILES=`ls ${current_dir}/${bin_dir}/*.bin`
  
  echo "Choose a flash file for SLOT 0..."
  COUNTER=-1
  for F in $FILES
  do
    let COUNTER=COUNTER+1
    FILE[$COUNTER]=$F
	PFILE=`basename $F`
    echo "($COUNTER) $PFILE"
  done
  read NUM
  
  echo "Choose a flash file for SLOT 1..."
  COUNTER=-1
  echo "(N) NO FILE"
  for F in $FILES
  do
    let COUNTER=COUNTER+1
    FILE[$COUNTER]=$F
	PFILE=`basename $F`
    echo "($COUNTER) $PFILE"
  done
  read NUM2
  
  # Check the size of the SD card for safety...
  if [[ "${sd_size}" == "2.0" ]]; then
	# do programming
	diskutil unmountDisk ${disk_name}
	
	echo "sudo dd bs=512 seek=131072 if=${FILE[$NUM]} of=${disk_name}"
	sudo dd bs=512 seek=131072 if=${FILE[$NUM]} of=${disk_name}
	
	# do slot2 programming
	if [[ "${NUM2}" == "N" ]]; then
	  echo "Writing nothing in SLOT 1..."
	else
	  #actually do it
	  sleep 1
	  diskutil unmountDisk ${disk_name}
	  echo "sudo dd bs=512 seek=163840 if=${FILE[$NUM2]} of=${disk_name}"
	  sudo dd bs=512 seek=163840 if=${FILE[$NUM2]} of=${disk_name}
	fi
	
  else
    # Don't program
    echo "WARNING: the default disk $disk_name is NOT 2.0 GB: $sd_size"
    echo "         Nothing will be done."
    exit 1
  fi

else
  echo "Exiting without doing anything..."
fi

# eject the sucker
echo "Ejecting SD card..."
sleep 1
diskutil eject $disk_name

echo "Done!"
exit 0


