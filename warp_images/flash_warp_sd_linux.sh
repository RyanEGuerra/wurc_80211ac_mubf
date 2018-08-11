#!/bin/bash
# sd_flasher.sh
# 
# Linux-only!
# Automatically flash SD cards for WARPv3. Detects any .bin files in
# the local directory called "warp_configs" and lets you choose one.
# Handles all the magic memory location numbers for a single-program
# SD card.
#
# Not complete: doesn't handle multiple SD program locations. See
# http://warp.rice.edu/trac/wiki/howto/SD_Config for more details.
#
# This is also greedy: it assumes the SD card appears as /dev/sdb
# which is NOT always true!! Especially on multiple-HD systems,
# this could totally overwrite one of your drives!!
#
# (c) ryan@guerra.rocks 2015
# http://www.apache.org/licenses/LICENSE-2.0

TARGET="/dev/sdc"
CUR_DIR="`dirname \"$0\"`"

echo "sd_flasher.sh utility v1.0"
echo "This will totally nuke everything on $TARGET. Are you sure? [Y/N]"
read RESPONSE

if [[ "$RESPONSE" == "y" || "$RESPONSE" == "Y" ]]; then
  FILES=`ls $CUR_DIR/warp_configs/*.bin`
  COUNTER=-1
  
  echo "Choose a flash file..."
  for F in $FILES
  do
    let COUNTER=COUNTER+1
    FILE[$COUNTER]=$F
	PFILE=`basename $F`
    echo "($COUNTER) $PFILE"
  done
  read NUM
  
  echo "sudo dd bs=512 seek=131072 if=${FILE[$NUM]} of=$TARGET"
  sudo dd bs=512 seek=131072 if=${FILE[$NUM]} of=$TARGET

  echo
  echo "Don't forget to unmount the SD card before removing!"
else
  echo "Exiting w/o doing anything..."
fi

echo "Done!"
exit 0

