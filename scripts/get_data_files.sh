#!/bin/sh
# Download all data files for the CECAM MDAnalysis tutorial
#
# Shared top-level dropbox folder:
# https://www.dropbox.com/sh/ln0klc9j7mhvxkg/AAB0gMcPPsrDhdVrM2PWmopXa?dl=0

# The curl trick for getting zipped folders is from
# http://stackoverflow.com/questions/21322614/use-curl-to-download-a-dropbox-folder-via-shared-link-not-public-link

# trajectories folder:
DROPBOX_FOLDER='https://www.dropbox.com/sh/am6y00kac8myihe/AABDiQI28fWnRZueQTT7W2s1a?dl=1'

# temporary zip file
ZIPNAME='mdatrj.zip'

echo "-- Download zipped trajectories/ folder from Dropbox..."
curl -o $ZIPNAME -L "$DROPBOX_FOLDER" || { echo "EE curl failed, stopping"; exit 1; }
unzip $ZIPNAME && rm $ZIPNAME




