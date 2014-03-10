#!/bin/bash

SCAF_DIR=/Users/bouchard/Documents/courses/stat547q-sp2013-14/exercises/4/
SRC_DIR=/Users/bouchard/w/conifer
PRJ=conifer
CUR=`pwd`

# backup old
cd $SCAF_DIR
zip -r ${PRJ}-scaffold ${PRJ}-scaffold > /dev/null
mv ${PRJ}-scaffold.zip "backups/`date +%s`.zip"

# back up .git
yes | rm -r git-temp-backup 2> /dev/null
cp -r ${PRJ}-scaffold/.git git-temp-backup

# remove the scaffold
yes | rm -r ${PRJ}-scaffold

# copy back
cp -r $SRC_DIR ${PRJ}-scaffold

# remove non-scaffold .git
yes | rm -r ${PRJ}-scaffold/.git

# move back our .git
mv git-temp-backup ${PRJ}-scaffold/.git

# go over files, removing flagged parts of the code
removeFlaggedCode ${PRJ}-scaffold y

# refresh tutorialj
cd ${PRJ}-scaffold
rm create-scaffold.bash
gradle tutorialj -x test

cd $CUR