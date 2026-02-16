#!/bin/csh

starver SL24a

# data check
set runnum=21035003a
root4star -b readPicoDst.C\(\"$runnum.list\",\"$runnum.root\"\)
