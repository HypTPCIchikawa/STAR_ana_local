#!/bin/tcsh

# spack 由来を無効化（最重要）
setenv LD_LIBRARY_PATH ""
unsetenv DYLD_LIBRARY_PATH
unsetenv ROOTSYS
unsetenv PYTHONPATH
unsetenv CMAKE_PREFIX_PATH

# STAR を読み直し
source /afs/rhic.bnl.gov/star/packages/setup.csh
starver SL24a
rehash
