#!/bin/bash
set -euo pipefail
qm_project=/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer
qm_root=/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/root-5.34.38-eyawaunrtayjlt3aqon5claakm3c7xl7
qm_star=/star/nfs4/AFS/star/packages/SL24y
qm_out="$qm_project/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_quickcheck_20261111_001"
/usr/bin/g++ -std=c++11 -O2 -fPIC -shared $("$qm_root/bin/root-config" --cflags) -I"$qm_project/include" -I"$qm_project/src/third_party/yaml-cpp/include" -I"$qm_star/.sl73_x8664_gcc485/include" -I"$qm_star/StRoot" -I"$qm_star/StRoot/StChain" -I"$qm_star/StRoot/St_base" -I"$qm_star/StRoot/StarClassLibrary" "$qm_out/source/QuickMergingCheck.cxx" -o "$qm_out/source/libQuickMergingCheck.so" -L"$qm_project/lib" -lStarAnaConfig -L"$qm_star/.sl73_x8664_gcc485/lib" -lStChain -lStPicoEvent -lStarClassLibrary $("$qm_root/bin/root-config" --libs) -Wl,-rpath,"$qm_project/lib" -Wl,-rpath,"$qm_star/.sl73_x8664_gcc485/lib"
