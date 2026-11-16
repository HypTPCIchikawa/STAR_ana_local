#!/bin/bash
set -euo pipefail
[[ $# -eq 3 ]] || { echo 'usage: build_runtime.sh RUNTIME ROOTSYS STAR_ROOT' >&2; exit 2; }
mc_runtime="$1";mc_root="$2";mc_star="$3"
/usr/bin/g++ -std=c++11 -O2 -fPIC -shared $("$mc_root/bin/root-config" --cflags) -I"$mc_runtime/include" -I"$mc_star/.sl73_x8664_gcc485/include" -I"$mc_star/StRoot" -I"$mc_star/StRoot/StChain" -I"$mc_star/StRoot/St_base" -I"$mc_star/StRoot/StarClassLibrary" "$mc_runtime/source/QuickMergingCheck.cxx" -o "$mc_runtime/lib/libQuickMergingCheck.so" -L"$mc_runtime/lib" -lStarAnaConfig -L"$mc_star/.sl73_x8664_gcc485/lib" -lStChain -lStPicoEvent -lStarClassLibrary $("$mc_root/bin/root-config" --libs) -Wl,-rpath,'$ORIGIN' -Wl,-rpath,"$mc_star/.sl73_x8664_gcc485/lib"
"$mc_star/.sl73_x8664_gcc485/bin/root4star" -l -b -q "$mc_runtime/source/BuildFrozenLoader.C(\"$mc_runtime\")"
