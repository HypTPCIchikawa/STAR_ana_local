#!/bin/bash
set -euo pipefail
cd /gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer
bash script/singularity_root_batch.sh config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_histOnlyNoJointSparse.yaml root4star -l -b -q '/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_audit_20261111_001/source/ExportMergingHistograms.C+("/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_1CE33994716F7896423A80AB994BAB14_merge.root","/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/track_merging_audit_20261111_001/histograms.json")'
