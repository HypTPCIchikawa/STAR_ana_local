# Λ–4He高k*ベースラインとME/k*調査

**高k*の低下は確認できたが、その原因は未特定。k*の静止系変換・電荷補正・単純なME fillの誤りは確認されなかった。0.4–0.5を物理的plateauとして認定しない。**

- [図付き調査PDF](/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/cf_mixing_kstar_audit_20261102_001/lambda_he4_mixing_kstar_audit_20261102.pdf)

## Finding and normalization status

The decreasing high-k* baseline in the restored report p.18 is confirmed. It occurs before purity fitting and remains after the adopted Poisson-L correction. The interval 0.40-0.50 GeV/c is only a reference scale for comparisons; this audit does not establish it as a physical normalization plateau.

No wrong-PRF, missing factor of two, wrong boost sign, He4 charge-scaling or simple mixed-pair filling bug was found in the actual submitted 04BF implementation. This is a bounded code/numerical audit, not proof that the mixed-event reference describes the physical uncorrelated baseline. The origin of the remaining trend is not identified.

The sample is the same 2021 P24iy production, SUMS 04BF2DAE47164B6F040E3A460E539289. Existing ROOT inputs, PID, Maker/configuration and previous reports are unchanged. No production job or new physics normalization is adopted.

## k* definition and frozen-production verification

The submitted routine forms q=pA-pB and P=pA+pB, boosts q with -P.BoostVector(), and returns half the magnitude of its spatial vector. In the pair rest frame, pA*=-pB*, so this equals the momentum magnitude of either particle, including unequal particle masses. Fixed-target laboratory input requires no preliminary CM shift.

SE and ME call the same ComputeKStar and use the same candidate four-vectors. Phi-deuteron in the shared StFemtoMaker uses this same pair routine. The separate StCandidateFemtoMaker and Legacy mapping helper implement the same PRF convention.

He4 starts from primary-track rigidity momentum and receives the charge-two momentum correction once. Its energy uses mass 3.72742 GeV/c^2. The KF Lambda pairing four-vector uses configured mass 1.115683 GeV/c^2; reconstructed invariant mass is stored separately for mass fits.

Six production functions, including both pair loops, pool binning, event storage and MakeLambdaEvent, are byte-identical to the current functions. The full Maker diff changes only generic non-Lambda bad-run/pileup guards; Lambda events return through their earlier dedicated branch.

## Mixed-event bookkeeping and data checks

Submitted mixing is bufferAll, both directions, at most20 buffered events, nine cent9 bins, ten Vz bins over198-202cm, and one event-plane bin. Every eligible candidate pair is attempted; no capped sampling or inverse-sampling correction applies in this run.

The current event is mixed before it is stored. Its eventIndex increases for each input event. Shared-track checks compare eventIndex and trackIndex together. Forward current-Lambda x buffered-nucleus and reverse buffered-Lambda x current-nucleus have independent eligibility. Base channels are processed once; their helper fans out signal/sideband classes without repeating mixing.

The absence of a dedicated event Vz cut made edge-bin clamping a candidate concern, but the recorded mixingVzClampedEvents counter is zero. The coarse Vz QA bin201-202.5cm crosses the mixing boundary: its contents must not be counted as events outside202cm. The current evidence does not support out-of-range clamping as the cause.

The pool key has cent9/Vz but no run, field-sign or measured event-plane matching. Whether those omissions cause a bias requires event-level comparisons; a missing key is not by itself a demonstrated bug. Completely candidate-empty events are not stored; this alone does not prove a shape bias.

## What the measured controls exclude

The decrease remains in individual cent9 samples. For10-20%, the raw0.70-0.80 interval is0.87625 after the0.40-0.50 reference scaling. At k*=0.75 its raw point is0.88517, and its Method5 L point is0.88137. Combining different cent9 bins cannot be the sole cause.

An explicitly diagnostic ME exposure match multiplies each fine-cent denominator by its own SE_norm/ME_norm before summation. The group20-80% raw point at0.75 changes only from0.88136 to0.87932. This does not recover a flat baseline and is not adopted as a physics correction.

Physical He3 and the same He3 tracks under the He4 mass hypothesis are additional raw controls. For0-20% at0.75 they are0.96128 and0.95840, while selected He4 is0.90766. A trend is also present in the control population, but changing the mass hypothesis alone does not reproduce the full He4 effect. These are not purity-corrected null measurements.

The L-extraction raw counts use the documented whole-mass-bin window, while species-control histograms use the exact Maker candidate gate. The maximum absolute difference of their He4 normalized curves over0.30-0.80 is0.00102, much smaller than the broad deficit. This gate difference is kept explicit.

Point errors are statistical and rely on the stored/count-only assumptions. L-derived CF errors retain the existing covariance treatment. Descriptive slope fits omit event/mixed-pair covariance; no calibrated discovery significance is claimed.

## Unresolved cause and the next discriminating measurements

A theoretical femtoscopic component near unity does not by itself make the experimental SE/ME ratio flat: the latter also tests whether the event ensembles and their conditional single-particle/pair distributions match. This audit establishes the trend, not its physical origin.

Next compare event-level candidate momentum distributions and ME closure within finer multiplicity and Vz classes, run/field classes, and measured event-plane classes. Test eventwise conservation/flow baselines and detector/TOF acceptance dependence using a controlled same-candidate remix. Hold PID, Lambda mass selection and event cuts fixed while changing one mixing condition at a time.

The stored projected histograms cannot recreate arbitrary event-level remixes. A candidate/event record that includes both Lambda and nuclear four-vectors and event-class identifiers is needed for that comparison. No full reprocessing is launched by this audit. Choosing a different normalization interval alone is not a remedy for a sloping baseline.
