#!/usr/bin/env python3
"""Report the frozen-production kinematics and measured baseline diagnostics."""
import argparse,csv,datetime,hashlib,json
from pathlib import Path
from build_method6_report import esc,pdf_pages
from build_purity_diagnostic_update import compile_pdf


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def measured(row):return '%.4f +/- %.4f'%(float(row['value']),float(row['error']))


def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ['output-dir','org-report','mainconf','runtime-bundle']:p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();o=a.output_dir.resolve(strict=True);stem='lambda_he4_mixing_kstar_audit_20261102'
 if (o/(stem+'.pdf')).exists():raise FileExistsError(o/(stem+'.pdf'))
 baseline=json.loads((o/'baseline_audit.json').read_text());qa=json.loads((o/'mixing_qa_summary.json').read_text());math=json.loads((o/'kstar_prf_math_audit_20261102_result.json').read_text());functions=json.loads((o/'runtime_functions_verified.json').read_text());controls=json.loads((o/'control_summary.json').read_text());rows=list(csv.DictReader((o/'baseline_points.tsv').open(),delimiter='\t'))
 assert math['allPassed'] and len(functions)==6 and qa['parameters']['mixingVzClampedEvents']==0
 labels=baseline['settings']['labels'];groups=['cent20_80','cent0_20','cent0_10','cent10_20','cent20_60']
 table=[]
 for gid in groups:
  r={x['series']:x for x in rows if x['kind']=='group' and x['group']==gid and abs(float(x['low'])-.74)<1.e-9}
  table.append([labels[gid],measured(r['raw']),measured(r['method5']),measured(r['centExposureMatched'])])
 sections=[dict(title='Finding and normalization status',paragraphs=[
  'The decreasing high-k* baseline in the restored report p.18 is confirmed. It occurs before purity fitting and remains after the adopted Poisson-L correction. The interval 0.40-0.50 GeV/c is only a reference scale for comparisons; this audit does not establish it as a physical normalization plateau.',
  'No wrong-PRF, missing factor of two, wrong boost sign, He4 charge-scaling or simple mixed-pair filling bug was found in the actual submitted 04BF implementation. This is a bounded code/numerical audit, not proof that the mixed-event reference describes the physical uncorrelated baseline. The origin of the remaining trend is not identified.',
  'The sample is the same 2021 P24iy production, SUMS 04BF2DAE47164B6F040E3A460E539289. Existing ROOT inputs, PID, Maker/configuration and previous reports are unchanged. No production job or new physics normalization is adopted.'],
  tables=[dict(columns=['Cent.','Raw at k*=0.75','Method5 L','Cent-exposure diagnostic'],rows=table)]),
 dict(title='k* definition and frozen-production verification',newpage=True,paragraphs=[
  'The submitted routine forms q=pA-pB and P=pA+pB, boosts q with -P.BoostVector(), and returns half the magnitude of its spatial vector. In the pair rest frame, pA*=-pB*, so this equals the momentum magnitude of either particle, including unequal particle masses. Fixed-target laboratory input requires no preliminary CM shift.',
  'SE and ME call the same ComputeKStar and use the same candidate four-vectors. Phi-deuteron in the shared StFemtoMaker uses this same pair routine. The separate StCandidateFemtoMaker and Legacy mapping helper implement the same PRF convention.',
  'He4 starts from primary-track rigidity momentum and receives the charge-two momentum correction once. Its energy uses mass 3.72742 GeV/c^2. The KF Lambda pairing four-vector uses configured mass 1.115683 GeV/c^2; reconstructed invariant mass is stored separately for mass fits.',
  'Six production functions, including both pair loops, pool binning, event storage and MakeLambdaEvent, are byte-identical to the current functions. The full Maker diff changes only generic non-Lambda bad-run/pileup guards; Lambda events return through their earlier dedicated branch.'],
  equations=[r'k^*=\frac{|\boldsymbol q^*|}{2}=|\boldsymbol p_A^*|=\frac{\sqrt{[s-(m_A+m_B)^2][s-(m_A-m_B)^2]}}{2\sqrt{s}}.'],
  tables=[dict(columns=['Actual ROOT5 test','Result'],rows=[['Runtime','ROOT 5.34/38, exact submitted function'],['Generated pairs','60,000: Lambda-He4, phi-d, Lambda-He3'],['Known PRF max error (Lambda-He4)','7.55e-15 GeV/c'],['Independent invariant formula max error','3.43e-12 GeV/c'],['Exchange / common-boost check','PASS'],['Float candidate storage max error','1.07e-6 GeV/c (0.00107 MeV/c)'],['Test domain','k*=0-1.5 GeV/c; common boost |beta|<=0.95']])]),
 dict(title='Mixed-event bookkeeping and data checks',newpage=True,paragraphs=[
  'Submitted mixing is bufferAll, both directions, at most20 buffered events, nine cent9 bins, ten Vz bins over198-202cm, and one event-plane bin. Every eligible candidate pair is attempted; no capped sampling or inverse-sampling correction applies in this run.',
  'The current event is mixed before it is stored. Its eventIndex increases for each input event. Shared-track checks compare eventIndex and trackIndex together. Forward current-Lambda x buffered-nucleus and reverse buffered-Lambda x current-nucleus have independent eligibility. Base channels are processed once; their helper fans out signal/sideband classes without repeating mixing.',
  'The absence of a dedicated event Vz cut made edge-bin clamping a candidate concern, but the recorded mixingVzClampedEvents counter is zero. The coarse Vz QA bin201-202.5cm crosses the mixing boundary: its contents must not be counted as events outside202cm. The current evidence does not support out-of-range clamping as the cause.',
  'The pool key has cent9/Vz but no run, field-sign or measured event-plane matching. Whether those omissions cause a bias requires event-level comparisons; a missing key is not by itself a demonstrated bug. Completely candidate-empty events are not stored; this alone does not prove a shape bias.'],
  tables=[dict(columns=['He4 channel, all Lambda mass candidates','Count'],rows=[[k,str(qa['he4Sampler'][k])] for k in ['eligible','attempted','filled','filledForward','filledReverse','skippedByCap','skippedByPairCut','skippedOverlap']])]),
 dict(title='What the measured controls exclude',newpage=True,paragraphs=[
  'The decrease remains in individual cent9 samples. For10-20%, the raw0.70-0.80 interval is0.87625 after the0.40-0.50 reference scaling. At k*=0.75 its raw point is0.88517, and its Method5 L point is0.88137. Combining different cent9 bins cannot be the sole cause.',
  'An explicitly diagnostic ME exposure match multiplies each fine-cent denominator by its own SE_norm/ME_norm before summation. The group20-80% raw point at0.75 changes only from0.88136 to0.87932. This does not recover a flat baseline and is not adopted as a physics correction.',
  'Physical He3 and the same He3 tracks under the He4 mass hypothesis are additional raw controls. For0-20% at0.75 they are0.96128 and0.95840, while selected He4 is0.90766. A trend is also present in the control population, but changing the mass hypothesis alone does not reproduce the full He4 effect. These are not purity-corrected null measurements.',
  'The L-extraction raw counts use the documented whole-mass-bin window, while species-control histograms use the exact Maker candidate gate. The maximum absolute difference of their He4 normalized curves over0.30-0.80 is0.00102, much smaller than the broad deficit. This gate difference is kept explicit.',
  'Point errors are statistical and rely on the stored/count-only assumptions. L-derived CF errors retain the existing covariance treatment. Descriptive slope fits omit event/mixed-pair covariance; no calibrated discovery significance is claimed.']),
 dict(title='Unresolved cause and the next discriminating measurements',paragraphs=[
  'A theoretical femtoscopic component near unity does not by itself make the experimental SE/ME ratio flat: the latter also tests whether the event ensembles and their conditional single-particle/pair distributions match. This audit establishes the trend, not its physical origin.',
  'Next compare event-level candidate momentum distributions and ME closure within finer multiplicity and Vz classes, run/field classes, and measured event-plane classes. Test eventwise conservation/flow baselines and detector/TOF acceptance dependence using a controlled same-candidate remix. Hold PID, Lambda mass selection and event cuts fixed while changing one mixing condition at a time.',
  'The stored projected histograms cannot recreate arbitrary event-level remixes. A candidate/event record that includes both Lambda and nuclear four-vectors and event-class identifiers is needed for that comparison. No full reprocessing is launched by this audit. Choosing a different normalization interval alone is not a remedy for a sloping baseline.'])]
 figures=[('baseline_group_L.pdf','Raw and L-derived high-k* baselines'),('baseline_fine_L.pdf','Fine-centrality check: raw and Method3 L'),('baseline_cent_exposure.pdf','Centrality-composition diagnostic; not an adopted correction'),('baseline_species_control.pdf','Species/mass-hypothesis control: raw, without purity correction')]
 now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9)))
 text=[r'''\documentclass[10pt,a4paper]{article}
\usepackage[utf8]{inputenc}\usepackage[T1]{fontenc}\usepackage{geometry,graphicx,amsmath,longtable,booktabs,hyperref}
\geometry{margin=1.7cm}\setlength{\parindent}{0pt}\setlength{\parskip}{.4em}\emergencystretch=3em
\begin{document}\begin{center}{\Large $\Lambda$--$^4$He: high-$k^*$ baseline and mixing audit}\end{center}
''',esc(now.strftime('%Y-%m-%d %H:%M JST'))+'\n\n']
 for sec in sections:
  if sec.get('newpage'):text.append('\\clearpage\n')
  text.append('\\section*{'+esc(sec['title'])+'}\n')
  for line in sec['paragraphs']:text.append(esc(line)+'\n\n')
  for eq in sec.get('equations',[]):text.append('\\['+eq+'\\]\n')
  for table in sec.get('tables',[]):
   text.append(r'{\small\setlength{\tabcolsep}{3pt}\begin{longtable}{'+'l'*len(table['columns'])+r'}\toprule'+'\n')
   text.append(' & '.join(esc(x) for x in table['columns'])+r'\\\midrule\endhead'+'\n')
   for row in table['rows']:text.append(' & '.join(esc(x) for x in row)+r'\\'+'\n')
   text.append(r'\bottomrule\end{longtable}}'+'\n')
 text.append('\\clearpage\n');figs=[]
 for name,title in figures:
  path=o/'plots'/name;assert pdf_pages(path)==1;figs.append(dict(path=str(path),sha256=sha(path),pages=1))
  text.append('\\section*{'+esc(title)+'}\n'+esc('0.40-0.50 GeV/c is a reference scale only; no physical plateau is established.')+'\n\\begin{center}\\includegraphics[width=\\textwidth,height=.8\\textheight,keepaspectratio]{'+str(path)+'}\\end{center}\\clearpage\n')
 text.append('\\end{document}\n');tex=o/(stem+'.tex');tex.write_text(''.join(text));report=compile_pdf(tex,o)
 sources=[o/'baseline_audit.json',o/'mixing_qa_summary.json',o/'control_summary.json',o/'runtime_functions_verified.json',o/'kstar_prf_math_audit_20261102_result.json',o/'baseline_settings.yaml',o/'plot_settings.yaml',a.org_report,a.mainconf,a.runtime_bundle]
 record=dict(phase='built',report=report,figures=figs,sourceRecords={str(x.resolve()):dict(sha256=sha(x),bytes=x.stat().st_size) for x in sources},sourceCode=dict(path=str(Path(__file__).resolve()),sha256=sha(__file__)),normalizationStatus='reference scale only; plateau not established',causeStatus='unresolved; no tested kstar or elementary ME bookkeeping bug found')
 save(o/'audit_report_metadata.json',record)
 (o/'audit_report.md').write_text('# Λ–4He高k*ベースラインとME/k*調査\n\n**高k*の低下は確認できたが、その原因は未特定。k*の静止系変換・電荷補正・単純なME fillの誤りは確認されなかった。0.4–0.5を物理的plateauとして認定しない。**\n\n'+ '- [図付き調査PDF]('+str(Path(report['path']))+')\n\n'+'\n\n'.join('## '+sec['title']+'\n\n'+'\n\n'.join(sec['paragraphs']) for sec in sections)+'\n')
 print(json.dumps(report))

if __name__=='__main__':main()
