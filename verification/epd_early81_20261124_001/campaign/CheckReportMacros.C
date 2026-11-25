void CheckReportMacros(const char* folder) {
  TString base(folder);
  TInterpreter::EErrorCode error=TInterpreter::kNoError;
  gROOT->ProcessLine(".L "+base+"/VerifyEarlyPhysicsPilot.C+",&error);
  if(error!=TInterpreter::kNoError) {gSystem->Exit(1);return;}
  gROOT->ProcessLine(".L "+base+"/CombineEarlyCF.C+",&error);
  if(error!=TInterpreter::kNoError) {gSystem->Exit(1);return;}
  std::cout<<"EARLY_REPORT_MACROS_COMPILE_PASS"<<std::endl;
}
