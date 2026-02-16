#include <TSystem>
#include "TStopwatch.h"

class StMaker;
class StChain;
class StPicoDstMaker;
class StBES2QaMaker;
StChain *chain;

int runMode = 1;
        // 0= pass0: make gain correction factors of BBC-ADC 
		// 1= pass1: make recentering parameters
		// 2= pass2: make flattnining parameters
		// 3= analysis mode
		

void readPicoDst(const Char_t *inputFile="test.list", const Char_t *outputFile="test.root")
{
	TStopwatch timer;
	timer.Start();

    Int_t nEvents = 10000000;
	//Int_t nEvents = 500;	

	TString sInput(outputFile);
	TString RunNum = sInput(0,8);

	//Load all the System libraries
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("StPicoEvent"); // to be called first
	gSystem->Load("StPicoDstMaker");
    gSystem->Load("StBES2QaMaker");
	gSystem->Load("QVMaker");
	gSystem->Load("StRefMultCorr");
	gSystem->Load("SgPEvent");
	gSystem->Load("StEpdUtil");


	chain = new StChain();
	StPicoDstMaker    *picoMaker = new StPicoDstMaker(StPicoDstMaker::IoRead,inputFile,"picoDst");
	picoMaker->SetStatus("*",0);
	picoMaker->SetStatus("Event",1);
	picoMaker->SetStatus("Track",1);
	picoMaker->SetStatus("BTofHit",1);
	picoMaker->SetStatus("BTofPidTraits",1);
	picoMaker->SetStatus("BbcHit",1);
	picoMaker->SetStatus("EpdHit",1);
	picoMaker->SetStatus("MtdHit",1);
	picoMaker->SetStatus("BTowHit",1);
	picoMaker->SetStatus("ETofPidTraits",1);
	//picoReader->SetStatus("TrackCovMatrix",1);


    StBES2QaMaker *anaMaker  = new StBES2QaMaker("ana",picoMaker,outputFile);
	anaMaker->SetRunNumber(RunNum.Atoi());
	anaMaker->SetCalibrationMode(runMode);
		// 0= gain correction/QA
		// 1= make recentering para.
		// 2= make flattening para.
		// 3= analysis mode
	anaMaker->SetData(14);
		//  0=run10
		//  1=run11
		//  2=run14
		//  3=run16
		//  4=run18
		//  5=run19 19GeV
		//  6=run19 14GeV
		//  7=run19 7GeV
		//  8=run19 9GeV
		//  9=run19 FXT
		// 10=run19 FXT 31 GeV 
		// 11=run19 200GeV
		// 12=run20 11.5GeV
		// 13=run20 9.2GeV
		// 14=run20 7.3GeV
	anaMaker->SetZVertexMax(70.0);
	//anaMaker->SetZVertexMax(50.0); // 200GeV
	anaMaker->SetVertexRMax(2.0);
	anaMaker->SetVpdVzDif  (20.0);
	//anaMaker->SetVpdVzDif  (5.0); // 200GeV

	// Copied from Prithwish's code
	//Accroding to the RunLog browser Au+Au 27 GeV (Run 18)
	//the first run is: 19130060
	//the last run is:  19268002
	//anaMaker->SetRunRange(19139960,19268002);
	
	// fast offline run19 19GeV
	//anaMaker->SetRunRange(20056020,20056050); // 056
	//anaMaker->SetRunRange(20057003,20057050); // 057
	//anaMaker->SetRunRange(20057003,20100000); // for production_19GeV_2019
	//anaMaker->SetRunRange(20094000,20170000); // for production_14GeV_2019
	//anaMaker->SetRunRange(20154047,20180000); // for production_7GeV_2019
	//anaMaker->SetRunRange(20179000,20185000); // for production_4p59GeV_FXT_2019
	//anaMaker->SetRunRange(20179000,20195000); // for production_9GeV_2019
	//anaMaker->SetRunRange(20189000,20195000); // for production_31GeV_FXT_2019
	//anaMaker->SetRunRange(20191000,20199000); // for production_200GeV_2019
	
	// run20 fast offline
	//anaMaker->SetRunRange(20344000,20380000); // until 12/31/2019 for production_200GeV_2019
	//anaMaker->SetRunRange(21000000,21040000); // from 1/1/2020 
	//anaMaker->SetRunRange(21028000,21040000); // FXT 2020 
	//anaMaker->SetRunRange(21037000,21057000); // 2020 collider for 9p2
	//anaMaker->SetRunRange(21044000,21046000); // run20 FXT 5.75 
	//anaMaker->SetRunRange(21040000,21090000); // run20 11.5 
	//anaMaker->SetRunRange(21050000,21070000); // run20 9p2b 
	//anaMaker->SetRunRange(21050000,21080000); // run20 9p2b 
	anaMaker->SetRunRange(21035000,21036020); // run20 9p2b 

	// These are for Fixed-Target, comment-out for collider mode
	anaMaker->SetFixedMode(1); // 0=collider,  1=fixed-target for vertex cuts
	anaMaker->SetVertexRMax(2.0); // relative to (-0.4,-2.0) for FXT


	if( chain->Init()==kStErr ){ 
		cout<<"chain->Init();"<<endl;
		return;
	}

	int total = picoMaker->chain()->GetEntries();
    cout << " Total entries = " << total << endl;
    if(nEvents>total) nEvents = total;

	for (Int_t i=0; i<nEvents; i++){

	  if(i%1000==0)
	  //if(i%(nEvents/20)==0)
		cout << "Working on eventNumber " << i << endl;
		
	  chain->Clear();
	  int iret = chain->Make(i);
		
	  if (iret) { cout << "Bad return code!" << iret << endl; break;}

	  total++;
	}
	
	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;
	timer.Stop();
	std::cout << "RealTime: " << timer.RealTime() << " CpuTime: " << timer.CpuTime() << std::endl;

	
	delete anaMaker;
	delete picoMaker;
	delete chain;
	
}
