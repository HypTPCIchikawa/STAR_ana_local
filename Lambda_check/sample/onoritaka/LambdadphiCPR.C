double sig[200][2][3][4][4], pol[200][2][3][4][4];
int ck, ce, pt, ra;
double fitFunction0(double *x, double *par) {
	return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.0)) + par[3] + par[4]*x[0] + par[5]*x[0]*x[0];
}
double fitFunction1(double *x, double *par) {
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}
double Purity(double *x, double *par) {
	return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.0))/(par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.0)) + par[3] +par[4]*x[0]);
}
double fitFunction2(double *x, double *par) {
	int i0 = (x[0]-1.07)/(0.1/200);
	if(i0 < 0){ cout << "error0" << endl; i0 = 0; }
	if(i0 > 199){ cout << "error1" << endl; i0 = 199; }
	double gauss = sig[i0][ck][ce][pt][ra];
	double pol0 = pol[i0][ck][ce][pt][ra];
	double f0 = 0.0;
	if (gauss+pol0 != 0){
		f0 = gauss/(gauss+pol0);
	} else {
		f0 = 0.0;
	}
	//return par[5]*f0+(par[6]*(x[0]-1.1157)+par[7])*(1.0-f0);
	return par[5]*f0+(par[6]*x[0]+par[7])*(1.0-f0);
}
double sinfit(double *x, double *par) {
	return par[0]+par[1]*sin(x[0])+par[2]*sin(2*x[0])+par[3]*sin(3*x[0])+par[4]*sin(4*x[0])+par[5]*sin(5*x[0]);
}
void LambdadphiCPR(){
	TCanvas *tmp = new TCanvas("Masspos","Masspos",50,50,1000,900);
	tmp->Divide(8,6);
	TCanvas *tpp = new TCanvas("Puritypos","Puritypos",50,50,1000,900);
	tpp->Divide(8,6);
	TCanvas *tvp = new TCanvas("v1pos","v1pos",50,50,1000,900);
	tvp->Divide(8,6);
	TCanvas *tgp = new TCanvas("gpolpos","gpolpos",50,50,1000,900);
	tgp->Divide(8,6);
	TCanvas *tmn = new TCanvas("Massneg","Massneg",50,50,1000,900);
	tmn->Divide(8,6);
	TCanvas *tpn = new TCanvas("Purityneg","Purityneg",50,50,1000,900);
	tpn->Divide(8,6);
	TCanvas *tvn = new TCanvas("v1neg","v1neg",50,50,1000,900);
	tvn->Divide(8,6);
	TCanvas *tgn = new TCanvas("gpolneg","gpolneg",50,50,1000,900);
	tgn->Divide(8,6);
	TCanvas *tvpn = new TCanvas("v1pn","v1pn",50,50,1000,900);
	tgn->Divide(16,6);
	
	TFile *tf = TFile::Open(Form("Out6.root"));
	cout << "open : " << Form ("Out6.root") <<endl;

	TH2D *MasAL[10][4][4];
	TProfile2D *v1AL[10][4][4];
	TProfile2D *gpolAL[10][4][4];
	TH2D *MasALC[3][4][4];
	TProfile2D *v1ALC[4][4][4];
	TProfile2D *gpolALC[4][4][4];
	TH1D *ll[2][3][4][4];
	TF1 *ff[5][2][3][4][4];
	TH1D *gg[2][3][4][4];
	TH1D *GpolALfit[2][3][4][4];
	TH1D *v1ALfit[2][3][4][4];

	float a[2][3][4][4], b[2][3][4][4], c[2][3][4][4], d[2][3][4][4], e[2][3][4][4];
	double tot[50][2][3][4][4], xx[50][2][3][4][4];
	float j[2][3][4][4], p[2][3][4][4], l[2][3][4][4], mm[2][3][4][4], nn[2][3][4][4];

	
	TH1D *hist[2][3][4][4];
	TH1D *histv[2][3][4][4];
	float pp[2][3][4][4], qq[2][3][4][4], rr[2][3][4][4], err[2][3][4][4],aerr[2][3][4][4],berr[2][3][4][4];
	float ppv[2][3][4][4], qqv[2][3][4][4], rrv[2][3][4][4], errv[2][3][4][4],aerrv[2][3][4][4],berrv[2][3][4][4];
	float Mt = 0;
	float mt = 0;

	int sigma = 3;//can choose 2 or 3
	if(sigma==3){
		Mt = 1.1215;
		mt = 1.1100;
	} else if(sigma== 2){
		Mt = 1.1195;
		mt = 1.1120;
	} else if(sigma== 1){
		Mt = 1.13;
		mt = 1.105;
	}

	int Sbin_min = (mt-1.07)/0.0005;
	int Sbin_max =(Mt-1.07)/0.0005;

	
	for( int im=0; im<10; im++ ){
		for( int iop=0; iop<4; iop++ ){
			for( int ior=0; ior<4; ior++ ){
				MasAL[im][iop][ior] = (TH2D*)tf->Get(Form("hLambdaMassCPR_%d_%d_%d",im,iop,ior));
				v1AL[im][iop][ior] = (TProfile2D*)tf->Get(Form("v1dphiALCPR_%d_%d_%d",im,iop,ior));
				gpolAL[im][iop][ior] = (TProfile2D*)tf->Get(Form("GlobalPoldphiALCPR_%d_%d_%d",im,iop,ior));
				MasAL[im][iop][ior]->RebinY(25);
				v1AL[im][iop][ior]->RebinY(25);
				gpolAL[im][iop][ior]->RebinY(25);

				if(im == 0){
					MasALC[0][iop][ior] = MasAL[im][iop][ior];
					v1ALC[0][iop][ior] = v1AL[im][iop][ior];
                                        gpolALC[0][iop][ior] =  gpolAL[im][iop][ior];
                                }else if(im == 1){
                                        MasALC[0][iop][ior]->Add(MasAL[im][iop][ior]);
                                        v1ALC[0][iop][ior]->Add(v1AL[im][iop][ior]);
                                        gpolALC[0][iop][ior]->Add(gpolAL[im][iop][ior]);
                                }else if(im == 2){
					MasALC[1][iop][ior] = MasAL[im][iop][ior];
					v1ALC[1][iop][ior] = v1AL[im][iop][ior];
                                        gpolALC[1][iop][ior] =  gpolAL[im][iop][ior];
                                }else if(im == 3){
                                        MasALC[1][iop][ior]->Add(MasAL[im][iop][ior]);
                                        v1ALC[1][iop][ior]->Add(v1AL[im][iop][ior]);
                                        gpolALC[1][iop][ior]->Add(gpolAL[im][iop][ior]);
                                }else if(im == 4){
                                        MasALC[1][iop][ior]->Add(MasAL[im][iop][ior]);
                                        v1ALC[1][iop][ior]->Add(v1AL[im][iop][ior]);
                                	gpolALC[1][iop][ior]->Add(gpolAL[im][iop][ior]);
                                }else if(im == 5 ){
					MasALC[2][iop][ior] = MasAL[im][iop][ior];
					v1ALC[2][iop][ior] = v1AL[im][iop][ior];
                                        gpolALC[2][iop][ior] =  gpolAL[im][iop][ior];
                                }else if(im>5){
                                        MasALC[2][iop][ior]->Add(MasAL[im][iop][ior]);
                                        v1ALC[2][iop][ior]->Add(v1AL[im][iop][ior]);
                                        gpolALC[2][iop][ior]->Add(gpolAL[im][iop][ior]);
                                }
			}
		}
	}
	for(int kk=0; kk<2; kk++){
		for( int im=0; im<3; im++ ){
			for( int iop=0; iop<4; iop++ ){
				for( int ior=0; ior<4; ior++ ){
					ll[kk][im][iop][ior] = new TH1D(Form("LambdaMass_%d_%d_%d_%d",kk,im,iop,ior), Form("LambdaMass_%d_%d_%d_%d",kk,im,iop,ior), 200, 1.07, 1.17);
					for(int n=0; n<200; n++){
						float y = MasALC[im][iop][ior]->GetBinContent(n+1,kk+1);
						ll[kk][im][iop][ior]->SetBinContent(n+1,y);
						y = 0;
					}
					if(ll[kk][im][iop][ior]){
						if(kk ==1){
							tmp->cd(16*im+4*ior+iop+1);
							ll[kk][im][iop][ior]->Draw();
						}else{
							tmn->cd(16*im+4*ior+iop+1);
							ll[kk][im][iop][ior]->Draw();
						}
						ff[0][kk][im][iop][ior] = new TF1(Form("ff_0_%d_%d_%d_%d",kk,im,iop,ior),fitFunction1,1.13,1.17,3);
						ff[1][kk][im][iop][ior] = new TF1(Form("ff_1_%d_%d_%d_%d",kk,im,iop,ior),fitFunction0,1.09,1.17,6);
						ff[2][kk][im][iop][ior] = new TF1(Form("ff_2_%d_%d_%d_%d",kk,im,iop,ior),fitFunction1,1.09,1.17,3);
						ff[0][kk][im][iop][ior]->SetParameter(0,1.0);
						ff[0][kk][im][iop][ior]->SetParameter(1,1.0);
						ff[0][kk][im][iop][ior]->SetParameter(2,0.0);
						ll[kk][im][iop][ior]->Fit(Form("ff_0_%d_%d_%d_%d",kk,im,iop,ior),"N","",1.13,1.17);
						a[kk][im][iop][ior] = ff[0][kk][im][iop][ior]->GetParameter(0);
				   		b[kk][im][iop][ior] = ff[0][kk][im][iop][ior]->GetParameter(1);
   						e[kk][im][iop][ior] = ff[0][kk][im][iop][ior]->GetParameter(2);
   						c[kk][im][iop][ior] = ff[0][kk][im][iop][ior]->Eval(1.1157);
   						d[kk][im][iop][ior] = ll[kk][im][iop][ior] ->GetBinContent(91);

						ff[1][kk][im][iop][ior]->SetParameter(0,d[kk][im][iop][ior]-c[kk][im][iop][ior]);
						ff[1][kk][im][iop][ior]->FixParameter(1,1.1155);
						ff[1][kk][im][iop][ior]->SetParameter(2,0.002);
				   		ff[1][kk][im][iop][ior]->SetParameter(3,a[kk][im][iop][ior]);
  						ff[1][kk][im][iop][ior]->SetParameter(4,b[kk][im][iop][ior]);
  						ff[1][kk][im][iop][ior]->SetParameter(5,e[kk][im][iop][ior]);
						ll[kk][im][iop][ior]->Fit(Form("ff_1_%d_%d_%d_%d",kk,im,iop,ior),"","",1.09,1.17);

						ff[2][kk][im][iop][ior]->FixParameter(0,ff[1][kk][im][iop][ior]->GetParameter(3));
						ff[2][kk][im][iop][ior]->FixParameter(1,ff[1][kk][im][iop][ior]->GetParameter(4));
						ff[2][kk][im][iop][ior]->FixParameter(2,ff[1][kk][im][iop][ior]->GetParameter(5));
						ff[2][kk][im][iop][ior]->SetLineColor(4);
						ff[2][kk][im][iop][ior]->Draw("same");
					}
					ck = kk;
					ce = im;
					pt = iop;
					ra = ior;
					gg[kk][im][iop][ior] = new TH1D(Form("purity_%d_%d_%d_%d",kk,im,iop,ior),"",200,1.07,1.17);
					for(int ii=0; ii<200; ii++){
						if(ii > 70 && ii <130){ 
							tot[ii][kk][im][iop][ior] = ll[kk][im][iop][ior]->GetBinContent(ii+1);
							xx[ii][kk][im][iop][ior] = ll[kk][im][iop][ior]->GetBinCenter(ii+1);
							pol[ii][kk][im][iop][ior] = ff[2][kk][im][iop][ior]->Eval(xx[ii][kk][im][iop][ior]);
							sig[ii][kk][im][iop][ior] = tot[ii][kk][im][iop][ior]-pol[ii][kk][im][iop][ior];
							if(sig[ii][kk][im][iop][ior]<0) sig[ii][kk][im][iop][ior]=0;
							if(tot[ii][kk][im][iop][ior]>0) gg[kk][im][iop][ior]->SetBinContent(ii+1,sig[ii][kk][im][iop][ior]/tot[ii][kk][im][iop][ior]);
						}else{
							gg[kk][im][iop][ior]->SetBinContent(ii+1,0);
						}
					}
					gg[kk][im][iop][ior]->SetMinimum(0.0);
					gg[kk][im][iop][ior]->SetMaximum(1.0);
					if(kk ==1){
						tpp->cd(16*im+4*ior+iop+1);
						gg[kk][im][iop][ior]->Draw();
					}else{
						tpn->cd(16*im+4*ior+iop+1);
						gg[kk][im][iop][ior]->Draw();
					}
					float N = 0;
					float aa = 0;
					for (int jj = Sbin_min; jj<=Sbin_max; jj++) {
						float xi = ll[kk][im][iop][ior]->GetBinCenter(jj);
						N += ff[2][kk][im][iop][ior]->Eval(xi);
						float cc = ll[kk][im][iop][ior]->GetBinContent(jj);
						aa += cc;
					}

					float S = aa-N;
					//float SN = S/aa
					//float E = (S/aa)*S;
					float Ssqrt = S/sqrt(aa);
					cout << "K=" << kk << ",Cent="<< im << ",pT=" << iop << ",rap=" << ior << endl;
					//cout << "signal(all)" << aa << endl;
					cout << "SN=" << S/N << endl;
					//cout << "N/a=" << N/aa << endl;
					//cout << "S/a=" << S/aa << endl;
					cout << "統計的有意性" << Ssqrt << endl;

					j[kk][im][iop][ior] = ff[1][kk][im][iop][ior]->GetParameter(0);
					p[kk][im][iop][ior] = ff[1][kk][im][iop][ior]->GetParameter(1);
					l[kk][im][iop][ior] = ff[1][kk][im][iop][ior]->GetParameter(2);
					mm[kk][im][iop][ior] = ff[1][kk][im][iop][ior]->GetParameter(3);
					nn[kk][im][iop][ior] = ff[1][kk][im][iop][ior]->GetParameter(4);

					GpolALfit[kk][im][iop][ior] = new TH1D(Form("GlobalPolALfit_%d_%d_%d_%d",kk,im,iop,ior), Form("GlobalPolALfit_%d_%d_%d_%d",kk,im,iop,ior), 100, 1.09, 1.14);
					for(int jh=0; jh<100; jh++){
								double yg = gpolALC[im][iop][ior]->GetBinContent(jh+1,kk+1);
								double yger = gpolALC[im][iop][ior]->GetBinError(jh+1,kk+1);
								GpolALfit[kk][im][iop][ior]->SetBinContent(jh+1,yg);
								GpolALfit[kk][im][iop][ior]->SetBinError(jh+1,yger);
					}
					if(GpolALfit[kk][im][iop][ior]){
						if(kk ==1){
							tgp->cd(16*im+4*ior+iop+1);
							GpolALfit[kk][im][iop][ior]->Draw();
						}else{
							tgn->cd(16*im+4*ior+iop+1);
							GpolALfit[kk][im][iop][ior]->Draw();
						}
						ff[3][kk][im][iop][ior] = new TF1(Form("ff_3_%d_%d_%d_%d",kk,im,iop,ior),fitFunction2,1.09,1.14,8);
						ff[3][kk][im][iop][ior]->FixParameter(0,j[kk][im][iop][ior]);
						ff[3][kk][im][iop][ior]->FixParameter(1,p[kk][im][iop][ior]);
						ff[3][kk][im][iop][ior]->FixParameter(2,l[kk][im][iop][ior]);
						ff[3][kk][im][iop][ior]->FixParameter(3,mm[kk][im][iop][ior]);
						ff[3][kk][im][iop][ior]->FixParameter(4,nn[kk][im][iop][ior]);
						ff[3][kk][im][iop][ior]->SetParameter(5,1.0);
						ff[3][kk][im][iop][ior]->SetParameter(6,1.0);
						ff[3][kk][im][iop][ior]->SetParameter(7,1.0);
						GpolALfit[kk][im][iop][ior]->Fit(Form("ff_3_%d_%d_%d_%d",kk,im,iop,ior),"nq","",1.09,1.14);
						hist[kk][im][iop][ior] = new TH1D(Form("hist_%d_%d_%d_%d",kk,im,iop,ior),"",200,1.07,1.17);
						for(int ib=0; ib<200; ib++){
		 					float xfg = ll[kk][im][iop][ior]->GetBinCenter(ib+1);
		 					float yfg = ff[3][kk][im][iop][ior]->Eval(xfg);
	             					hist[kk][im][iop][ior]->SetBinContent(ib+1,yfg);
						}
						hist[kk][im][iop][ior]->SetLineColor(2);
						hist[kk][im][iop][ior]->Draw("same");
//						ff[3][kk][im][iop][ior]->SetLineColor(2);
//						ff[3][kk][im][iop][ior]->Draw("same");
						//pp[kk][im][iop][ior] = ff[3][kk][im][iop][ior]->GetParameter(5);
						//qq[kk][im][iop][ior] = ff[3][kk][im][iop][ior]->GetParameter(6);
						//rr[kk][im][iop][ior] = ff[3][kk][im][iop][ior]->GetParameter(7);
						//cout << "A=" << pp[kk][im][iop][ior] << endl;
						//err[kk][im][iop][ior] = ff[3][kk][im][iop][ior]->GetParError(5);
						//cout << "err=" << err[kk][im][iop][ior] << endl;
						//aerr[kk] = ff[3][kk]->GetParError(6);
						//berr[kk] = ff[3][kk]->GetParError(7);
					}
					v1ALfit[kk][im][iop][ior] = new TH1D(Form("v1ALfit_%d_%d_%d_%d",kk,im,iop,ior), Form("v1ALfit_%d_%d_%d_%d",kk,im,iop,ior), 100, 1.09, 1.14);
					for(int jh=0; jh<100; jh++){
								double xv = v1ALC[im][iop][ior]->GetBinContent(jh+1,kk+1);
								double xver = v1ALC[im][iop][ior]->GetBinError(jh+1,kk+1);
								v1ALfit[kk][im][iop][ior]->SetBinContent(jh+1,xv);
								v1ALfit[kk][im][iop][ior]->SetBinError(jh+1,xver);
					}
					if(v1ALfit[kk][im][iop][ior]){
						v1ALfit[kk][im][iop][ior]->SetMaximum(1.0);
						v1ALfit[kk][im][iop][ior]->SetMinimum(-0.2);
						if(kk ==1){
							tvp->cd(16*im+4*ior+iop+1);
							v1ALfit[kk][im][iop][ior]->Draw();
						}else{
							tvn->cd(16*im+4*ior+iop+1);
							v1ALfit[kk][im][iop][ior]->Draw();
						}
						ff[4][kk][im][iop][ior] = new TF1(Form("ff_4_%d_%d_%d_%d",kk,im,iop,ior),fitFunction2,1.09,1.14,8);
						ff[4][kk][im][iop][ior]->FixParameter(0,j[kk][im][iop][ior]);
						ff[4][kk][im][iop][ior]->FixParameter(1,p[kk][im][iop][ior]);
						ff[4][kk][im][iop][ior]->FixParameter(2,l[kk][im][iop][ior]);
						ff[4][kk][im][iop][ior]->FixParameter(3,mm[kk][im][iop][ior]);
						ff[4][kk][im][iop][ior]->FixParameter(4,nn[kk][im][iop][ior]);
						ff[4][kk][im][iop][ior]->SetParameter(5,1.0);
						ff[4][kk][im][iop][ior]->SetParameter(6,1.0);
						ff[4][kk][im][iop][ior]->SetParameter(7,1.0);
						v1ALfit[kk][im][iop][ior]->Fit(Form("ff_4_%d_%d_%d_%d",kk,im,iop,ior),"nq","",1.09,1.14);
						histv[kk][im][iop][ior] = new TH1D(Form("histv_%d_%d_%d_%d",kk,im,iop,ior),"",200,1.07,1.17);
						histv[kk][im][iop][ior]->SetLineColor(2);
						for(int ib=0; ib<200; ib++){
		 					float x1 = ll[kk][im][iop][ior]->GetBinCenter(ib+1);
		 					float y1 = ff[4][kk][im][iop][ior]->Eval(x1);
                					histv[kk][im][iop][ior]->SetBinContent(ib+1,y1);
						}
						histv[kk][im][iop][ior]->Draw("same");
						ppv[kk][im][iop][ior] = ff[4][kk][im][iop][ior]->GetParameter(5);
						//qq[kk][im][iop][ior] = ff[4][kk][im][iop][ior]->GetParameter(6);
						//rr[kk][im][iop][ior] = ff[4][kk][im][iop][ior]->GetParameter(7);
						//cout << "A=" << ppv[kk][im][iop][ior] << endl;
						//errv[kk][im][iop][ior] = ff[4][kk][im][iop][ior]->GetParError(5);
						//cout << "err=" << errv[kk][im][iop][ior] << endl;
						//aerr[kk] = ff[3][kk]->GetParError(6);
						//berr[kk] = ff[3][kk]->GetParError(7);
					}
					if(v1ALfit[kk][im][iop][ior]){
						tvp->cd(32*im+8*ior+2*iop+kk+1);
						v1ALfit[kk][im][iop][ior]->SetMaximum(0.2);
						v1ALfit[kk][im][iop][ior]->SetMinimum(-0.2);
						v1ALfit[kk][im][iop][ior]->Draw();
					}
						
				}
			}
		}
	}
}
