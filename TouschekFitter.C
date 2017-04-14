#include <math.h>
#include <TH1F.h>
#include <TMatrix.h>
#include <TTree.h>

#include "fitter.h"

/*
 *     A description of what this script does
 */



const int nSubRunLER=15;
bool sim=false;

//the Unix time when each subrun starts
UInt_t subrunStartLER[] = {1463469869,//81
			   1463470115,//65
			   1463470361,//51
			   1463470690,//38
			   1463471041,//32

			   1463471443,//95
			   1463471688,//72
			   1463472068,//62
			   1463472308,//58
			   1463472578,//51

			   1463473067,//147
			   1463473412,//146
			   1463474064,//145
			   1463473852,//144
			   1463473643 //141
};

//the Unix time when each subrun ends
UInt_t subrunEndLER[] = {1463470049,//81
			 1463470295,//65
			 1463470541,//51
			 1463470870,//38
			 1463471221,//32

			 1463471623,//95
			 1463471868,//72
			 1463472248,//62
			 1463472488,//58
			 1463472758,//51

			 1463473187,//147
			 1463473532,//146
			 1463474184,//145
			 1463473972,//144
			 1463473817 //141
};





void TouschekFitter(int ch=0, bool draw=true){

  //load the ntuple
  TFile *f;
  if(!sim) f= new TFile("BEAST_run300789.root");    //real data
  else f = new TFile("mc_300789.root");             //simulated data
  TTree *tout = (TTree*)f->Get("tout");

  
  UInt_t ts;
  vector<double> *he3;
  vector<double> *Current;
  vector<double> *Pressure;
  vector<double> *BeamSize;
  vector<int> *injFlag;
  vector<double> *nBunch;
  vector<double> *Zeff;

  TString Ring = "LER";

  tout->SetBranchAddress("ts", &ts);
  tout->SetBranchAddress("HE3_rate", &he3);
  tout->SetBranchAddress("SKB_"+Ring+"_current", &Current);
  tout->SetBranchAddress("SKB_"+Ring+"_pressure_average", &Pressure);
  tout->SetBranchAddress("SKB_"+Ring+"_correctedBeamSize_xray_Y", &BeamSize);
  tout->SetBranchAddress("SKB_"+Ring+"_injectionNumberOfBunches", &nBunch);
  tout->SetBranchAddress("SKB_"+Ring+"_injectionFlag_safe", &injFlag);
  tout->SetBranchAddress("SKB_LER_Zeff_D02", &Zeff);

  //reserve variables
  double he3Mean[nSubRunLER][4] = {{0}};  //the mean rate for each of the four tubes during each subrun
  double he3n[nSubRunLER][4] = {{0}};     //for calculation of RMS
  double M2[nSubRunLER][4] = {{0}};       //for calculation of RMS

  double Cur[nSubRunLER] = {0};           //Beam Current during each subrun
  double BS[nSubRunLER] = {0};            //Beam size during each subrun
  double Pre[nSubRunLER] = {0};           //Pressure during each subrun
  double Zef[nSubRunLER] = {0};           //Z_eff during each subrun
  double nBun[nSubRunLER] = {0};          //number of bunches during each subrun
  int n[nSubRunLER] = {0};

  //loop over the ntuple
  int nentries = tout->GetEntries();
  for(int i=0; i<nentries; i++){
    tout->GetEntry(i);

    //data quality cuts
    if(he3->size()!=4) continue;
    if(Pressure->size()!=1) continue;
    if(Current->size()!=1) continue;
    if(BeamSize->size()!=1) continue;
    if(injFlag->size()!=1) continue;

    if(Pressure->at(0)!=Pressure->at(0))continue;
    if(Current->at(0)!=Current->at(0))continue;
    if(BeamSize->at(0)!=BeamSize->at(0))continue;
    
    if(BeamSize->at(0)>200||BeamSize->at(0)<15) continue;

    //ignore data during beam injection
    if(injFlag->at(0)!=0) continue;

    //add data from each subrun to the arrays
    for(int j=0; j<nSubRunLER; j++){
      if(ts<subrunEndLER[j]&&ts>subrunStartLER[j]){
	Cur[j] += Current->at(0);
	Pre[j] += Pressure->at(0);
	BS[j] += BeamSize->at(0);
	Zef[j] += Zeff->at(0);
	nBun[j] += nBunch->at(0);

	//calculate mean and RMS using Welford's algorithm.
	for(int k=0; k<4; k++){
	  he3n[j][k]++;
	  double delta = he3->at(k) - he3Mean[j][k];
	  he3Mean[j][k] += delta/he3n[j][k];
	  double delta2=he3->at(k) - he3Mean[j][k];
	  M2[j][k] +=delta*delta2;

	}
	n[j]++;
      }
    }
  }

  /*
   *  Now the data will be collected into TMatrixes
   */
  
  //one matrix for each he3tube channel
  TMatrixD ydata0(nSubRunLER,1);
  TMatrixD ydata1(nSubRunLER,1);
  TMatrixD ydata2(nSubRunLER,1);
  TMatrixD ydata3(nSubRunLER,1);

  //one matrix for the accelerator variables
  TMatrixD X(nSubRunLER,2);

  double error[nSubRunLER][4] = {{0}};

  for(int i=0; i<nSubRunLER; i++){

    
    X[i][0] = (Cur[i]/n[i])*(Pre[i]/n[i])*(Zef[i]/n[i])**2;  //The beam gas component: IxPxZ^2
    X[i][1] = (Cur[i]/n[i])**2/(BS[i]/n[i]*nBun[i]/n[i]);    //The Touschek component: I^2/(sigma*nBunches)

    //copy the he3tube data into the TMatrixes
    ydata0[i][0] = he3Mean[i][0];
    ydata1[i][0] = he3Mean[i][1];
    ydata2[i][0] = he3Mean[i][2];
    ydata3[i][0] = he3Mean[i][3];


    for(int j=0; j<4; j++){
      error[i][j] = sqrt(M2[i][j]/(he3n[i][j]-1));
    }


  }
  
  //use the fitting function in fitter.h to solve the least squares fit
  TMatrixD soln0 = fitter(ydata0, X);
  TMatrixD appliedSoln0 = MatrixMultiply(X, soln0);

  TMatrixD soln1 = fitter(ydata1, X);
  TMatrixD appliedSoln1 = MatrixMultiply(X, soln1);

  TMatrixD soln2 = fitter(ydata2, X);
  TMatrixD appliedSoln2 = MatrixMultiply(X, soln2);

  TMatrixD soln3 = fitter(ydata3, X);
  TMatrixD appliedSoln3 = MatrixMultiply(X, soln3);


  /*
   * Display visualizations of the fits.
   */

  
  //choose which of the 4 channels to display
  TMatrixD displayData=ydata0;
  TMatrixD displaySoln=soln0;
  if(ch==1){
    displayData=ydata1;
    displaySoln=soln1;
  }else if(ch==2){
    displayData=ydata2;
    displaySoln=soln2;
  }else if(ch==3){
    displayData=ydata3;
    displaySoln=soln3;
  }

  //The real data
  TGraphErrors *gr_real = new TGraphErrors();
  for(int i=0; i<ydata0.GetNrows(); i++){
    
    gr_real->SetPoint(i, i+0.5,displayData[i][0]);
    gr_real->SetPointError(i, 0, error[i][ch]);

  }
  
  //the Touschek and beam-gas components of the fit
  TH1F* Touschek = new TH1F("Touschek", "Touschek contribution", 1, -0, 1);
  TH1F* Coulomb = new TH1F("Coulomb", "Coulomb contribution", 1, 0, 1);

  
  stringstream ss;
  for(int i=0; i<ydata0.GetNrows(); i++){

    //label each bin with the average beam size during that subrun
    char label[100];
    sprintf(label, "%.0f", BS[i]/n[i]);

    TString l = label;
    l = l +" #mum";

    for(int j=0; j<i; j++) l = l+" ";
    
    Touschek->Fill(l, displaySoln[1][0]*X[i][1]);
    Coulomb->Fill(l,displaySoln[0][0]*X[i][0]);
   
    
  }
  
  //divide the display into sections for each run
  TGraph *grDivider = new TGraph();

  // a vertical line at x=5
  grDivider->SetPoint(0, 5, 100); 
  grDivider->SetPoint(1, 5, -10);

  //a vertical line at x=10
  grDivider->SetPoint(2, 10, -10);
  grDivider->SetPoint(3, 10, 100);
   
  //set drawing options (axis labels and formatting)
  if(draw){
    THStack *hs = new THStack("hs", "");
    Coulomb->SetFillColor(1003);
    Touschek->SetFillColor(1004);
    hs->Add(Coulomb);
    hs->Add(Touschek);
    hs->SetMaximum(90);
    hs->Draw();
    hs->GetYaxis()->SetTitle("Helium-3 tube rate (Hz)");
    hs->GetXaxis()->SetTitle("360mA                                    540mA                                    720mA            ");
    hs->GetXaxis()->SetRangeUser(0, nSubRunLER);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->GetXaxis()->SetLabelSize(0.05);
    gr_real->SetMarkerStyle(20);
    gr_real->Draw("sameP");
    grDivider->Draw("sameL");
    
    if(!sim){
      TPaveText *text = new TPaveText(1.1, 80, 2.9, 90);
      text->AddText("Data");
      text->SetFillColor(kWhite);
      text->SetBorderSize(1);
      text->Draw();
      TPaveText *Beam = new TPaveText(1.1, 60, 2.9, 70);
      Beam->AddText("LER")->SetTextColor(1003);
      Beam->SetFillColor(kWhite);
      Beam->SetBorderSize(1);
      Beam->Draw();
    }else if(sim){
      TPaveText *text = new TPaveText(1.1, 80, 3.9, 90);
      text->AddText("Simulation");
      text->SetFillColor(kWhite);
      text->SetBorderSize(1);
      text->Draw();
    }
    
    
    
  }

  //calculate uncertainties
  TMatrixD Var0 = uncertainty(ydata0, appliedSoln0, X);
  TMatrixD Var1 = uncertainty(ydata1, appliedSoln1, X);
  TMatrixD Var2 = uncertainty(ydata2, appliedSoln2, X);
  TMatrixD Var3 = uncertainty(ydata3, appliedSoln3, X);

  //calculcate chi^2 values
  double chiSq[4] = {0};
  for(int i=0; i<nSubRunLER; i++){
    chiSq[0] += (ydata0[i][0] - appliedSoln0[i][0])**2/(error[i][0])**2;
    chiSq[1] += (ydata1[i][0] - appliedSoln1[i][0])**2/(error[i][1])**2;
    chiSq[2] += (ydata2[i][0] - appliedSoln2[i][0])**2/(error[i][2])**2;
    chiSq[3] += (ydata3[i][0] - appliedSoln3[i][0])**2/(error[i][3])**2;
    //cout<<error[i][1]<<endl;
  }

  int ndf = nSubRunLER-1-X.GetNcols();

  //print out fit parameters with uncertainties
  cout.precision(5);
  cout<<"\tBeamGas\tError\tTous\t\tError\t\tChi^2\tndf\tPval\n";
  cout<<0<<"\t"<<soln0[0][0]<<"\t"<<sqrt(Var0[0][0])<<"\t"<<soln0[1][0]<<"\t"<<sqrt(Var0[1][1])<<"\t"<<chiSq[0]<<"\t"<<ndf<<"\t"<<1-ROOT::Math::ROOT::Math::chisquared_cdf (chiSq[0], ndf)<<endl;
  cout<<1<<"\t"<<soln1[0][0]<<"\t"<<sqrt(Var1[0][0])<<"\t"<<soln1[1][0]<<"\t"<<sqrt(Var1[1][1])<<"\t"<<chiSq[1]<<"\t"<<ndf<<"\t"<<1-ROOT::Math::ROOT::Math::chisquared_cdf (chiSq[1], ndf)<<endl;
  cout<<2<<"\t"<<soln2[0][0]<<"\t"<<sqrt(Var2[0][0])<<"\t"<<soln2[1][0]<<"\t"<<sqrt(Var2[1][1])<<"\t"<<chiSq[2]<<"\t"<<ndf<<"\t"<<1-ROOT::Math::ROOT::Math::chisquared_cdf (chiSq[2], ndf)<<endl;
  cout<<3<<"\t"<<soln3[0][0]<<"\t"<<sqrt(Var3[0][0])<<"\t"<<soln3[1][0]<<"\t"<<sqrt(Var3[1][1])<<"\t"<<chiSq[3]<<"\t"<<ndf<<"\t"<<1-ROOT::Math::ROOT::Math::chisquared_cdf (chiSq[3], ndf)<<endl;


}


