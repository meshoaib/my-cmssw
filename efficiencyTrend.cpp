#include "TROOT.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "rpcTreeFrame.h"
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
using namespace std;

struct theRoll {
  int wheel;
  int station;
  int sublayer;
  int sector;
  int roll;
}; 
  

int main() {
  gROOT->SetStyle("Plain");
  cout << "Open new root tree file" << endl;
  TFile *f;
  TTree *tree;
  // load root file
  f = new TFile("/tmp/piccolo/dummy.root");
  f->cd("rpcRecHitTree");
  tree = (TTree*)gDirectory->Get("rpcTree");
  if (!tree) {
    cout << "Tree does not exist " << endl;
    return 0;
  }
  cout << "Tree opened" << endl;

  // load names of the rolls from the external ascii file
  ifstream myfile("mapRoll.ascii",ios::in);
  int irolls=0;
  int mapId;
  string mapName;
  map<int, string> theRollName;
  map<int,string>::iterator itname;
  while (!myfile.eof()) {
    myfile >> mapId >> mapName;
    irolls++; 
    theRollName[mapId]=mapName;
  }  

  // load the pressure from ascii file
  ifstream pressureFile("myPressureTable.rtf",ios::in);
  int irun;
  float pressure;
  map<int,float> pressureByRun;
  map<int,float>::iterator itp;
  while (!pressureFile.eof()) {
    pressureFile >> irun >> pressure;
    pressureByRun[irun]=pressure;
  }  


  // init the tree
  Init(tree);
  Long64_t nentries = tree->GetEntriesFast();
  //  nentries = 1000;
  cout << "Tree has " << nentries << " entries" << endl;

  // define maps for counters and efficiency
  map<int, int> nDTGoodPointsByRoll;
  map<int, int> nRPCfoundByRoll;
  map<int, double> efficiencyByRoll;
  map<int, double> effErrByRoll;
  map<int,int>::iterator itn;

  map<int, int> nCSCGoodPointsByRoll;
  map<int, int> nEndcapRPCfoundByRoll;
  map<int, double> efficiencyEndcapByRoll;
  map<int, double> effErrEndcapByRoll;

  // open output file
  TFile *outfile;
  outfile = new TFile("efficiencyTrend.root","RECREATE");
  // define root directories
  TDirectory *TDgeneral = outfile->mkdir("general");
  TDirectory *TDpositionStudies = outfile->mkdir("positionStudies");

  // General histograms

  TH1D *hBarrelEfficiency = new TH1D("hBarrelEfficiency","Barrel Efficiency",101, -.5,100.5);
  TH1D *hRB1inEfficiency = new TH1D("hRB1inEfficiency","RB1in Efficiency",101, -.5,100.5);
  TH1D *hRB1outEfficiency = new TH1D("hRB1outEfficiency","RB1out Efficiency",101, -.5,100.5);
  TH1D *hRB2inEfficiency = new TH1D("hRB2inEfficiency","RB2in Efficiency",101, -.5,100.5);
  TH1D *hRB2outEfficiency = new TH1D("hRB2outEfficiency","RB2out Efficiency",101, -.5,100.5);
  TH1D *hRB3Efficiency = new TH1D("hRB3Efficiency","RB3 Efficiency",101, -.5,100.5);
  TH1D *hRB4Efficiency = new TH1D("hRB4Efficiency","RB4 Efficiency",101, -.5,100.5);

  TH1D *hEndcapEfficiency = new TH1D("hEndcapEfficiency","Endcap Efficiency",101, -.5,100.5);
  TH1D *hREm3Efficiency = new TH1D("hREm3Efficiency","RE-3 Efficiency",101, -.5,100.5);
  TH1D *hREm2Efficiency = new TH1D("hREm2Efficiency","RE-2 Efficiency",101, -.5,100.5);
  TH1D *hREm1Efficiency = new TH1D("hREm1Efficiency","RE-1 Efficiency",101, -.5,100.5);
  TH1D *hREp3Efficiency = new TH1D("hREp3Efficiency","RE+3 Efficiency",101, -.5,100.5);
  TH1D *hREp2Efficiency = new TH1D("hREp2Efficiency","RE+2 Efficiency",101, -.5,100.5);
  TH1D *hREp1Efficiency = new TH1D("hREp1Efficiency","RE+1 Efficiency",101, -.5,100.5);

  TH1D *hrecoveryBarrelXYDist = new TH1D("hrecoveryBarrelXYDist","xy dist in Barrel for recovery hits",200, 0.,400.);
  TH1D *hrecoveryBarrelZDist = new TH1D("hrecoveryBarrelZDist","z dist in Barrel for recovery hits",200, 0.,400.);
  TH2D *hrecoveryBarrelZvsXYDist = new TH2D("hrecoveryBarrelZvsXYDist","z vs xy dist in Barrel for recovery hits",100, 0., 400., 100, 0.,400.);

  TH1D *hrecoveryEcapDeltaPhi = new TH1D("hrecoveryEcapDeltaPhi","Delta phi Endcap for recovery hits",160, -40.,40.);
  TH1D *hrecoveryEcapDeltaR = new TH1D("hrecoveryEcapDelatR","Delta R Endcap for recovery hits",200, 0.,400.);
  TH2D *hrecoveryEcapDeltaRvsDeltaPhi = new TH2D("hrecoveryEcapDeltaRvsDetalPhi","Delta R vs Detal phi Endcap for recovery hits",160, -40., 40., 200, 0.,400.);

  TH1D *hEcapDeltaPhi = new TH1D("hEcapDeltaPhi","Delta phi Endcap",160, -20.,20.);
  TH1D *hEcapDeltaR = new TH1D("hEcapDelatR","Delta R Endcap",200, 0.,400.);
  TH2D *hEcapDeltaPhivsPhi_REm1_R2 = new TH2D("hEcapDetalPhivsPhi_REm1_R2","Delta Phi vs phi RE-1 ring 2",360, -180., 180., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsPhi_REm1_R3 = new TH2D("hEcapDetalPhivsPhi_REm1_R3","Delta Phi vs phi RE-1 ring 3",360, -180., 180., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsPhi_REm2_R2 = new TH2D("hEcapDetalPhivsPhi_REm2_R2","Delta Phi vs phi RE-2 ring 2",360, -180., 180., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsPhi_REm2_R3 = new TH2D("hEcapDetalPhivsPhi_REm2_R3","Delta Phi vs phi RE-2 ring 3",360, -180., 180., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsPhi_REm3_R2 = new TH2D("hEcapDetalPhivsPhi_REm3_R2","Delta Phi vs phi RE-3 ring 2",360, -180., 180., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsPhi_REm3_R3 = new TH2D("hEcapDetalPhivsPhi_REm3_R3","Delta Phi vs phi RE-3 ring 3",360, -180., 180., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsPhi_REp1_R2 = new TH2D("hEcapDetalPhivsPhi_REp1_R2","Delta Phi vs phi RE+1 ring 2",360, -180., 180., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsPhi_REp1_R3 = new TH2D("hEcapDetalPhivsPhi_REp1_R3","Delta Phi vs phi RE+1 ring 3",360, -180., 180., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsPhi_REp2_R2 = new TH2D("hEcapDetalPhivsPhi_REp2_R2","Delta Phi vs phi RE+2 ring 2",360, -180., 180., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsPhi_REp2_R3 = new TH2D("hEcapDetalPhivsPhi_REp2_R3","Delta Phi vs phi RE+2 ring 3",360, -180., 180., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsPhi_REp3_R2 = new TH2D("hEcapDetalPhivsPhi_REp3_R2","Delta Phi vs phi RE+3 ring 2",360, -180., 180., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsPhi_REp3_R3 = new TH2D("hEcapDetalPhivsPhi_REp3_R3","Delta Phi vs phi RE+3 ring 3",360, -180., 180., 160, -20.,20.);

  TH2D *hEcapDeltaPhivsLocalx_REm1_R2 =  new TH2D("hEcapDeltaPhivsLocalx_REm1_R2","Delta Phi vs Local x RE-1 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocalx_REm1_R3 =  new TH2D("hEcapDeltaPhivsLocalx_REm1_R3","Delta Phi vs Local x RE-1 ring 3",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocalx_REm2_R2 =  new TH2D("hEcapDeltaPhivsLocalx_REm2_R2","Delta Phi vs Local x RE-2 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocalx_REm2_R3 =  new TH2D("hEcapDeltaPhivsLocalx_REm2_R3","Delta Phi vs Local x RE-2 ring 3",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocalx_REm3_R2 =  new TH2D("hEcapDeltaPhivsLocalx_REm3_R2","Delta Phi vs Local x RE-3 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocalx_REm3_R3 =  new TH2D("hEcapDeltaPhivsLocalx_REm3_R3","Delta Phi vs Local x RE-3 ring 3",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocalx_REp1_R2 =  new TH2D("hEcapDeltaPhivsLocalx_REp1_R2","Delta Phi vs Local x RE+1 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocalx_REp1_R3 =  new TH2D("hEcapDeltaPhivsLocalx_REp1_R3","Delta Phi vs Local x RE+1 ring 3",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocalx_REp2_R2 =  new TH2D("hEcapDeltaPhivsLocalx_REp2_R2","Delta Phi vs Local x RE+2 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocalx_REp2_R3 =  new TH2D("hEcapDeltaPhivsLocalx_REp2_R3","Delta Phi vs Local x RE+2 ring 3",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocalx_REp3_R2 =  new TH2D("hEcapDeltaPhivsLocalx_REp3_R2","Delta Phi vs Local x RE+3 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocalx_REp3_R3 =  new TH2D("hEcapDeltaPhivsLocalx_REp3_R3","Delta Phi vs Local x RE+3 ring 3",100, -50., 50., 160, -20.,20.);

  TH2D *hEcapDeltaPhivsLocaly_REm1_R2 =  new TH2D("hEcapDeltaPhivsLocaly_REm1_R2","Delta Phi vs Local y RE-1 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocaly_REm1_R3 =  new TH2D("hEcapDeltaPhivsLocaly_REm1_R3","Delta Phi vs Local y RE-1 ring 3",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocaly_REm2_R2 =  new TH2D("hEcapDeltaPhivsLocaly_REm2_R2","Delta Phi vs Local y RE-2 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocaly_REm2_R3 =  new TH2D("hEcapDeltaPhivsLocaly_REm2_R3","Delta Phi vs Local y RE-2 ring 3",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocaly_REm3_R2 =  new TH2D("hEcapDeltaPhivsLocaly_REm3_R2","Delta Phi vs Local y RE-3 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocaly_REm3_R3 =  new TH2D("hEcapDeltaPhivsLocaly_REm3_R3","Delta Phi vs Local y RE-3 ring 3",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocaly_REp1_R2 =  new TH2D("hEcapDeltaPhivsLocaly_REp1_R2","Delta Phi vs Local y RE+1 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocaly_REp1_R3 =  new TH2D("hEcapDeltaPhivsLocaly_REp1_R3","Delta Phi vs Local y RE+1 ring 3",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocaly_REp2_R2 =  new TH2D("hEcapDeltaPhivsLocaly_REp2_R2","Delta Phi vs Local y RE+2 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocaly_REp2_R3 =  new TH2D("hEcapDeltaPhivsLocaly_REp2_R3","Delta Phi vs Local y RE+2 ring 3",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocaly_REp3_R2 =  new TH2D("hEcapDeltaPhivsLocaly_REp3_R2","Delta Phi vs Local y RE+3 ring 2",100, -50., 50., 160, -20.,20.);
  TH2D *hEcapDeltaPhivsLocaly_REp3_R3 =  new TH2D("hEcapDeltaPhivsLocaly_REp3_R3","Delta Phi vs Local y RE+3 ring 3",100, -50., 50., 160, -20.,20.);
  // trend Analysis
  int runRangeMin, runRangeMax, runRangeBins;
  runRangeMin = 146000;
  runRangeMax = 147600;
  runRangeBins = runRangeMax-runRangeMin;
  //  runRangeBins = 500;

  TH1F *hpressureByRun = new TH1F("hpressureByRun","Environmental pressure vs Run",runRangeBins,runRangeMin,runRangeMax);

  TH2D *hEffVsPressure_RB1in = new TH2D("hEffVsPressure_RB1in","Eff vs pressure RB1in",100,930,980,100,0.8,1.);
  TH2D *hEffVsPressure_RB1out = new TH2D("hEffVsPressure_RB1out","Eff vs pressure RB1out",100,930,980,100,0.8,1.);
  TH2D *hEffVsPressure_RB2in = new TH2D("hEffVsPressure_RB2in","Eff vs pressure RB2in",100,930,980,100,0.8,1.);
  TH2D *hEffVsPressure_RB2out = new TH2D("hEffVsPressure_RB2out","Eff vs pressure RB2out",100,930,980,100,0.8,1.);
  TH2D *hEffVsPressure_RB3 = new TH2D("hEffVsPressure_RB3","Eff vs pressure RB3",100,930,980,100,0.8,1.);
  TH2D *hEffVsPressure_RB4 = new TH2D("hEffVsPressure_RB4","Eff vs pressure RB4",100,930,980,100,0.8,1.);

  TH2D *hEffVsPressure_REm1 = new TH2D("hEffVsPressure_REm1","Eff vs pressure RE-1",100,930,980,100,0.8,1.);
  TH2D *hEffVsPressure_REm2 = new TH2D("hEffVsPressure_REm2","Eff vs pressure RE-2",100,930,980,100,0.8,1.);
  TH2D *hEffVsPressure_REm3 = new TH2D("hEffVsPressure_REm3","Eff vs pressure RE-3",100,930,980,100,0.8,1.);
  TH2D *hEffVsPressure_REp1 = new TH2D("hEffVsPressure_REp1","Eff vs pressure RE+1",100,930,980,100,0.8,1.);
  TH2D *hEffVsPressure_REp2 = new TH2D("hEffVsPressure_REp2","Eff vs pressure RE+2",100,930,980,100,0.8,1.);
  TH2D *hEffVsPressure_REp3 = new TH2D("hEffVsPressure_REp3","Eff vs pressure RE+3",100,930,980,100,0.8,1.);

  // barrel histos
  TH1F *hExpRB1inVsRun = new TH1F("hExpRB1inVsRun","expected RB1in vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hExpRB1outVsRun = new TH1F("hExpRB1outVsRun","expected RB1out vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hExpRB2inVsRun = new TH1F("hExpRB2inVsRun","expected RB2in Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hExpRB2outVsRun = new TH1F("hExpRB2outVsRun","expected RB2out vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hExpRB3VsRun = new TH1F("hExpRB3VsRun","expected RB3 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hExpRB4VsRun = new TH1F("hExpRB4VsRun","expected RB4 vs Run",runRangeBins,runRangeMin,runRangeMax);
  
  TH1F *hFoundRB1inVsRun = new TH1F("hFoundRB1inVsRun","found RB1in vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hFoundRB1outVsRun = new TH1F("hFoundRB1outVsRun","found RB1out vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hFoundRB2inVsRun = new TH1F("hFoundRB2inVsRun","found RB2in vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hFoundRB2outVsRun = new TH1F("hFoundRB2outVsRun","found RB2out vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hFoundRB3VsRun = new TH1F("hFoundRB3VsRun","found RB3 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hFoundRB4VsRun = new TH1F("hFoundRB4VsRun","found RB4 vs Run",runRangeBins,runRangeMin,runRangeMax);
  
  TH1F *hEffRB1inVsRun = new TH1F("hEffRB1inVsRun","Eff RB1in vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hEffRB1outVsRun = new TH1F("hEffRB1outVsRun","Eff RB1out vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hEffRB2inVsRun = new TH1F("hEffRB2inVsRun","Eff RB2in vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hEffRB2outVsRun = new TH1F("hEffRB2outVsRun","Eff RB2out vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hEffRB3VsRun = new TH1F("hEffRB3VsRun","Eff RB3 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hEffRB4VsRun = new TH1F("hEffRB4VsRun","Eff RB4 vs Run",runRangeBins,runRangeMin,runRangeMax);
  
  // endcap histos
  TH1F *hExpREp1VsRun = new TH1F("hExpREp1VsRun","expected RE+1 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hExpREp2VsRun = new TH1F("hExpREp2VsRun","expected RE+2 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hExpREp3VsRun = new TH1F("hExpREp3VsRun","expected RE+3 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hExpREm1VsRun = new TH1F("hExpREm1VsRun","expected RE-1 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hExpREm2VsRun = new TH1F("hExpREm2VsRun","expected RE-2 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hExpREm3VsRun = new TH1F("hExpREm3VsRun","expected RE-3 vs Run",runRangeBins,runRangeMin,runRangeMax);
     
  TH1F *hFoundREp1VsRun = new TH1F("hFoundREp1VsRun","found RE+1 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hFoundREp2VsRun = new TH1F("hFoundREp2VsRun","found RE+2 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hFoundREp3VsRun = new TH1F("hFoundREp3VsRun","found RE+3 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hFoundREm1VsRun = new TH1F("hFoundREm1VsRun","found RE-1 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hFoundREm2VsRun = new TH1F("hFoundREm2VsRun","found RE-2 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hFoundREm3VsRun = new TH1F("hFoundREm3VsRun","found RE-3 vs Run",runRangeBins,runRangeMin,runRangeMax);

  TH1F *hEffREp1VsRun = new TH1F("hEffREp1VsRun","Eff RE+1 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hEffREp2VsRun = new TH1F("hEffREp2VsRun","Eff RE+2 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hEffREp3VsRun = new TH1F("hEffREp3VsRun","Eff RE+3 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hEffREm1VsRun = new TH1F("hEffREm1VsRun","Eff RE-1 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hEffREm2VsRun = new TH1F("hEffREm2VsRun","Eff RE-2 vs Run",runRangeBins,runRangeMin,runRangeMax);
  TH1F *hEffREm3VsRun = new TH1F("hEffREm3VsRun","Eff RE-3 vs Run",runRangeBins,runRangeMin,runRangeMax);

  // Global 2D Endcap Efficiency Map
  TH2D *h2DExpectedREm1 = new TH2D("h2DExpectedREm1","Expected REm1",400,-800,800,400,-800,800);
  TH2D *h2DExpectedREm2 = new TH2D("h2DExpectedREm2","Expected REm2",400,-800,800,400,-800,800);
  TH2D *h2DExpectedREm3 = new TH2D("h2DExpectedREm3","Expected REm3",400,-800,800,400,-800,800);
  TH2D *h2DExpectedREp1 = new TH2D("h2DExpectedREp1","Expected REp1",400,-800,800,400,-800,800);
  TH2D *h2DExpectedREp2 = new TH2D("h2DExpectedREp2","Expected REp2",400,-800,800,400,-800,800);
  TH2D *h2DExpectedREp3 = new TH2D("h2DExpectedREp3","Expected REp3",400,-800,800,400,-800,800);
  
  TH2D *h2DFoundREm1 = new TH2D("h2DFoundREm1","Found REm1",400,-800,800,400,-800,800);
  TH2D *h2DFoundREm2 = new TH2D("h2DFoundREm2","Found REm2",400,-800,800,400,-800,800);
  TH2D *h2DFoundREm3 = new TH2D("h2DFoundREm3","Found REm3",400,-800,800,400,-800,800);
  TH2D *h2DFoundREp1 = new TH2D("h2DFoundREp1","Found REp1",400,-800,800,400,-800,800);
  TH2D *h2DFoundREp2 = new TH2D("h2DFoundREp2","Found REp2",400,-800,800,400,-800,800);
  TH2D *h2DFoundREp3 = new TH2D("h2DFoundREp3","Found REp3",400,-800,800,400,-800,800);

  TH2D *h2DEffREm1 = new TH2D("h2DEffREm1","Eff REm1",400,-800,800,400,-800,800);
  TH2D *h2DEffREm2 = new TH2D("h2DEffREm2","Eff REm2",400,-800,800,400,-800,800);
  TH2D *h2DEffREm3 = new TH2D("h2DEffREm3","Eff REm3",400,-800,800,400,-800,800);
  TH2D *h2DEffREp1 = new TH2D("h2DEffREp1","Eff REp1",400,-800,800,400,-800,800);
  TH2D *h2DEffREp2 = new TH2D("h2DEffREp2","Eff REp2",400,-800,800,400,-800,800);
  TH2D *h2DEffREp3 = new TH2D("h2DEffREp3","Eff REp3",400,-800,800,400,-800,800);


  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //    Long64_t ientry = LoadTree(jentry);
    //    if (ientry < 0) break;
    int nb = tree->GetEntry(jentry);  
    //    cout << "Entry " << jentry << "Run " << run << "  event " << evt << "   RPC clusters " << nClusters << endl;

    // Barrel Efficiency
    for (int i=0;i<nDTPoints; i++) {
      int tstation = DTPoint_Station[i];
      int tsublayer = DTPoint_SubLayer[i];
      int twheel = DTPoint_Wheel[i];
      int tsector = DTPoint_Sector[i];
      int trawid = DTPoint_RawId[i];
      int ilay = (tstation-1)*2+tsublayer;  // convert station .. sublayer in a number between 1 and 6
      if (tstation==3) ilay=5;
      if (tstation==4) ilay=6;

     // this part is just to correlate with DT segments
      int tdimension = 0;
      int segmentPointer = -1;
      for (int j=0; j<nDTSegments; j++) {
	int tDTsector = DTSector[j];
	if (tDTsector==13) tDTsector=4; // handle exceptions in DT sector numbering
	if (tDTsector==14) tDTsector=10;  // handle exceptions in DT sector numbering
	if (tstation==DTStation[j] && twheel==DTWheel[j] && tsector==tDTsector) {
	  if (tstation<4) tdimension = DTdimension[j];
	  if (tstation==4) tdimension=3;
	  segmentPointer=j;
	} // end comparison point - segment
      } // end loop on DT segments

      // selection and Histograms filling
      bool DTsegmentGood = false;
      bool BarrelFiducialCut = false;
      bool BarrelTightFiducialCut = false;
      if (segmentPointer>-1 && tdimension>2) DTsegmentGood=true; // quality of associted segment
      if (DTPoint_XdistToBorder[i]>0 && DTPoint_YdistToBorder[i]>0) BarrelFiducialCut=true; // extrapolated point inside RPC region
      if (DTPoint_XdistToBorder[i]>5 && DTPoint_YdistToBorder[i]>5) BarrelTightFiducialCut=true; // extrapolated point inside Tight RPC region

      if (DTsegmentGood && BarrelFiducialCut) { // extrapolation used for analysis
	nDTGoodPointsByRoll[trawid]+=1;
	if (ilay==1) hExpRB1inVsRun->Fill(run);
	if (ilay==2) hExpRB1outVsRun->Fill(run);
	if (ilay==3) hExpRB2inVsRun->Fill(run);
	if (ilay==4) hExpRB2outVsRun->Fill(run);
	if (ilay==5) hExpRB3VsRun->Fill(run);
	if (ilay==6) hExpRB4VsRun->Fill(run);
	
	bool rpcHit=false;
	if (DTPoint_associatedClusterMulti[i]>0) {
	  rpcHit=true;
	} else {
	  for (int r=0; r<nClusters; r++) {
	    if (rpcStation[r]==tstation && rpcSubLayer[r]==tsublayer) {
	      float xyDist = sqrt(pow(xRPC[r]-xDTPoint[i],2)+pow(yRPC[r]-yDTPoint[i],2));
	      float zDist = fabs(zRPC[r]-zDTPoint[i]);
	      if (xyDist<5 && zDist<100) rpcHit=true;
	      hrecoveryBarrelXYDist->Fill(xyDist);
	      hrecoveryBarrelZDist->Fill(zDist);
	      hrecoveryBarrelZvsXYDist->Fill(xyDist,zDist);
	      //	      if (rpcHit) cout << "Run " << run << "  event " << evt << "  Sector/wheel/station/sublayer " << tsector << " " << twheel << " " << tstation << " " << tsublayer << endl;
	    }
	  } 
	} 

	if (rpcHit) { // RPC found
	  nRPCfoundByRoll[trawid]+=1;
	  if (ilay==1) hFoundRB1inVsRun->Fill(run);
	  if (ilay==2) hFoundRB1outVsRun->Fill(run);
	  if (ilay==3) hFoundRB2inVsRun->Fill(run);
	  if (ilay==4) hFoundRB2outVsRun->Fill(run);
	  if (ilay==5) hFoundRB3VsRun->Fill(run);
	  if (ilay==6) hFoundRB4VsRun->Fill(run);
	} // end RPC found
      } // end if is extrapolation used in analysis
    } // end loop on DTPoints 

    // Endcap Efficiency
    for (int i=0;i<nCSCPoints; i++) {
      int tregion = CSCPoint_Region[i];
      int tstation = CSCPoint_Station[i];
      int tring = CSCPoint_Ring[i];
      int trawid = CSCPoint_RawId[i];
      int troll = CSCPoint_Roll[i];
      int tchamber = CSCPoint_Chamber[i];
      bool EndcapFiducialCut = false;
      bool EndcapTightFiducialCut = false;
      if (CSCPoint_XdistToBorder[i]>0 && CSCPoint_YdistToBorder[i]>0) EndcapFiducialCut=true; // extrapolated point inside RPC region
      if (CSCPoint_XdistToBorder[i]>0 && CSCPoint_YdistToBorder[i]>0) EndcapTightFiducialCut=true; // extrapolated point inside RPC region

      if (EndcapFiducialCut) { // extrapolation used for analysis
      //      if (EndcapFiducialCut) { // extrapolation used for analysis
	//	if (EndcapTightFiducialCut && run>146807)  nCSCGoodPointsByRoll[trawid]+=1;
	if (EndcapTightFiducialCut)  nCSCGoodPointsByRoll[trawid]+=1;
	if (tregion*tstation==-3) {
	  if (EndcapTightFiducialCut) hExpREm3VsRun->Fill(run);
	  h2DExpectedREm3->Fill(xCSCPoint[i],yCSCPoint[i]);
	}
	if (tregion*tstation==-2) {
	  if (EndcapTightFiducialCut) hExpREm2VsRun->Fill(run);
	  h2DExpectedREm2->Fill(xCSCPoint[i],yCSCPoint[i]);
	}
	if (tregion*tstation==-1) {
	  if (EndcapTightFiducialCut) hExpREm1VsRun->Fill(run);
	  h2DExpectedREm1->Fill(xCSCPoint[i],yCSCPoint[i]);
	}
	if (tregion*tstation==3) {
	  if (EndcapTightFiducialCut) hExpREp3VsRun->Fill(run);
	  h2DExpectedREp3->Fill(xCSCPoint[i],yCSCPoint[i]);
	}
	if (tregion*tstation==2) {
	  if (EndcapTightFiducialCut) hExpREp2VsRun->Fill(run);
	  h2DExpectedREp2->Fill(xCSCPoint[i],yCSCPoint[i]);
	}
	if (tregion*tstation==1){
	  if (EndcapTightFiducialCut) hExpREp1VsRun->Fill(run);
	  h2DExpectedREp1->Fill(xCSCPoint[i],yCSCPoint[i]);
	}
	// check angular displacement
	for (int r=0; r<nClusters; r++) {
	  if (rpcStation[r]==tstation && rpcRegion[r]==tregion) {
	    
	    int sideCSC = 0;
	    if (yCSCPoint[i]>=0 && xCSCPoint[i]>=0) sideCSC=1;
	    if (yCSCPoint[i]>0 && xCSCPoint[i]<0) sideCSC=2;
	    if (yCSCPoint[i]<=0 && xCSCPoint[i]<=0) sideCSC=3;
	    if (yCSCPoint[i]<0 && xCSCPoint[i]>0) sideCSC=4;
	    int sideRPC = 0;
	    if (yRPC[r]>=0 && xRPC[r]>=0) sideRPC=1;
	    if (yRPC[r]>0 && xRPC[r]<0) sideRPC=2;
	    if (yRPC[r]<=0 && xRPC[r]<=0) sideRPC=3;
	    if (yRPC[r]<0 && xRPC[r]>0) sideRPC=4;
	    float rCSC= sqrt(pow(yCSCPoint[i],2)+pow(xCSCPoint[i],2));
	    float rRPC= sqrt(pow(yRPC[r],2)+pow(xRPC[r],2));
	    
	    float phiCSC = asin(yCSCPoint[i]/rCSC)*180./3.1415;
	    float phiRPC = asin(yRPC[r]/rRPC)*180./3.1415;
	    if (sideCSC==3) phiCSC=-180.-phiCSC;
	    if (sideCSC==2) phiCSC=180.-phiCSC;
	    if (sideRPC==3) phiRPC=-180.-phiRPC;
	    if (sideRPC==2) phiRPC=180.-phiRPC;
	    float deltaPhi = phiCSC-phiRPC;
	    float deltaR = rCSC-rRPC;
	    hEcapDeltaPhi->Fill(deltaPhi);
	    if (tregion==1) {
	      if (tstation==1 && tring==2) {
		hEcapDeltaPhivsPhi_REp1_R2->Fill(phiCSC,deltaPhi); 
		hEcapDeltaPhivsLocalx_REp1_R2->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REp1_R2->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	      if (tstation==1 && tring==3) {
		hEcapDeltaPhivsPhi_REp1_R3->Fill(phiCSC,deltaPhi);
		hEcapDeltaPhivsLocalx_REp1_R3->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REp1_R3->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	      if (tstation==2 && tring==2) {
		hEcapDeltaPhivsPhi_REp2_R2->Fill(phiCSC,deltaPhi);
		hEcapDeltaPhivsLocalx_REp2_R2->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REp2_R2->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	      if (tstation==2 && tring==3) {
		hEcapDeltaPhivsPhi_REp2_R3->Fill(phiCSC,deltaPhi);
		hEcapDeltaPhivsLocalx_REp2_R3->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REp2_R3->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	      if (tstation==3 && tring==2) {
		hEcapDeltaPhivsPhi_REp3_R2->Fill(phiCSC,deltaPhi);
		hEcapDeltaPhivsLocalx_REp3_R2->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REp3_R2->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	      if (tstation==3 && tring==3) {
		hEcapDeltaPhivsPhi_REp3_R3->Fill(phiCSC,deltaPhi);
		hEcapDeltaPhivsLocalx_REp3_R3->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REp3_R3->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	    }
	    hEcapDeltaPhi->Fill(deltaPhi);
	    if (tregion==-1) {
	      if (tstation==1 && tring==2) {
		hEcapDeltaPhivsPhi_REm1_R2->Fill(phiCSC,deltaPhi);
		hEcapDeltaPhivsLocalx_REm1_R2->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REm1_R2->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	      if (tstation==1 && tring==3) {
		hEcapDeltaPhivsPhi_REm1_R3->Fill(phiCSC,deltaPhi);
		hEcapDeltaPhivsLocalx_REm1_R3->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REm1_R3->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	      if (tstation==2 && tring==2) {
		hEcapDeltaPhivsPhi_REm2_R2->Fill(phiCSC,deltaPhi);
		hEcapDeltaPhivsLocalx_REm2_R2->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REm2_R2->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	      if (tstation==2 && tring==3) {
		hEcapDeltaPhivsPhi_REm2_R3->Fill(phiCSC,deltaPhi);
		hEcapDeltaPhivsLocalx_REm2_R3->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REm2_R3->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	      if (tstation==3 && tring==2) {
		hEcapDeltaPhivsPhi_REm3_R2->Fill(phiCSC,deltaPhi);
		hEcapDeltaPhivsLocalx_REm3_R2->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REm3_R2->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	      if (tstation==3 && tring==3) {
		hEcapDeltaPhivsPhi_REm3_R3->Fill(phiCSC,deltaPhi);
		hEcapDeltaPhivsLocalx_REm3_R3->Fill(xLocalCSCPoint[i],deltaPhi); 
		hEcapDeltaPhivsLocaly_REm3_R3->Fill(yLocalCSCPoint[i],deltaPhi);
	      }
	    }
	    hEcapDeltaR->Fill(deltaR);
	  }
	  
	}
	bool ErpcHit=false;
	bool EndcapRecovery=true;
	if (CSCPoint_associatedClusterMulti[i]>0) {
	  ErpcHit=true;
	} else { // try to recover 
	  for (int r=0; r<nClusters; r++) {
	    if (rpcStation[r]==tstation && rpcRegion[r]==tregion) {

	      int sideCSC = 0;
	      if (yCSCPoint[i]>=0 && xCSCPoint[i]>=0) sideCSC=1;
	      if (yCSCPoint[i]>0 && xCSCPoint[i]<0) sideCSC=2;
	      if (yCSCPoint[i]<=0 && xCSCPoint[i]<=0) sideCSC=3;
	      if (yCSCPoint[i]<0 && xCSCPoint[i]>0) sideCSC=4;
	      int sideRPC = 0;
	      if (yRPC[r]>=0 && xRPC[r]>=0) sideRPC=1;
	      if (yRPC[r]>0 && xRPC[r]<0) sideRPC=2;
	      if (yRPC[r]<=0 && xRPC[r]<=0) sideRPC=3;
	      if (yRPC[r]<0 && xRPC[r]>0) sideRPC=4;
	      float rCSC= sqrt(pow(yCSCPoint[i],2)+pow(xCSCPoint[i],2));
	      float rRPC= sqrt(pow(yRPC[r],2)+pow(xRPC[r],2));

	      float phiCSC = asin(yCSCPoint[i]/rCSC)*180./3.1415;
	      float phiRPC = asin(yRPC[r]/rRPC)*180./3.1415;
	      if (sideCSC==3) phiCSC=-180.-phiCSC;
	      if (sideCSC==4) phiCSC=180.-phiCSC;
	      if (sideRPC==3) phiRPC=-180.-phiRPC;
	      if (sideRPC==4) phiRPC=180.-phiRPC;
	      float deltaPhi = phiCSC-phiRPC;
	      float deltaR = rCSC-rRPC;
	      hrecoveryEcapDeltaPhi->Fill(deltaPhi);
	      hrecoveryEcapDeltaR->Fill(deltaR);
	      hrecoveryEcapDeltaRvsDeltaPhi->Fill(deltaPhi,deltaR);
	      if (EndcapRecovery && fabs(deltaPhi)<1. && fabs(deltaR)<50) ErpcHit=true;
	    }
	  } 
	  
	} 

	if (ErpcHit) { // RPC found
	  //	  if (EndcapTightFiducialCut && run>146807) nEndcapRPCfoundByRoll[trawid]+=1;
	  if (EndcapTightFiducialCut) nEndcapRPCfoundByRoll[trawid]+=1;
	  if (tregion*tstation==-3) {
	    if (EndcapTightFiducialCut) hFoundREm3VsRun->Fill(run);
	    h2DFoundREm3->Fill(xCSCPoint[i],yCSCPoint[i]);
	  }
	  if (tregion*tstation==-2) {
	    if (EndcapTightFiducialCut) hFoundREm2VsRun->Fill(run);
	    h2DFoundREm2->Fill(xCSCPoint[i],yCSCPoint[i]);
	  }
	  if (tregion*tstation==-1) {
	    if (EndcapTightFiducialCut) hFoundREm1VsRun->Fill(run);
	    h2DFoundREm1->Fill(xCSCPoint[i],yCSCPoint[i]);
	  }
	  if (tregion*tstation==3) {
	    if (EndcapTightFiducialCut) hFoundREp3VsRun->Fill(run);
	    h2DFoundREp3->Fill(xCSCPoint[i],yCSCPoint[i]);
	  }
	  if (tregion*tstation==2) {
	    if (EndcapTightFiducialCut) hFoundREp2VsRun->Fill(run);
	    h2DFoundREp2->Fill(xCSCPoint[i],yCSCPoint[i]);
	  }
	  if (tregion*tstation==1) {
	    if (EndcapTightFiducialCut) hFoundREp1VsRun->Fill(run);
	    h2DFoundREp1->Fill(xCSCPoint[i],yCSCPoint[i]);
	  }
	} // end RPC found
      } // end if is extrapolation used in analysis
    } // end loop on CSC Points
  } // end loop on events
  
  // *********  analysis of efficiency  *********

  // pressure plot
  for (itp=pressureByRun.begin(); itp!=pressureByRun.end(); itp++) {
    int irun = (*itp).first;
    float pressure = (*itp).second;
    hpressureByRun->Fill(irun,pressure);
  } 

  // barrel
  for ( itn=nDTGoodPointsByRoll.begin() ; itn != nDTGoodPointsByRoll.end(); itn++ ) {
    double eff = 0;
    int n = nDTGoodPointsByRoll[(*itn).first];
    int nRPC = nRPCfoundByRoll[(*itn).first];
    if (n>0) eff = 100.*float(nRPC)/n;
    hBarrelEfficiency->Fill(eff);
    int ilay=0;
    string a = theRollName[(*itn).first];
    bool isRB1in = false;
    bool isRB1out = false;
    bool isRB2in = false;
    bool isRB2out = false;
    bool isRB3 = false;
    bool isRB4 = false;
    if (a.find("RB1in",0)!=string::npos) isRB1in=true;
    if (a.find("RB1out",0)!=string::npos) isRB1out=true;
    if (a.find("RB2in",0)!=string::npos) isRB2in=true;
    if (a.find("RB2out",0)!=string::npos) isRB2out=true;
    if (a.find("RB3",0)!=string::npos) isRB3=true;
    if (a.find("RB4",0)!=string::npos) isRB4=true;

    //   cout << theRollName[(*itn).first] << " " << n << " " << eff << endl;
    if (isRB1in) hRB1inEfficiency->Fill(eff);
    if (isRB1out) hRB1outEfficiency->Fill(eff);
    if (isRB2in) hRB2inEfficiency->Fill(eff);
    if (isRB2out) hRB2outEfficiency->Fill(eff);
    if (isRB3) hRB3Efficiency->Fill(eff);
    if (isRB4) hRB4Efficiency->Fill(eff);

  }
  hFoundRB1inVsRun->Sumw2();
  hExpRB1inVsRun->Sumw2();
  hFoundRB1outVsRun->Sumw2();
  hExpRB1outVsRun->Sumw2();
  hFoundRB2inVsRun->Sumw2();
  hExpRB2inVsRun->Sumw2();
  hFoundRB2outVsRun->Sumw2();
  hExpRB2outVsRun->Sumw2();
  hFoundRB3VsRun->Sumw2();
  hExpRB3VsRun->Sumw2();
  hFoundRB4VsRun->Sumw2();
  hExpRB4VsRun->Sumw2();
  
  hEffRB1inVsRun->Divide(hFoundRB1inVsRun,hExpRB1inVsRun,1,1,"B");
  hEffRB1outVsRun->Divide(hFoundRB1outVsRun,hExpRB1outVsRun,1,1,"B");
  hEffRB2inVsRun->Divide(hFoundRB2inVsRun,hExpRB2inVsRun,1,1,"B");
  hEffRB2outVsRun->Divide(hFoundRB2outVsRun,hExpRB2outVsRun,1,1,"B");
  hEffRB3VsRun->Divide(hFoundRB3VsRun,hExpRB3VsRun,1,1,"B");
  hEffRB4VsRun->Divide(hFoundRB4VsRun,hExpRB4VsRun,1,1,"B");

  // endcap
  for ( itn=nCSCGoodPointsByRoll.begin() ; itn != nCSCGoodPointsByRoll.end(); itn++ ) {
    double eff = -1.;
    int n = nCSCGoodPointsByRoll[(*itn).first];
    int nRPC = nEndcapRPCfoundByRoll[(*itn).first];
    if (n>10) {
      eff = 100.*float(nRPC)/n;
      hEndcapEfficiency->Fill(eff);
    }
    int ilay=0;
    string a = theRollName[(*itn).first];
    bool isREm3 = false;
    bool isREm2 = false;
    bool isREm1 = false;
    bool isREp3 = false;
    bool isREp2 = false;
    bool isREp1 = false;
    if (a.find("RE-3",0)!=string::npos) isREm3=true;
    if (a.find("RE-2",0)!=string::npos) isREm2=true;
    if (a.find("RE-1",0)!=string::npos) isREm1=true;
    if (a.find("RE+3",0)!=string::npos) isREp3=true;
    if (a.find("RE+2",0)!=string::npos) isREp2=true;
    if (a.find("RE+1",0)!=string::npos) isREp1=true;

    cout << theRollName[(*itn).first] << " " << n << " " << eff << endl;
    if (isREm3) hREm3Efficiency->Fill(eff);
    if (isREm2) hREm2Efficiency->Fill(eff);
    if (isREm1) hREm1Efficiency->Fill(eff);
    if (isREp3) hREp3Efficiency->Fill(eff);
    if (isREp2) hREp2Efficiency->Fill(eff);
    if (isREp1) hREp1Efficiency->Fill(eff);     
  }  
  hFoundREm1VsRun->Sumw2();
  hExpREm1VsRun->Sumw2();
  hFoundREm2VsRun->Sumw2();
  hExpREm2VsRun->Sumw2();
  hFoundREm3VsRun->Sumw2();
  hExpREm3VsRun->Sumw2();
  hFoundREp1VsRun->Sumw2();
  hExpREp1VsRun->Sumw2();
  hFoundREp2VsRun->Sumw2();
  hExpREp2VsRun->Sumw2();
  hFoundREp3VsRun->Sumw2();
  hExpREp3VsRun->Sumw2();
  
  hEffREm1VsRun->Divide(hFoundREm1VsRun,hExpREm1VsRun,1,1,"B");
  hEffREm2VsRun->Divide(hFoundREm2VsRun,hExpREm2VsRun,1,1,"B");
  hEffREm3VsRun->Divide(hFoundREm3VsRun,hExpREm3VsRun,1,1,"B");
  hEffREp1VsRun->Divide(hFoundREp1VsRun,hExpREp1VsRun,1,1,"B");
  hEffREp2VsRun->Divide(hFoundREp2VsRun,hExpREp2VsRun,1,1,"B");
  hEffREp3VsRun->Divide(hFoundREp3VsRun,hExpREp3VsRun,1,1,"B");

  // 2D Endcap Eff plot (ratio)
  h2DEffREm1->Divide(h2DFoundREm1,h2DExpectedREm1,1,1,"B");
  h2DEffREm2->Divide(h2DFoundREm2,h2DExpectedREm2,1,1,"B");
  h2DEffREm3->Divide(h2DFoundREm3,h2DExpectedREm3,1,1,"B");
  h2DEffREp1->Divide(h2DFoundREp1,h2DExpectedREp1,1,1,"B");
  h2DEffREp2->Divide(h2DFoundREp2,h2DExpectedREp2,1,1,"B");
  h2DEffREp3->Divide(h2DFoundREp3,h2DExpectedREp3,1,1,"B");

  // correlation efficiency - pressure
  int nbins = hEffRB1inVsRun->GetXaxis()->GetNbins();
  cout << " N x bins " << nbins << endl;
  for (int xbin=1; xbin<=nbins; xbin++) {
    int tmpRun = int(hEffRB1inVsRun->GetBinCenter(xbin));
    float tmpPressure = pressureByRun[tmpRun];

    float tmpEff =  hEffRB1inVsRun->GetBinContent(xbin);
    if (tmpEff>0) cout << tmpRun << " " << tmpPressure <<endl;
    if (tmpEff>0) hEffVsPressure_RB1in->Fill(tmpPressure,tmpEff);

    tmpEff =  hEffRB1outVsRun->GetBinContent(xbin);
    if (tmpEff>0) hEffVsPressure_RB1out->Fill(tmpPressure,tmpEff);

    tmpEff =  hEffRB2inVsRun->GetBinContent(xbin);
    if (tmpEff>0) hEffVsPressure_RB2in->Fill(tmpPressure,tmpEff);

    tmpEff =  hEffRB2outVsRun->GetBinContent(xbin);
    if (tmpEff>0) hEffVsPressure_RB2out->Fill(tmpPressure,tmpEff);

    tmpEff =  hEffRB3VsRun->GetBinContent(xbin);
    if (tmpEff>0) hEffVsPressure_RB3->Fill(tmpPressure,tmpEff);

    tmpEff =  hEffRB4VsRun->GetBinContent(xbin);
    if (tmpEff>0) hEffVsPressure_RB4->Fill(tmpPressure,tmpEff);

    tmpEff =  hEffREm1VsRun->GetBinContent(xbin);
    if (tmpEff>0 && tmpRun<146807) hEffVsPressure_REm1->Fill(tmpPressure,tmpEff);
    tmpEff =  hEffREm2VsRun->GetBinContent(xbin);
    if (tmpEff>0 && tmpRun<146807) hEffVsPressure_REm2->Fill(tmpPressure,tmpEff);
    tmpEff =  hEffREm3VsRun->GetBinContent(xbin);
    if (tmpEff>0 && tmpRun<146807) hEffVsPressure_REm3->Fill(tmpPressure,tmpEff);
    tmpEff =  hEffREp1VsRun->GetBinContent(xbin);
    if (tmpEff>0 && tmpRun<146807) hEffVsPressure_REp1->Fill(tmpPressure,tmpEff);
    tmpEff =  hEffREp2VsRun->GetBinContent(xbin);
    if (tmpEff>0 && tmpRun<146807) hEffVsPressure_REp2->Fill(tmpPressure,tmpEff);
    tmpEff =  hEffREp3VsRun->GetBinContent(xbin);
    if (tmpEff>0 && tmpRun<146807) hEffVsPressure_REp3->Fill(tmpPressure,tmpEff);

  } 



  // Write Histos
  TDpositionStudies->cd();

  hrecoveryBarrelXYDist->Write();
  hrecoveryBarrelZDist->Write();
  hrecoveryBarrelZvsXYDist->Write();
  hrecoveryEcapDeltaPhi->Write();
  hrecoveryEcapDeltaR->Write();
  hrecoveryEcapDeltaRvsDeltaPhi->Write();
  hEcapDeltaPhi->Write();
  hEcapDeltaR->Write();
  hEcapDeltaPhivsPhi_REm1_R2->Write();
  hEcapDeltaPhivsPhi_REm1_R3->Write();
  hEcapDeltaPhivsPhi_REm2_R2->Write();
  hEcapDeltaPhivsPhi_REm2_R3->Write();
  hEcapDeltaPhivsPhi_REm3_R2->Write();
  hEcapDeltaPhivsPhi_REm3_R3->Write();
  hEcapDeltaPhivsPhi_REp1_R2->Write();
  hEcapDeltaPhivsPhi_REp1_R3->Write();
  hEcapDeltaPhivsPhi_REp2_R2->Write();
  hEcapDeltaPhivsPhi_REp2_R3->Write();
  hEcapDeltaPhivsPhi_REp3_R2->Write();
  hEcapDeltaPhivsPhi_REp3_R3->Write();

  hEcapDeltaPhivsLocalx_REm1_R2->Write();
  hEcapDeltaPhivsLocalx_REm1_R3->Write();
  hEcapDeltaPhivsLocalx_REm2_R2->Write();
  hEcapDeltaPhivsLocalx_REm2_R3->Write();
  hEcapDeltaPhivsLocalx_REm3_R2->Write();
  hEcapDeltaPhivsLocalx_REm3_R3->Write();
  hEcapDeltaPhivsLocalx_REp1_R2->Write();
  hEcapDeltaPhivsLocalx_REp1_R3->Write();
  hEcapDeltaPhivsLocalx_REp2_R2->Write();
  hEcapDeltaPhivsLocalx_REp2_R3->Write();
  hEcapDeltaPhivsLocalx_REp3_R2->Write();
  hEcapDeltaPhivsLocalx_REp3_R3->Write();
  
  hEcapDeltaPhivsLocaly_REm1_R2->Write();
  hEcapDeltaPhivsLocaly_REm1_R3->Write();
  hEcapDeltaPhivsLocaly_REm2_R2->Write();
  hEcapDeltaPhivsLocaly_REm2_R3->Write();
  hEcapDeltaPhivsLocaly_REm3_R2->Write();
  hEcapDeltaPhivsLocaly_REm3_R3->Write();
  hEcapDeltaPhivsLocaly_REp1_R2->Write();
  hEcapDeltaPhivsLocaly_REp1_R3->Write();
  hEcapDeltaPhivsLocaly_REp2_R2->Write();
  hEcapDeltaPhivsLocaly_REp2_R3->Write();
  hEcapDeltaPhivsLocaly_REp3_R2->Write();
  hEcapDeltaPhivsLocaly_REp3_R3->Write();
  
  TDgeneral->cd();
  hpressureByRun->Write();
  hEffVsPressure_RB1in->Write();
  hEffVsPressure_RB1out->Write();
  hEffVsPressure_RB2in->Write();
  hEffVsPressure_RB2out->Write();
  hEffVsPressure_RB3->Write();
  hEffVsPressure_RB4->Write();
  hEffVsPressure_REm1->Write();
  hEffVsPressure_REm2->Write();
  hEffVsPressure_REm3->Write();
  hEffVsPressure_REp1->Write();
  hEffVsPressure_REp2->Write();
  hEffVsPressure_REp3->Write();

  hBarrelEfficiency->Write();
  hRB1inEfficiency->Write();
  hRB1outEfficiency->Write();
  hRB2inEfficiency->Write();
  hRB2outEfficiency->Write();
  hRB3Efficiency->Write();
  hRB4Efficiency->Write();
  hEndcapEfficiency->Write();
  hREm3Efficiency->Write();
  hREm2Efficiency->Write();
  hREm1Efficiency->Write();
  hREp3Efficiency->Write();
  hREp2Efficiency->Write();
  hREp1Efficiency->Write();

  hExpRB1inVsRun->Write();
  hExpRB1outVsRun->Write();
  hExpRB2inVsRun->Write();
  hExpRB2outVsRun->Write();
  hExpRB3VsRun->Write();
  hExpRB4VsRun->Write();
  hFoundRB1inVsRun->Write();
  hFoundRB1outVsRun->Write();
  hFoundRB2inVsRun->Write();
  hFoundRB2outVsRun->Write();
  hFoundRB3VsRun->Write();
  hFoundRB4VsRun->Write();
  hEffRB1inVsRun->Write();
  hEffRB1outVsRun->Write();
  hEffRB2inVsRun->Write();
  hEffRB2outVsRun->Write();
  hEffRB3VsRun->Write();
  hEffRB4VsRun->Write();

  hExpREm3VsRun->Write();
  hExpREm2VsRun->Write();
  hExpREm1VsRun->Write();
  hExpREp3VsRun->Write();
  hExpREp2VsRun->Write();
  hExpREp1VsRun->Write();
  hFoundREm3VsRun->Write();
  hFoundREm2VsRun->Write();
  hFoundREm1VsRun->Write();
  hFoundREp3VsRun->Write();
  hFoundREp2VsRun->Write();
  hFoundREp1VsRun->Write();
  hEffREm3VsRun->Write();
  hEffREm2VsRun->Write();
  hEffREm1VsRun->Write();
  hEffREp3VsRun->Write();
  hEffREp2VsRun->Write();
  hEffREp1VsRun->Write();

  h2DEffREm1->Write();
  h2DEffREm2->Write();
  h2DEffREm3->Write();
  h2DEffREp1->Write();
  h2DEffREp2->Write();
  h2DEffREp3->Write();
  h2DExpectedREm1->Write();
  h2DExpectedREm2->Write();
  h2DExpectedREm3->Write();
  h2DExpectedREp1->Write();
  h2DExpectedREp2->Write();
  h2DExpectedREp3->Write();

  outfile->Close();
}

