// compile with : g++ writeLHEToTree.cc  -o writeLHEToTree.exe `root-config --cflags --glibs`
// last one includes root config file, that basically contains where to find the root libraries

#include "LHEF.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/root/6.02.00-odfocd5/include/TTree.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/root/6.02.00-odfocd5/include/TFile.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/root/6.02.00-odfocd5/include/TMath.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/root/6.02.00-odfocd5/include/TLorentzVector.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <vector>

// function takes cmsgrid_final.lhe file for G->WW->lvjj and G->WW/ZZ->jjjj and writes important properties to a root tree
//writes branches with:
// neve : event number
// mWW : invariant mass of WW or ZZ
// Wlep/Whad pt,eta,phi,e
// W+/W- pt,eta,phi,e for all hadronic channel
// for semileptonic channel:
// lepton pt,eta,phi,e
// tauon pt,eta,phi,e if there is a tau in the event
// leptonID = 13 or 11 depending on lepton in final state
// isTauEvent: True if event had tau False otherwise
// neutrino: pt,eta,phi,e
// important : pdgId of G =39 at the moment


int main(int argc, char ** argv)
{
  
  if(argc<3)
  {
    std::cout << " too few input parameters \n program needs: outputdirectory/output_filename decayMode = WW_lvjj or WW_jjjj or ZZ_jjjj input_directory  " <<std::endl;
    return 1;
  }
  std::string Name = argv[1];
  std::string decayMode = argv[2];
  
  
  std::string dir = "/home/dschaefer/root/genMC/"; 
  if(argc==4)
  {
   dir = argv[3];   
  }
  dir = dir+Name+"/";
  std::string outdir = "/storage/jbod/dschaefer/";
  std::string filename = "cmsgrid_final.lhe";
  std::cout << dir+filename << std::endl;
  LHEF::Reader reader(dir+filename);
  
  int pdgIdGraviton = 39;
  
  long neve = 0;
  
  TFile* output = new TFile((outdir+Name+".root").c_str(),"RECREATE");
  TTree* tree = new TTree("tree","tree");
  
  if(decayMode.find("lvjj")!=std::string::npos)
  {
      int NumberOfTauEvents =0;
      int NumberOfTauToEleEvents =0;
      int NumberOfTauToMuEvents =0;
      
  float Wlep_px;
  float Wlep_py;
  float Wlep_pz;
  float Wlep_E ;
  float Wlep_M ;
  
  float Whad_px ;
  float Whad_py ;
  float Whad_pz ;
  float Whad_E  ;
  float Whad_M  ;
  
  float Wlep_pt ;
  float Wlep_eta;
  float Wlep_phi;
  
  float Whad_pt ;
  float Whad_eta;
  float Whad_phi;
  int   leptonID;
  bool  isTauEvent;
  
  float lepton_px;
  float lepton_py;
  float lepton_pz;
  float lepton_E ;
  float lepton_M ;
  float lepton_pt ;
  float lepton_eta;
  float lepton_phi;

  float neutrino_E ;
  float neutrino_M ;
  float neutrino_pt ;
  float neutrino_eta;
  float neutrino_phi;
  
  float tau_pt;
  float tau_eta;
  float tau_phi;
  float tau_E;
  
  
  
  float mWW;
  
    tree->Branch("neve",&neve,"neve/L");
    tree->Branch("mWW" ,&mWW ,"mWW/F" );
    tree->Branch("leptonID",&leptonID,"leptonID/I");
    tree->Branch("isTauEvent",&isTauEvent,"isTauEvent/O");
    tree->Branch("Wlep_pt", &Wlep_pt, "Wlep_pt/F");
    tree->Branch("Wlep_eta",&Wlep_eta,"Wlep_eta/F");
    tree->Branch("Wlep_phi",&Wlep_phi,"Wlep_phi/F");
    tree->Branch("Wlep_E",  &Wlep_E,  "Wlep_E/F");

    tree->Branch("lepton_pt", &lepton_pt, "lepton_pt/F");
    tree->Branch("lepton_eta",&lepton_eta,"lepton_eta/F");
    tree->Branch("lepton_phi",&lepton_phi,"lepton_phi/F");
    tree->Branch("lepton_E",  &lepton_E,  "lepton_E/F");
    
    tree->Branch("neutrino_pt", &neutrino_pt, "neutrino_pt/F");
    tree->Branch("neutrino_eta",&neutrino_eta,"neutrino_eta/F");
    tree->Branch("neutrino_phi",&neutrino_phi,"neutrino_phi/F");
    tree->Branch("neutrino_E",  &neutrino_E,  "neutrino_E/F");
    
    tree->Branch("Whad_pt", &Whad_pt, "Whad_pt/F");
    tree->Branch("Whad_eta",&Whad_eta,"Whad_eta/F");
    tree->Branch("Whad_phi",&Whad_phi,"Whad_phi/F");
    tree->Branch("Whad_E",  &Whad_E  ,"Whad_E/F");
    
    tree->Branch("tau_pt",&tau_pt,"tau_pt/F");
    tree->Branch("tau_eta",&tau_eta,"tau_eta/F");
    tree->Branch("tau_phi",&tau_phi,"tau_phi/F");
    tree->Branch("tau_E",&tau_E,"tau_E/F");
    
    int index_Wlep =-1;
    int index_Whad =-1;
    int i_q1 =-1;
    int i_q2 =-1;
    int i_l1 =-1;
    int i_l2 =-1;
    TLorentzVector l;
    int whichLepton=0;

  
  while ( reader.readEvent() ) {
    ++neve;
    if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
    LHEF::HEPEUP eventInfo = reader.hepeup;
    
    int N = eventInfo.NUP;
    
    std::vector<long> pdgId = eventInfo.IDUP;
    std::vector< std::vector<double> > list_P4 = eventInfo.PUP;
    std::vector<std::pair<int,int> > mothers = eventInfo.MOTHUP;
    std::vector<double> P4_had;
    std::vector<double> P4_lep;
    std::vector<double> P4_lepton;
    std::vector<double> P4_neutrino;
    std::vector<double> P4_tau;
    TLorentzVector L4_neutrino(0.,0.,0.,0.);
    
    
    i_q1 =-1;
    i_q2 =-1;
    i_l1 =-1;
    i_l2 =-1;
    whichLepton=0;
    isTauEvent=0;
    
    
    for(int i=0;i<N;i++)
    {
    
      std::pair<int,int> pair = mothers.at(i);
      int first = pair.first-1;
      int second= pair.second-1;
      if(first<0 or second<0) continue;
      
     
      
      if(TMath::Abs(pdgId.at(i))==15)
      {
       P4_tau = list_P4.at(i);   
      }
      if(TMath::Abs(pdgId.at(i))==11 or  TMath::Abs(pdgId.at(i))==13)
      {
	P4_lepton = list_P4.at(i);
      }
      if(TMath::Abs(pdgId.at(i))==12 or TMath::Abs(pdgId.at(i))==14 or TMath::Abs(pdgId.at(i))==16)
      {
       P4_neutrino = list_P4.at(i);
       TLorentzVector tmp;
       tmp.SetPxPyPzE(P4_neutrino.at(0),P4_neutrino.at(1),P4_neutrino.at(2),P4_neutrino.at(3));
       L4_neutrino = L4_neutrino + tmp;
      }
      
      
      if(TMath::Abs(pdgId.at(i)) <=8)
      {
	
	if(TMath::Abs(pdgId.at(first))==24 or TMath::Abs(pdgId.at(first))==23)
	{
	  P4_had = list_P4.at(first);
	  index_Whad = first;
	}
	if(TMath::Abs(pdgId.at(second))==24 or TMath::Abs(pdgId.at(second))==23)
	{
	  P4_had = list_P4.at(second); 
	  index_Whad = second;
	}
	if(TMath::Abs(pdgId.at(first))==pdgIdGraviton and TMath::Abs(pdgId.at(second))==pdgIdGraviton )
	{
	  if(i_q1==-1)
	  {
	    i_q1 = i;
	  }
	  else
	  {
	    i_q2 =i; 
	  }
	}
      }
      
      if(TMath::Abs(pdgId.at(i)) >10 and TMath::Abs(pdgId.at(i))<17)
      {
	
	if(TMath::Abs(pdgId.at(i))==11){whichLepton=11;}
	if(TMath::Abs(pdgId.at(i))==13){whichLepton=13;}
	if(TMath::Abs(pdgId.at(i))==15){isTauEvent=1; NumberOfTauEvents+=1;}
	if(TMath::Abs(pdgId.at(first))==24)
	{
	  P4_lep = list_P4.at(first);
	  index_Wlep = first;
	}
	if(TMath::Abs(pdgId.at(second))==24)
	{
	  P4_lep = list_P4.at(second); 
	  index_Wlep = second;
	}
	if(TMath::Abs(pdgId.at(first))==pdgIdGraviton and TMath::Abs(pdgId.at(second))==pdgIdGraviton )
	{
	  if(i_l1==-1)
	  {
	    i_l1 = i;
	  }
	  else
	  {
	    i_l2 =i; 
	  }
	}
      }
      if(i_q1!=-1 and i_q2!=-1)
      {
	TLorentzVector q1;
	TLorentzVector q2;
	TLorentzVector W;
	q1.SetPxPyPzE(list_P4.at(i_q1).at(0),list_P4.at(i_q1).at(1),list_P4.at(i_q1).at(2),list_P4.at(i_q1).at(3));
	q2.SetPxPyPzE(list_P4.at(i_q2).at(0),list_P4.at(i_q2).at(1),list_P4.at(i_q2).at(2),list_P4.at(i_q2).at(3));
	
	W = q1+q2;
	P4_had.push_back(W.Px());
	P4_had.push_back(W.Py());
	P4_had.push_back(W.Pz());
	P4_had.push_back(W.E());
	P4_had.push_back(W.M());
	
      }
      if(i_l1!=-1 and i_l2!=-1)
      {
	TLorentzVector l1;
	TLorentzVector l2;
	TLorentzVector W;
	l1.SetPxPyPzE(list_P4.at(i_l1).at(0),list_P4.at(i_l1).at(1),list_P4.at(i_l1).at(2),list_P4.at(i_l1).at(3));
	l2.SetPxPyPzE(list_P4.at(i_l2).at(0),list_P4.at(i_l2).at(1),list_P4.at(i_l2).at(2),list_P4.at(i_l2).at(3));
	
	W = l1+l2;
	P4_lep.push_back(W.Px());
	P4_lep.push_back(W.Py());
	P4_lep.push_back(W.Pz());
	P4_lep.push_back(W.E());
	P4_lep.push_back(W.M());
	
      }
    }
    
    if(isTauEvent)
    {
     if(leptonID==11)
     {
      NumberOfTauToEleEvents+=1;   
     }
     if(leptonID==13)
     {
      NumberOfTauToMuEvents+=1;   
     }
    }
    
    
    
    if(isTauEvent)
    {
     TLorentzVector tau; tau.SetPxPyPzE(P4_tau.at(0),P4_tau.at(1),P4_tau.at(2),P4_tau.at(3));
     tau_pt = tau.Pt();
     tau_eta = tau.Eta();
     tau_phi = tau.Phi();
     tau_E   = tau.E();
    }
    else
    {
     tau_pt  = 0;
     tau_eta = 0;
     tau_phi = 0;
     tau_E   = 0;
    }
    
    
    TLorentzVector Whad;
    Whad_px = P4_had.at(0);
    Whad_py = P4_had.at(1);
    Whad_pz = P4_had.at(2);
    Whad_E  = P4_had.at(3);
    Whad_M  = P4_had.at(4);
    Whad.SetPxPyPzE(Whad_px,Whad_py,Whad_pz,Whad_E);
    Whad_pt  = Whad.Pt();
    Whad_eta = Whad.Eta();
    Whad_phi = Whad.Phi();
    Whad_E   = Whad.E();
    
    neutrino_E = L4_neutrino.E();
    neutrino_eta = L4_neutrino.Eta();
    neutrino_phi = L4_neutrino.Phi();
    neutrino_pt  = L4_neutrino.Pt();
    
    
    TLorentzVector Wlep;
    Wlep_px = P4_lep.at(0);
    Wlep_py = P4_lep.at(1);
    Wlep_pz = P4_lep.at(2);
    Wlep_E  = P4_lep.at(3);
    Wlep_M  = P4_lep.at(4);
    Wlep.SetPxPyPzE(Wlep_px,Wlep_py,Wlep_pz,Wlep_E);
    Wlep_pt  =Wlep.Pt();
    Wlep_eta =Wlep.Eta();
    Wlep_phi =Wlep.Phi();
    Wlep_E   =Wlep.E();
    
    lepton_px = P4_lepton.at(0);
    lepton_py = P4_lepton.at(1);
    lepton_pz = P4_lepton.at(2);
    lepton_E  = P4_lepton.at(3);
    lepton_M  = P4_lepton.at(4);
    l.SetPxPyPzE(lepton_px,lepton_py,lepton_pz,lepton_E);
    lepton_pt  = l.Pt();
    lepton_eta = l.Eta();
    lepton_phi = l.Phi();
    lepton_E   = l.E();
    
    leptonID = whichLepton;
    TLorentzVector G = Wlep+Whad;
    mWW = G.M();

    tree->Fill();
  }
  
  tree->Write();
  output->Close();
  return 0;
  }
  if(decayMode.find("jjjj")!=std::string::npos and (decayMode.find("WW")!=std::string::npos or decayMode.find("WZ")!=std::string::npos))
  {
    float Wplus_px;
    float Wplus_py;
    float Wplus_pz;
    float Wplus_E ;
    float Wplus_M ;
  
    float Wminus_px ;
    float Wminus_py ;
    float Wminus_pz ;
    float Wminus_E  ;
    float Wminus_M  ;
  
    float Wplus_pt ;
    float Wplus_eta;
    float Wplus_phi;
  
    float Wminus_pt ;
    float Wminus_eta;
    float Wminus_phi;
    float mWW;
  
    tree->Branch("neve",&neve,"neve/L");
    tree->Branch("mWW",&mWW,"mWW/F");
    tree->Branch("Wplus_pt", &Wplus_pt, "Wplus_pt/F");
    tree->Branch("Wplus_eta",&Wplus_eta,"Wplus_eta/F");
    tree->Branch("Wplus_phi",&Wplus_phi,"Wplus_phi/F");
    tree->Branch("Wplus_E",  &Wplus_E,  "Wplus_E/F");

    tree->Branch("Wminus_pt", &Wminus_pt, "Wminus_pt/F");
    tree->Branch("Wminus_eta",&Wminus_eta,"Wminus_eta/F");
    tree->Branch("Wminus_phi",&Wminus_phi,"Wminus_phi/F");
    tree->Branch("Wminus_E",  &Wminus_E  ,"Wminus_E/F");
    
    int index_Wplus =-1;
    int index_Wminus =-1;
    int i_q1 =-1;
    int i_q2 =-1;
    int i_q3 =-1;
    int i_q4 =-1;
    TLorentzVector l;

  
  while ( reader.readEvent() ) {
    ++neve;
    if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
    LHEF::HEPEUP eventInfo = reader.hepeup;
 
    int N = eventInfo.NUP;
    std::vector<long> pdgId = eventInfo.IDUP;
    std::vector< std::vector<double> > list_P4 = eventInfo.PUP;
    std::vector<std::pair<int,int> > mothers = eventInfo.MOTHUP;
    std::vector<double> P4_minus;
    std::vector<double> P4_plus;

    i_q1 =-1;
    i_q2 =-1;
    i_q3 =-1;
    i_q4 =-1;
    
    for(int i=0;i<N;i++)
    {
      std::pair<int,int> pair = mothers.at(i);
      int first = pair.first-1;
      int second= pair.second-1;
      if(first<0 or second<0) continue;
      
      if(TMath::Abs(pdgId.at(i)) <=8)
      {
        if(decayMode.find("Z")!=std::string::npos)
        {
        if(TMath::Abs(pdgId.at(first))==23)
	{
	  P4_plus = list_P4.at(first);
	  index_Wplus = first;
	}
	if(TMath::Abs(pdgId.at(second))==23)
	{
	  P4_plus = list_P4.at(second); 
	  index_Wplus = second;
	} 
        }
        else
        {
	if(pdgId.at(first)==24)
	{
	  P4_plus = list_P4.at(first);
	  index_Wplus = first;
	}
	if(pdgId.at(second)==24)
	{
	  P4_plus = list_P4.at(second); 
	  index_Wplus = second;
	}
        }
	if(TMath::Abs(pdgId.at(first))==pdgIdGraviton and TMath::Abs(pdgId.at(second))==pdgIdGraviton )
	{
	  if(i_q1==-1)
	  {
	    i_q1 = i;
	  }
	  else
	  {
	    i_q2 =i; 
	  }
	}
      }
      
      if(TMath::Abs(pdgId.at(i)) <=8)
      {
        if(decayMode.find("Z")!=std::string::npos)
        {
        if(TMath::Abs(pdgId.at(first))==24)
	{
	  P4_plus = list_P4.at(first);
	  index_Wplus = first;
	}
	if(TMath::Abs(pdgId.at(second))==24)
	{
	  P4_plus = list_P4.at(second); 
	  index_Wplus = second;
	} 
        }
        else
        {
	if(pdgId.at(first)==-24)
	{
	  P4_minus = list_P4.at(first);
	  index_Wminus = first;
	}
	if(pdgId.at(second)==-24)
	{
	  P4_minus = list_P4.at(second); 
	  index_Wminus = second;
	}
        }
	if(TMath::Abs(pdgId.at(first))==pdgIdGraviton and TMath::Abs(pdgId.at(second))==pdgIdGraviton )
	{
	  if(i_q3==-1)
	  {
	    i_q3 = i;
	  }
	  else
	  {
	    i_q4 =i; 
	  }
	}
      }
      if(i_q1!=-1 and i_q2!=-1)
      {
	TLorentzVector q1;
	TLorentzVector q2;
	TLorentzVector W;
	q1.SetPxPyPzE(list_P4.at(i_q1).at(0),list_P4.at(i_q1).at(1),list_P4.at(i_q1).at(2),list_P4.at(i_q1).at(3));
	q2.SetPxPyPzE(list_P4.at(i_q2).at(0),list_P4.at(i_q2).at(1),list_P4.at(i_q2).at(2),list_P4.at(i_q2).at(3));
	
	W = q1+q2;
	P4_plus.push_back(W.Px());
	P4_plus.push_back(W.Py());
	P4_plus.push_back(W.Pz());
	P4_plus.push_back(W.E());
	P4_plus.push_back(W.M());
	
      }
      if(i_q3!=-1 and i_q4!=-1)
      {
	TLorentzVector q3;
	TLorentzVector q4;
	TLorentzVector W;
	q3.SetPxPyPzE(list_P4.at(i_q3).at(0),list_P4.at(i_q3).at(1),list_P4.at(i_q3).at(2),list_P4.at(i_q3).at(3));
	q4.SetPxPyPzE(list_P4.at(i_q4).at(0),list_P4.at(i_q4).at(1),list_P4.at(i_q4).at(2),list_P4.at(i_q4).at(3));
	
	W = q3+q4;
	P4_minus.push_back(W.Px());
	P4_minus.push_back(W.Py());
	P4_minus.push_back(W.Pz());
	P4_minus.push_back(W.E());
	P4_minus.push_back(W.M());
	
      }
    }
    
    Wplus_px = P4_plus.at(0);
    Wplus_py = P4_plus.at(1);
    Wplus_pz = P4_plus.at(2);
    Wplus_E  = P4_plus.at(3);
    Wplus_M  = P4_plus.at(4);
    l.SetPxPyPzE(Wplus_px,Wplus_py,Wplus_pz,Wplus_E);
    Wplus_pt  = l.Pt();
    Wplus_eta = l.Eta();
    Wplus_phi = l.Phi();
    Wplus_E   = l.E();
    
 
    TLorentzVector Wminus;
    Wminus_px = P4_minus.at(0);
    Wminus_py = P4_minus.at(1);
    Wminus_pz = P4_minus.at(2);
    Wminus_E  = P4_minus.at(3);
    Wminus_M  = P4_minus.at(4);
    Wminus.SetPxPyPzE(Wminus_px,Wminus_py,Wminus_pz,Wminus_E);
    Wminus_pt  = Wminus.Pt();
    Wminus_eta = Wminus.Eta();
    Wminus_phi = Wminus.Phi();
    Wminus_E   = Wminus.E();
    
    TLorentzVector G = l+Wminus;
    mWW = G.M();

    tree->Fill();
  }
  tree->Write();
  output->Close();
    
  return 0; 
  }
  if(decayMode.find("jjjj")!=std::string::npos and decayMode.find("ZZ")!=std::string::npos)
  {
    float Z1_px;
    float Z1_py;
    float Z1_pz;
    float Z1_E ;
    float Z1_M ;
  
    float Z2_px ;
    float Z2_py ;
    float Z2_pz ;
    float Z2_E  ;
    float Z2_M  ;
  
    float Z1_pt ;
    float Z1_eta;
    float Z1_phi;
  
    float Z2_pt ;
    float Z2_eta;
    float Z2_phi;
    float mWW;
  
    tree->Branch("neve",&neve,"neve/L");
    tree->Branch("mWW",&mWW,"mWW/F");
    tree->Branch("Z1_pt", &Z1_pt, "Z1_pt/F");
    tree->Branch("Z1_eta",&Z1_eta,"Z1_eta/F");
    tree->Branch("Z1_phi",&Z1_phi,"Z1_phi/F");
    tree->Branch("Z1_E",  &Z1_E,  "Z1_E/F");

    tree->Branch("Z2_pt", &Z2_pt, "Z2_pt/F");
    tree->Branch("Z2_eta",&Z2_eta,"Z2_eta/F");
    tree->Branch("Z2_phi",&Z2_phi,"Z2_phi/F");
    tree->Branch("Z2_E",  &Z2_E  ,"Z2_E/F");
    
    int index_Z1 =-1;
    int index_Z2 =-1;
    int i_q1 =-1;
    int i_q2 =-1;
    int i_q3 =-1;
    int i_q4 =-1;
    TLorentzVector l;

  
  while ( reader.readEvent() ) {
    ++neve;
    if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
    LHEF::HEPEUP eventInfo = reader.hepeup;
 
    int N = eventInfo.NUP;
    std::vector<long> pdgId = eventInfo.IDUP;
    std::vector< std::vector<double> > list_P4 = eventInfo.PUP;
    std::vector<std::pair<int,int> > mothers = eventInfo.MOTHUP;
    std::vector<double> P4_Z2;
    std::vector<double> P4_Z1;
    std::vector<int> Zcandidate;

    i_q1 =-1;
    i_q2 =-1;
    i_q3 =-1;
    i_q4 =-1;
    
    std::cout  << N << std::endl;
    for(int i=0;i<N;i++)
    {
      std::pair<int,int> pair = mothers.at(i);
      int first = pair.first-1;
      int second= pair.second-1;
      if(first<0 or second<0) continue;
     
      if(TMath::Abs(pdgId.at(i)) <=8)
      {

	if(TMath::Abs(pdgId.at(first))==23)
	{
	  index_Z1 = first;
          Zcandidate.push_back(index_Z1);
        }
	if(TMath::Abs(pdgId.at(second))==23)
	{
            index_Z1 = second;
            Zcandidate.push_back(index_Z1);
        }
	if(TMath::Abs(pdgId.at(first))==pdgIdGraviton and TMath::Abs(pdgId.at(second))==pdgIdGraviton)
	{
          if(i_q1==-1)
	  {
	    i_q1 = i;
	  }
	  else
	  {
	    i_q2 =i; 
	  }
	}
      }
      
      if(TMath::Abs(pdgId.at(i)) <=8)
      {
	
	if(TMath::Abs(pdgId.at(first))==23)
	{
	  index_Z2 = first;
          Zcandidate.push_back(index_Z2);
	}
	if(TMath::Abs(pdgId.at(second))==23)
	{
	  index_Z2 = second;
          Zcandidate.push_back(index_Z2);
	}
	if(TMath::Abs(pdgId.at(first))==pdgIdGraviton and TMath::Abs(pdgId.at(second))==pdgIdGraviton)
	{
	  if(i_q3==-1)
	  {
	    i_q3 = i;
	  }
	  else
	  {
	    i_q4 =i; 
	  }
	}
      }
    }
      if(i_q1!=-1 and i_q2!=-1)
      {
	TLorentzVector q1;
	TLorentzVector q2;
	TLorentzVector W;
	q1.SetPxPyPzE(list_P4.at(i_q1).at(0),list_P4.at(i_q1).at(1),list_P4.at(i_q1).at(2),list_P4.at(i_q1).at(3));
	q2.SetPxPyPzE(list_P4.at(i_q2).at(0),list_P4.at(i_q2).at(1),list_P4.at(i_q2).at(2),list_P4.at(i_q2).at(3));
	
	W = q1+q2;
	P4_Z1.push_back(W.Px());
	P4_Z1.push_back(W.Py());
	P4_Z1.push_back(W.Pz());
	P4_Z1.push_back(W.E());
	P4_Z1.push_back(W.M());
	
      }
      if(i_q3!=-1 and i_q4!=-1)
      {
	TLorentzVector q3;
	TLorentzVector q4;
	TLorentzVector W;
	q3.SetPxPyPzE(list_P4.at(i_q3).at(0),list_P4.at(i_q3).at(1),list_P4.at(i_q3).at(2),list_P4.at(i_q3).at(3));
	q4.SetPxPyPzE(list_P4.at(i_q4).at(0),list_P4.at(i_q4).at(1),list_P4.at(i_q4).at(2),list_P4.at(i_q4).at(3));
	
	W = q3+q4;
	P4_Z2.push_back(W.Px());
	P4_Z2.push_back(W.Py());
	P4_Z2.push_back(W.Pz());
	P4_Z2.push_back(W.E());
	P4_Z2.push_back(W.M());
	
      }
    if(Zcandidate.size()>0)
    {

     if(Zcandidate.size()==16 and Zcandidate.at(0)!=Zcandidate.at(15))
     {
        P4_Z2 = list_P4.at(Zcandidate.at(0));
        P4_Z1 = list_P4.at(Zcandidate.at(15));
        if(P4_Z1.at(0)<P4_Z2.at(0))
        {
            P4_Z2 = list_P4.at(Zcandidate.at(15));
            P4_Z1 = list_P4.at(Zcandidate.at(0));
        }
     }
     if((Zcandidate.size()==8 and Zcandidate.at(0)==Zcandidate.at(7)))
     {
       if(i_q1==-1)
       {
        P4_Z1 = list_P4.at(Zcandidate.at(0));   
       }
       else if(i_q3==-1)
       {
        P4_Z2 = list_P4.at(Zcandidate.at(0));   
       }
     }
     }
      
      
    Z1_px = P4_Z1.at(0);
    Z1_py = P4_Z1.at(1);
    Z1_pz = P4_Z1.at(2);
    Z1_E  = P4_Z1.at(3);
    Z1_M  = P4_Z1.at(4);
    l.SetPxPyPzE(Z1_px,Z1_py,Z1_pz,Z1_E);
    Z1_pt  = l.Pt();
    Z1_eta = l.Eta();
    Z1_phi = l.Phi();
    Z1_E   = l.E();
    
 
    TLorentzVector Z2;
    Z2_px = P4_Z2.at(0);
    Z2_py = P4_Z2.at(1);
    Z2_pz = P4_Z2.at(2);
    Z2_E  = P4_Z2.at(3);
    Z2_M  = P4_Z2.at(4);
    Z2.SetPxPyPzE(Z2_px,Z2_py,Z2_pz,Z2_E);
    Z2_pt  = Z2.Pt();
    Z2_eta = Z2.Eta();
    Z2_phi = Z2.Phi();
    Z2_E   = Z2.E();
    
    TLorentzVector G = l+Z2;
    mWW = G.M();

    tree->Fill();
  }
  tree->Write();
  output->Close();
   return 0;
  }
  return 1;
}
