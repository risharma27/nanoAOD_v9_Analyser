#define nano9Ana_cxx
// The class definition in nano9Ana.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.



#include "nano9Ana.h"
#include <TH2.h>
#include <TStyle.h>

void nano9Ana::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
}

void nano9Ana::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  //Initialization of the counters:
  nEvtRan        = 0;
  nEvtTotal      = 0;
  //Other custom counters can be initialized here.
  
  _HstFile = new TFile(_HstFileName,"recreate");
  BookHistograms();
}

void nano9Ana::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
   _HstFile->Write();
   _HstFile->Close();

  //The following lines are displayed on the root prompt.
  cout<<"Total events ran = "<<nEvtRan<<endl;
  cout<<"Total good events = "<<nEvtTotal<<endl;

  //The following lines are written on the sum_<process name>.txt file
  ofstream fout(_SumFileName);
  fout<<"Total events ran = "<<nEvtRan<<endl;
  fout<<"Total good events  = "<<nEvtTotal<<endl;
}

void nano9Ana::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file. 
}

Bool_t nano9Ana::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  fReader.SetLocalEntry(entry);
  if(_data == 0)
    fReader_MC.SetLocalEntry(entry);
  if(_data == 1)
    fReader_Data.SetLocalEntry(entry);

  //Verbosity determines the number of processed events after which the root prompt is supposed to display a status update.
  if(_verbosity==0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
  
  //The following flags throws away some events based on unwanted properties (such as detector problems)
  GoodEvt2018 = (_year==2018 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  
   GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;
  
  nEvtRan++;                             //Total number of events containing everything (including the trash events).
  
  if(GoodEvt){
    nEvtTotal++;                         //Total number of events containing goodEvents
                                         //The analysis is done for these good events.


    //Construction of the arrays:
    
    int nlep = 0;                        //This counts the number of leptons (excluding taus) in each event.
    goodLepton.clear();                   
    //goodMu array :
    int nmu = 0;                         // This counts the number of muons in each event.
    goodMu.clear();                      // Make sure to empty the array from previous event.
    for(unsigned int i=0; i<(*nMuon); i++){
                                         // This loop runs over all the muon candidates. Some of them will pass our selection criteria.
                                         // These will be stored in the goodMu array.
      Lepton temp;                       // 'temp' is the i-th candidate.
      temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105); //the muon mass in GeV is 0.105
      temp.id = -13*Muon_charge[i];    //pdgID for mu- = 13, pdgID for mu+ = -13  
      temp.ind = i; 

      //These are the flags the 'temp' object i.e. the muon candidate has to pass.
      bool passCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i];
      //passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
      //passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;
      
      if(passCuts){
	goodMu.push_back(temp); // If 'temp' satisfies all the conditions, it is pushed back into goodMu
	goodLepton.push_back(temp);
	
      }
    }                                  // This 'for' loop has created a goodMu array.
    
    //Now we sort the goodMu in decreasing pT order
    //Sort(0);
    
    Sortpt(goodMu);

    
    //met
 
    float metpt = *MET_pt;
    h.met[0]->Fill(metpt);

    float metphi = *MET_phi;
    h.met[1]->Fill(metphi);

    

    
   //Other arrays, such as RecoEle, GenMu, GenEle can be constructed here.


   //goodElectron array :
    int nelec = 0;
    goodElectron.clear();
    for(unsigned int i=0; i<(*nElectron); i++){

    Lepton temp;
      temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],Electron_mass[i]);
      temp.id = -11*Electron_charge[i];
      temp.ind = i;
      
      bool passCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4;
      
      if(passCuts){
	goodElectron.push_back(temp);
	goodLepton.push_back(temp);
      }
    }

    Sortpt(goodElectron);

    Sortpt(goodLepton);  //Sorts the goodLepton array.
    
    
    //goodPhoton array :
    int nph = 0;
    goodPhoton.clear();
    for(unsigned int i=0; i<(*nPhoton); i++){

    Lepton temp;
      temp.v.SetPtEtaPhiM(Photon_pt[i],Photon_eta[i],Photon_phi[i],Photon_mass[i]);

      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4;
      
      if(passCuts){
	goodPhoton.push_back(temp);
      }
    }

    Sortpt(goodPhoton);



    //goodJet array :
    
    int njet =0;
    goodJet.clear();
    btagged_jet.clear();
    for(unsigned int i=0; i<(*nJet); i++){

    Lepton temp;
      temp.v.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);

      bool passCuts = temp.v.Pt()>30 && fabs(temp.v.Eta())<3.0;
      
      if(passCuts){
	goodJet.push_back(temp);
      }

       if(passCuts && Jet_btagDeepFlavB[i]>0.4184){
	 btagged_jet.push_back(temp);
	 
       }
    }
   
    Sortpt(goodJet);
    Sortpt(btagged_jet);

   
    //genparticle array :
    genMu.clear();
    genElec.clear();
    
    for(unsigned int i=0; i<(*nGenPart); i++){
      Lepton temp; temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],GenPart_mass[i]);
      temp.status = GenPart_status[i]; temp.ind = i; temp.pdgid = GenPart_pdgId[i]; temp.momid=MotherID(i,GenPart_genPartIdxMother[i]);
      
      bool muPasscut = abs(temp.pdgid)==13 && temp.status==1 && temp.v.Pt()>5 && fabs(temp.v.Eta())<2.4;
      if(muPasscut){
	genMu.push_back(temp);
      }
      
      bool elPasscut = abs(temp.pdgid)==11 && temp.status==1 && temp.v.Pt()>5 && fabs(temp.v.Eta())<2.17;
      if(elPasscut){
	genElec.push_back(temp);
      }
    }
    
    Sortpt(genMu);
    Sortpt(genElec);

    
    
    //##############
    // Analysis
    //##############

    //For muons
    
    h.nmuons->Fill((int)goodMu.size());   //Fill the size of goodMu array
    h.nelectrons->Fill((int)goodElectron.size());
    h.nleptons->Fill((int)goodLepton.size());
    h.nphotons->Fill((int)goodPhoton.size());
    h.njets->Fill((int)goodJet.size());
    h.bjets->Fill((int)btagged_jet.size());

   
    //Plot pT, eta, phi of all muons in same plot
    for(int i=0; i<(int)goodMu.size(); i++){
      h.mu[0]->Fill(goodMu.at(i).v.Pt());
      h.mu[1]->Fill(goodMu.at(i).v.Eta());
      h.mu[2]->Fill(goodMu.at(i).v.Phi());
    
    }

    
    //Plotting the leading muon pT in each event.
    if((int)goodMu.size()>0){             
      h.muprop[0]->Fill(goodMu.at(0).v.Pt());
      h.muprop[1]->Fill(goodMu.at(0).v.Eta());
      h.muprop[2]->Fill(goodMu.at(0).v.Phi());
      h.muprop[3]->Fill(Muon_pfRelIso04_all[goodMu.at(0).ind]);
    }


    //Plotting the subleading muon pT and Eta in each event.
    if((int)goodMu.size()>1){
      h.submuprop[0]->Fill(goodMu.at(1).v.Pt());
      h.submuprop[1]->Fill(goodMu.at(1).v.Eta());
      h.submuprop[2]->Fill(goodMu.at(1).v.Phi());
    }
    
    //Plotting dphi, dR and  M_{#mu#mu} between the leading and subleading muon
    if((int)goodMu.size()>1){
    float dphi = goodMu.at(0).v.DeltaPhi(goodMu.at(1).v);
    float dR = goodMu.at(0).v.DeltaR(goodMu.at(1).v);
    h.dimuon[0]->Fill(dphi);
    h.dimuon[1]->Fill(dR);
    float dimu_invmass = (goodMu.at(0).v+goodMu.at(1).v).M();
    h.dimuon[2]->Fill(dimu_invmass);
    }  


    //##########################################################################

    
    //For Electrons

    
    //Plotting the leading electron pT, Eta and Phi in each event.
    if((int)goodElectron.size()>0){             
      h.eprop[0]->Fill(goodElectron.at(0).v.Pt());
      h.eprop[1]->Fill(goodElectron.at(0).v.Eta());
      h.eprop[2]->Fill(goodElectron.at(0).v.Phi());
      
    }
    


    //Plotting the subleading electron pT, Eta and Phi in each event.
    if((int)goodElectron.size()>1){
      h.subeprop[0]->Fill(goodElectron.at(1).v.Pt());
      h.subeprop[1]->Fill(goodElectron.at(1).v.Eta());
      h.subeprop[2]->Fill(goodElectron.at(1).v.Phi());

    }
    

    //Plotting dphi, dR and  M_{#e#e} between the leading and subleading electron
    if((int)goodElectron.size()>1){
      float edphi = goodElectron.at(0).v.DeltaPhi(goodElectron.at(1).v);
      float edR = goodElectron.at(0).v.DeltaR(goodElectron.at(1).v);
      h.dielectron[0]->Fill(edphi);
      h.dielectron[1]->Fill(edR);
      float diel_invmass = (goodElectron.at(0).v+goodElectron.at(1).v).M();
      h.dielectron[2]->Fill(diel_invmass);
    } 
    
    

     //#######################################################################



    //Plotting transverse mass

    //For dimuon events

    if(goodMu.size()>1){
      float dimu_phi = delta_phi(float (metphi), float (goodMu.at(0).v.Phi()));
      float dimuon_dphi = dimu_phi;
      float dimuon_mT = transv_mass(float (goodMu.at(0).v.Pt()), float (metpt), float (dimuon_dphi));
      h.transvmass[0]->Fill(dimuon_mT);
    }

    //For dielectron events
    
    if(goodElectron.size()>1){
      float diel_phi = delta_phi(float (metphi), float (goodElectron.at(0).v.Phi()));
      float diel_dphi = diel_phi;
      float dielec_mT = transv_mass(float (goodElectron.at(0).v.Pt()), float (metpt), float (diel_dphi));
      h.transvmass[1]->Fill(dielec_mT);
    }
    
    //For e#mu events and leading lepton is muon

    if(goodMu.size()>0 && goodElectron.size()>0){
      float emu_phi = delta_phi(float (metphi), float (goodLepton.at(0).v.Phi()));
      float emu_dphi = emu_phi;
      float emu_mT = transv_mass(float (goodLepton.at(0).v.Pt()), float (metpt), float (emu_dphi));
      h.transvmass[2]->Fill(emu_mT);
    }

    

    //#######################################################################
    

    //For Photons


    
    //Plotting the leading photon pT, Eta and Phi in each event.
    if((int)goodPhoton.size()>0){             
      h.ph[0]->Fill(goodPhoton.at(0).v.Pt());
      h.ph[1]->Fill(goodPhoton.at(0).v.Eta());
      h.ph[2]->Fill(goodPhoton.at(0).v.Phi());
      
    }
    


    //Plotting the subleading photon pT, Eta and Phi in each event.
    if((int)goodPhoton.size()>1){
      h.subph[0]->Fill(goodPhoton.at(1).v.Pt());
      h.subph[1]->Fill(goodPhoton.at(1).v.Eta());
      h.subph[2]->Fill(goodPhoton.at(1).v.Phi());

    }
    

    //Plotting dphi, dR and  M_{#gamma#gamma} between the leading and subleading photon
    if((int)goodPhoton.size()>1){
      float pdphi = goodPhoton.at(0).v.DeltaPhi(goodPhoton.at(1).v);
      float pdR = goodPhoton.at(0).v.DeltaR(goodPhoton.at(1).v);
      h.diphoton[0]->Fill(pdphi);
      h.diphoton[1]->Fill(pdR);
      float diph_invmass = (goodPhoton.at(0).v+goodPhoton.at(1).v).M();
      h.diphoton[2]->Fill(diph_invmass);
    } 
    
    
    //###########################################################################

    //For Jets



    //Plotting the leading jet pT in each event.
    if((int)goodJet.size()>0){             
      h.jet[0]->Fill(goodJet.at(0).v.Pt());
    }
    

    //Plotting dphi, dR and  M_{jj} between the leading and subleading jet
    if((int)goodJet.size()>1){
      float jdphi = goodJet.at(0).v.DeltaPhi(goodJet.at(1).v);
      float jdR = goodJet.at(0).v.DeltaR(goodJet.at(1).v);
      h.dijet[0]->Fill(jdphi);
      h.dijet[1]->Fill(jdR);
      float dijet_invmass = (goodJet.at(0).v+goodJet.at(1).v).M();
      h.dijet[2]->Fill(dijet_invmass);
    } 
    

     
    //Finding match of each electron to its closest genelectron

    std::pair<vector<int>, vector<float>> result_el = dR_matching(goodElectron, genElec);
    vector<int> el_matchto_genel = result_el.first;
    vector<float> el_delRmin_genel = result_el.second;

    for(int i=0; i<(int)el_matchto_genel.size(); i++){
      int elmatch=el_matchto_genel.at(i);
      if(elmatch>-1){
	int elmomid = genElec.at(elmatch).momid;
	float elmatchdR = el_delRmin_genel.at(i);
	h.motherID[0]->Fill(elmomid);
	h.dr[0]->Fill(elmatchdR);
      }
    }

 //Finding match of each muon to its closest genmuon

    std::pair<vector<int>, vector<float>> result_mu = dR_matching(goodMu, genMu);
    vector<int> mu_matchto_genmu = result_mu.first;
    vector<float> mu_delRmin_genmu = result_mu.second;

    for(int i=0; i<(int)mu_matchto_genmu.size(); i++){
      int mumatch=mu_matchto_genmu.at(i);
      if(mumatch>-1){
	int mumomid = genMu.at(mumatch).momid;
	float mumatchdR = mu_delRmin_genmu.at(i);
	h.motherID[1]->Fill(mumomid);
	h.dr[1]->Fill(mumatchdR);
      }
    }
    
 
    
    
    //###########################################################################
 
    //###########################################################################

    
    //For gen particles :

    if((int)genMu.size()>0){
      h.genMuProp[0]->Fill(genMu.at(0).v.Pt());
      h.genMuProp[1]->Fill(genMu.at(0).v.Eta());
      h.genMuProp[2]->Fill(genMu.at(0).v.Eta());

      }

    

    if((int)genMu.size()>1){
      h.subgenMuProp[0]->Fill(genMu.at(1).v.Pt());
      h.subgenMuProp[1]->Fill(genMu.at(1).v.Eta());
      h.subgenMuProp[2]->Fill(genMu.at(1).v.Phi());
    }

    
    if((int)genMu.size()>1){
      //cout<<"test"<<endl;
      float digenmu_invmass = (genMu.at(0).v+genMu.at(1).v).M();
      h.gen_invmass[0]->Fill(digenmu_invmass);
    }

    if((int)genElec.size()>0){
      h.genElProp[0]->Fill(genElec.at(0).v.Pt());
      h.genElProp[1]->Fill(genElec.at(0).v.Eta());
      h.genElProp[2]->Fill(genElec.at(0).v.Phi());
    }

    if((int)genElec.size()>1){
      h.subgenElProp[0]->Fill(genElec.at(1).v.Pt());
      h.subgenElProp[1]->Fill(genElec.at(1).v.Eta());
      h.subgenElProp[2]->Fill(genElec.at(1).v.Phi());
    }


    if((int)genElec.size()>1){
      float digenel_invmass = (genElec.at(0).v+genElec.at(1).v).M();
      h.gen_invmass[1]->Fill(digenel_invmass);
    }

     
    
    //########### ANALYSIS ENDS HERE ##############
  }//GoodEvt
  
  
  return kTRUE;
}


//######################################
//        USER DEFINED FUNCTIONS
//######################################


void nano9Ana::Sortpt(vector<Lepton> vec)
{
  
  for(int i=0; i<(int)vec.size()-1; i++){
    for(int j=i+1; j<(int)vec.size(); j++){
      if( vec[i].v.Pt() < vec[j].v.Pt() ) swap(vec.at(i), vec.at(j));
    }
  }
}

int nano9Ana::MotherID(int partindex, int momindex)
{
  int parid = GenPart_pdgId[partindex];
  int momid = GenPart_pdgId[momindex];
  while(parid==momid){
    partindex=momindex;
    momindex=GenPart_genPartIdxMother[momindex];
    parid =GenPart_pdgId[partindex];
    momid = GenPart_pdgId[momindex];
  }
  return momid;
}


float nano9Ana::delta_phi(float phi1, float phi2)
{
  //The correct deltaPhi falls in the interval [0 , pi]
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}

float nano9Ana::transv_mass(float E_lep, float MET, float dphi)
{
  //The inputs are the Energy of the lepton, MET and dPhi between the lepton and MET
  float mT = sqrt(2* E_lep * MET *(1-cos(dphi)));
  return mT;
}


std::pair<vector<int>, vector<float>> nano9Ana::dR_matching(vector<Lepton> vec1, vector<Lepton> vec2)
{
  float delR_min = 999; int match = -1;
  vector<int> foundMatch;
  vector<float> delRmin;
  
  if(vec1.size()>0 && vec2.size()>0){
    for(int i=0; i<(int)vec1.size(); i++){
    for(int j=0; j<(int)vec2.size(); j++){
      float delR = (vec1.at(i).v).DeltaR(vec2.at(j).v);
      if(delR<delR_min){
	delR_min=delR; match=j;
      }
    }
    foundMatch.push_back(match);
    delRmin.push_back(delR_min);
    }
  }
  return std::make_pair(foundMatch, delRmin);
}

void nano9Ana::BookHistograms()
{
  //The histograms are booked here.
  //Binning etc are done here.
  //These histograms are stored in the hst_<process name>.root file in the same order.

  //Example : new TH1F ("hst_name", "hst title", NBins, startVal, EndVal);
  
  h.nmuons = new TH1F("nmuons", "number of muons", 10, 0, 10);
  h.nelectrons = new TH1F("nelectrons", "number of electrons", 10, 0, 10);
  h.nleptons = new TH1F("nlepton", "number of leptons", 10, 0, 10);
  h.nphotons = new TH1F("nphoton", "number of photons", 10, 0, 10);
  h.njets = new TH1F("njet", "number of jets",10, 0, 10);
  h.bjets = new TH1F("bjet", "number of btagged jets", 10, 0, 10);
    
 
  h.mu[0] = new TH1F("mu_pt", "Muon p_{T}", 200, 0, 200);
  h.mu[1] = new TH1F("mu_eta", "Muon eta", 120, -3., 3.);
  h.mu[2] = new TH1F("mu_phi", "Muon phi", 64, -3.2, 3.2);
  
  h.muprop[0] = new TH1F("mu0_pt","Leading muon p_{T}",200,0,200);
  h.muprop[1] = new TH1F("mu0_eta","Leading muon #eta",120,-3.,3.);
  h.muprop[2] = new TH1F("mu0_phi", "Leading muon #phi", 64, -3.2, 3.2);
  h.muprop[3] = new TH1F("mu0_isol","Leading muon PFRelIso",100,0,2);

  h.submuprop[0] = new TH1F("mu1_pt", "Subleading muon p_{T}",200, 0, 100);
  h.submuprop[1] = new TH1F("mu1_eta", "Subleading muon #eta",120, -3., 3.);
  h.submuprop[2] = new TH1F("mu1_phi", "Subleading muon #phi", 64, -3.2, 3.2);

  h.dimuon[0] = new TH1F("dimu_deltaphi", "#Delta#phi_{#mu#mu}", 64, -3.2, 3.2);
  h.dimuon[1] = new TH1F("dimu_deltaR", "#DeltaR_{#mu#mu}", 200, 0, 10);
  h.dimuon[2] = new TH1F("dimu_invmass", "M_{#mu#mu}", 200, 0, 200);

  h.eprop[0] = new TH1F("e0_pt", "Leading electron p_{T}",250,0,250);
  h.eprop[1] = new TH1F("e0_eta", "Leading electron #eta",140,-3.5,3.5);
  h.eprop[2] = new TH1F("e0_phi", "Leading electron #phi",64,-3.2,3.2);
  
  h.subeprop[0] = new TH1F("e1_pt", "Subleading electron p_{T}",200,0,200);
  h.subeprop[1] = new TH1F("e1_eta", "Subleading electron #eta",140,-3.5,3.5);
  h.subeprop[2] = new TH1F("e1_phi", "Subleading electron #phi",64,-3.2,3.2);

  h.dielectron[0] = new TH1F("diel_deltaphi", "#Delta#phi_{ee}",64, -3.2, 3.2);
  h.dielectron[1] = new TH1F("diel_deltaR", "#DeltaR_{ee}",200, 0, 10);
  h.dielectron[2] = new TH1F("diel_invmass", "M_{ee}",200, 0, 200);

  h.ph[0] = new TH1F("gamma0_pt", "Leading photon p_{T}",250,0,250);
  h.ph[1] = new TH1F("gamma0_eta", "Leading photon #eta",140,-3.5,3.5);
  h.ph[2] = new TH1F("gamma0_phi", "Leading photon #phi",64,-3.2,3.2);

  h.subph[0] = new TH1F("gamma1_pt", " #gamma1 p_{T}",250,0,250);
  h.subph[1] = new TH1F("gamma1_eta", "#gamma1 #eta",140,-3.5,3.5);
  h.subph[2] = new TH1F("gamma1_phi", "#gamma1 #phi",64,-3.2,3.2); 
 
  h.diphoton[0] = new TH1F("diph_deltaphi", "d#phi_{#gamma#gamma}", 64, -3.2, 3.2);
  h.diphoton[1] = new TH1F("diph_deltaR", "#DeltaR_{#gamma#gamma}", 200, 0, 10);
  h.diphoton[2] = new TH1F("diph_invmass", "M_{#gamma#gamma}", 200, 0, 200);
  
  h.jet[0] = new TH1F("jet_pT", "jet p_{T}",250,0,250);
  
  
  h.dijet[0] = new TH1F("jet_deltaphi", "#Delta#phi_{jj}",64, -3.2, 3.2);
  h.dijet[1] = new TH1F("jet_deltaR", "#DeltaR_{jj}",200, 0, 10);
  h.dijet[2] = new TH1F("M_{jj}", "M_{jj}",200, 0, 200);
  
  h.btag = new TH1F("btag", "btag", 400, -2, 2);
  
  //h.dr_egamma = new TH1F("dR_e0_gamma", "#DeltaR_{e#gamma}",200, 0, 10);
  //h.dr_ejet = new TH1F("dR_e0_jet", "#DeltaR_{ej}",200, 0, 10);
  
  h.genMuProp[0] = new TH1F("genmu0_pt", "#mu0^{(gen)} p_{T}",250,0,250);
  h.genMuProp[1] = new TH1F("genmu0_eta","#mu0^{(gen)} #eta",140,-3.5,3.5);
  h.genMuProp[2] = new TH1F("genmu0_phi","#mu0^{(gen)} #phi",64,-3.2,3.2);

  h.subgenMuProp[0] = new TH1F("genmu1_pt", "#mu1^{(gen)} p_{T}",250,0,250);
  h.subgenMuProp[1] = new TH1F("genmu1_eta","#mu1^{(gen)} #eta",140,-3.5,3.5);
  h.subgenMuProp[2] = new TH1F("genmu1_phi","#mu1^{(gen)} #phi",64,-3.2,3.2);

  h.gen_invmass[0] = new TH1F("genmu_invmass","M_{#mu#mu}^{gen}",200, 0, 200);

  h.genElProp[0] = new TH1F("genel0_pt", "e0^{(gen)} p_{T}",250,0,250);
  h.genElProp[1] = new TH1F("genel0_eta","e0^{(gen)} #eta",140,-3.5,3.5);
  h.genElProp[2] = new TH1F("genel0_phi","e0^{(gen)} #phi",64,-3.2,3.2);

  h.subgenElProp[0] = new TH1F("genel1_pt", "e1^{(gen)} p_{T}",250,0,250);
  h.subgenElProp[1] = new TH1F("genel1_eta","e1^{(gen)} #eta",140,-3.5,3.5);
  h.subgenElProp[2] = new TH1F("genel1_phi","e1^{(gen)} #phi",64,-3.2,3.2);

  h.gen_invmass[1] = new TH1F("genel_invmass","M_{ee}^{gen}",200, 0, 200);
  
  h.met[0] = new TH1F("met_pt", "Met p_{T}", 200, 0, 200);
  h.met[1] = new TH1F("met_phi", "Met #phi", 64, -3.2, 3.2);

  h.motherID[0] = new TH1F("elecron_mom_id", "electron mom_id", 1200, -600, 600);
  h.motherID[1] = new TH1F("muon_mom_id", "muon mom_id", 1200, -600, 600);

  h.transvmass[0] = new TH1F("dimu_transv_mass", "dimuon m_{T}", 200, 0, 200);
  h.transvmass[1] = new TH1F("diel_transv_mass", "dielectron m_{T}", 200, 0, 200);
  h.transvmass[2] = new TH1F("elmu_transv_mass", "e#mu m_{T}", 200, 0, 200);
  
  h.dr[0] = new TH1F("el_drmin", "dR electron-genelectron", 100, 0, 100);
  h.dr[1] = new TH1F("mu_drmin", "dR electron-genmuon", 100, 0, 100);
  


  
  //############################################################################################################################
  
  
}
