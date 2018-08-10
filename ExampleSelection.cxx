#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <TH2D.h>
#include <TH1D.h>
#include <TMath.h>
#include <THStack.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "ExampleSelection.h"
#include "ExampleTools.h"
#include "TPad.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include <fstream>
#include "nusimdata/SimulationBase/MCParticle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include <stdio.h>
#include <time.h>
#include <string>
#include <cstring>
#include <TLorentzVector.h>
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include <random>
#include "TVector3.h"

namespace ana {
  namespace ExampleAnalysis {

ExampleSelection::ExampleSelection() : SelectionBase(), fNuCount(0), fEventCounter(0) {}

std::vector<TH1D*> ColorFill(std::vector<TH1D*> hists, std::mt19937 rng){
  //random number distribution in range [0,11]
  std::uniform_int_distribution<std::mt19937::result_type> color12(0,11);
  for(int i=0;i<hists.size();i++){
    bool color = true; //no repeated colors!
    do{
       auto num = color12(rng);
       if(num==0) hists[i]->SetFillColor(kOrange);
       else if(num==1) hists[i]->SetFillColor(kRed+2);
       else if(num==2) hists[i]->SetFillColor(kPink+1);
       else if(num==3) hists[i]->SetFillColor(kMagenta-4);
       else if(num==4) hists[i]->SetFillColor(kViolet-5);
       else if(num==5) hists[i]->SetFillColor(kBlue);
       else if(num==6) hists[i]->SetFillColor(kAzure+1);
       else if(num==7) hists[i]->SetFillColor(kCyan+2);
       else if(num==8) hists[i]->SetFillColor(kTeal-2);
       else if(num==9) hists[i]->SetFillColor(kGreen+2);
       else if(num==10) hists[i]->SetFillColor(kSpring+8);
       else hists[i]->SetFillColor(kYellow+1);
       if(i!=0){
         for(int k=i-1; k>=0; k--){
           if(hists[i]->GetFillColor()==hists[k]->GetFillColor()){
             color=false;
             break;
           }
           else color=true;           
         }
       }
    } while(color==false);
  }
   return hists;
}

std::vector<TH1D*> InitializeHists(int nbins, double lowX, double highX, int size, std::string baseTitle, std::mt19937 rng){
   std::vector<TH1D*> hists;
   char buffer[20];
   //srand(time(NULL));
   for(int i=0;i<size;i++){
      char* cstr = new char [baseTitle.length()+1];
      std::strcpy (cstr, baseTitle.c_str());
      std::sprintf(buffer, "%s%d", cstr, i); //in this case, using i to denote 1 (true = events that passed the cut) vs 0 (false)
      hists.push_back(new TH1D(buffer,"",nbins,lowX,highX));
   }
   hists = ColorFill(hists, rng);
   return hists;
}

std::vector<TH1D*> CombineHists(std::vector<TH1D*> h1, std::vector<TH1D*> h2, std::mt19937 rng){
  std::vector<TH1D*> total;
  for(size_t i=0; i<h1.size();i++){
    total.push_back(h1.at(i));
  }
  for(size_t j=0;j<h2.size();j++){
    total.push_back(h2.at(j));
  }
  total = ColorFill(total, rng);
  return total;
}

void WriteHists(std::vector<TH1D*> vec, THStack* stack){
  for(size_t i=0; i< vec.size(); i++){
    stack->Add(vec[i]);
  }
  stack->Write();
}

void ExampleSelection::Initialize(Json::Value* config) {
  familyConfusion.open("showerPDGcodes.txt");  

  rng.seed(std::random_device()());

  // Make a histogram
  fNuVertexXZHist = new TH2D("nu_vtx_XZ", "",100, -1000, 1000, 100, -1000, 1000);

  nuFlux = new TH1D("flux", "#nu flux;E_#nu;#nu count", 30,0, 5);
  E_hadron = new TH2D("E_hadron", "Charged hadronic activity;E_#nu (GeV);E_hadron (GeV)", 30, 0, 5, 15, 0, 1); //all charged hadron activity that came from a neutrino interaction

  fCut1 = InitializeHists(30, 0, 5, 2, "#nu_e CC", rng);  
  fnueCC = new THStack("fnueCC", "Intrinsic/signal #nu_e CC;E_#nu (GeV);count");
  fCut1_reco = InitializeHists(30,0,5,2,"e",rng);
  fnueCC_reco = new THStack("fnueCC_reco", "Cut 1 with reco E_#nu;Reconstructed E_#nu (GeV);count");
  fig11 = InitializeHists(30,0,5, 5,"thing", rng);
  fig11Stack = new THStack("fig11Stack",";Reconstructed Energy (GeV);Events/GeV");

  E_shower_nu = new TH2D("E_shower_nu", "Shower energy vs #nu_e energy;E_#nu (GeV);E_shower (GeV)",30,0,5,36,0,6);

  fCut2 = InitializeHists(30,0,5,2, "events", rng);
  fGammaNC = new THStack("fGammaNC", "NC #gamma production;E_#nu (GeV);count");
  E_gamma = InitializeHists(30,0,500,5,"#gamma", rng); //0 is all photons; 1 is all primary photons; 2 is all primary photons with >200 MeV; 3 is all secondary photons; 4 is all secondary photons with >100 MeV.
  gammaStack = new THStack("gammaStack", "Photon energies;E_#gamma (MeV);count");
  hPass = InitializeHists(30,0,5,3, "chHadron", rng);
  hadStack = new THStack("hadStack", "Charged hadron activity;E_#nu (GeV);count");
  pPass = InitializeHists(30,0,500,3, "photon", rng);
  phoStack = new THStack("phoStack", "Photon/interaction vertex cuts;E_#nu (GeV);count");
  dEdx = new TH1D("dEdx","Shower dE/dx;dE/dx (MeV/cm);particle count", 30, 0, 5);
  dEdx_2 = InitializeHists(30,0,5,2,"particle",rng);
  dEdx_2_stack = new THStack("showerStack","Shower dE/dx;dE/dx(MeV/cm);count");

  fCut3 = InitializeHists(30,0,5,2,"#nu_mu CC", rng);
  fnumuCC = new THStack("fnumuCC","CC #nu_#mu;E_#nu (GeV);count");

  fnumuCut = InitializeHists(30,0,5,2,"muons", rng);
  fnumuCandidates = new THStack("fnumuCandidates","#nu_#mu CC candidates");

  fmumu = InitializeHists(30,0,5,3,"mu_passing", rng);
  fpie = InitializeHists(30,0,5,3,"pi0_passing", rng);
  fDiscriminate_numu_mupi = new THStack("fDiscriminate_numu_mupi","#nu_#mu CC candidates;E_#nu (GeV);count");   
  reco_compare = InitializeHists(30,0,5,2,"mupi",rng);
  reco_compare_stack = new THStack("reco_compare_stack","#nu_#mu CC Candidates;E_#nu (GeV);count");

  fnumu_reco1 = new TH1D("fnumu_reco1","Reconstructed #nu_#mu cut;Reconstructed Energy (GeV);count",30,0,5);
  fnumu_reco2 = InitializeHists(30,0,5,2,"mupi",rng);
  fnumu_reco_stack = new THStack("fnumu_reco_stack","Reconstructed #nu_#mu CC cand;Reconstructed Energy (GeV);count");

  // Load configuration parameters
  fMyParam = 0;
  fTruthTag = { "generator" };
  fGenieTag = {"genieeventweight"};
  fFluxTag = {"fluxeventweight"};
  fParticleTag = {"largeant"};
  fTrackTag = {"mcreco"};

  if (config) {
    fMyParam = (*config)["ExampleAnalysis"].get("parameter", 0).asInt();
    fTruthTag = { (*config)["ExampleAnalysis"].get("MCTruthTag", "generator").asString() };
    fGenieTag = { (*config)["ExampleAnalysis"].get("MCEventWeightTag", "genieeventweight").asString() };
    fFluxTag = { (*config)["ExampleAnalysis"].get("MCEventWeightTag", "fluxeventweight").asString() };
    fParticleTag = { (*config)["ExampleAnalysis"].get("MCParticleTag", "largeant").asString() };
    fTrackTag = { (*config)["ExampleAnalysis"].get("MCParticleTag", "mcreco").asString() };
  }

  // Add custom branches
  AddBranch("nucount", &fNuCount);
  AddBranch("myvar", &fMyVar);
  AddBranch("interactionType", &fInteraction);
  AddBranch("parent", &fPParticle);
  AddBranch("nuenergy", &fEnergy);
  AddBranch("evtwts", &fEventWt);
  AddBranch("totFluxWts", &fTotFluxWt);

  // Use some library code
  hello();
}

/**
 * Find total flux weights
 */
void FindFluxWeights(evwgh::MCEventWeight fluxwt, std::vector<std::vector<double>> fTotFluxWt, std::vector<double> fWts){
  auto size = fluxwt.fWeight.at("horncurrent_FluxUnisim").size();
    for(size_t universe=0; universe<size; universe++){
      double weight = 1.0;
      for(auto entry: fluxwt.fWeight){
        //entry is a pair of <string, vector<double>> 
        //want to always skip anything to do with the first entry/parameter
        //after universe 0, since only the first parameter has only 1 
        //universe, whereas the rest all have 1000. Otherwise will 
        //seg fault since I'm going past its max index/size
        if(universe != 0 && entry.second.size() < size) continue;
        weight *= entry.second[universe];
      }
      fWts.push_back(weight);
    }
    fTotFluxWt.push_back(fWts);
    fWts.clear();
}

bool AnybodyHome_SBND(TLorentzVector pos){
   if(((-199.15 < pos.X() && pos.X() < -2.65) || (2.65 < pos.X() && pos.X() < 199.15)) && (-200 < pos.Y() && pos.Y() < 200) && (0 < pos.Z() && pos.Z() < 500)) return true;
   else return false;
}

bool AnybodyHome_muBooNE(TLorentzVector pos){
   if((-1.55 < pos.X() && pos.X() < 254.8) && (-115.53 < pos.Y() && pos.Y() < 117.47) && (0.1 < pos.Z() && pos.Z() < 1036.9)) return true;
   else return false;
}
 
bool AnybodyHome_ICARUS(TLorentzVector pos){
   if(((-364.49 < pos.X() && pos.X() < -216.29) || (-216.14 < pos.X() && pos.X() < -67.94)) && (-173.41 < pos.Y() && pos.Y() < 143.41) && (-909.950652 < pos.Z() && pos.Z() < 879.950652)) return true;
   else if(((67.94 < pos.X() && pos.X() < 216.14) || (216.29 < pos.X() && pos.X() < 364.49)) && (-173.41 < pos.Y() && pos.Y() < 143.41) && (-909.950652 < pos.Z() && pos.Z() < 879.950652)) return true;
   else return false;
}

double FindDistance(TLorentzVector a, TLorentzVector b){
  auto distance = TMath::Sqrt(TMath::Power(a.X()-b.X(),2) + TMath::Power(a.Y()-b.Y(),2) + TMath::Power(a.Z()-b.Z(),2));
  return distance;
}

double ActiveTrackLength(sim::MCTrack mctrack){
  double tracklength = 0.0;
  for(size_t i=0; i<mctrack.size();i++){
    sim::MCStep& step = mctrack[i];
    if(i>0){
      if(AnybodyHome_SBND(step.Position())==true && AnybodyHome_SBND(mctrack[i-1].Position())==true){
        tracklength += FindDistance(step.Position(),mctrack[i-1].Position());
      }
    }
  }
  return tracklength;
}

std::vector<sim::MCTrack> FindRelevantTracks(simb::MCTruth mctruth, std::vector<sim::MCTrack> mctracks){
  std::vector<sim::MCTrack> relTracks;
  auto nuVertex = mctruth.GetNeutrino().Nu().EndPosition();
  //iterate through only each neutrino's "associated" tracks + showers (within 5 cm of neutrino interaction vertex)
  //(we're cheating - it's okay)
  for(size_t r=0;r<mctracks.size();r++){
    auto const& mctrack = mctracks.at(r);
    if(FindDistance(mctrack.Start().Position(),nuVertex)<=5){
      relTracks.push_back(mctrack);
    }
  }
  return relTracks;
}

std::vector<sim::MCShower> FindRelevantShowers(simb::MCTruth mctruth, std::vector<sim::MCShower> mcshowers){
  std::vector<sim::MCShower> relShowers;
  auto nuVertex = mctruth.GetNeutrino().Nu().EndPosition();
  //iterate through only each neutrino's "associated" tracks + showers (within 5 cm of neutrino interaction vertex)
  for(size_t s=0;s<mcshowers.size();s++){
    auto const& mcshower = mcshowers.at(s);
    if(FindDistance(mcshower.Start().Position(),nuVertex)<=5){
      relShowers.push_back(mcshower);
    }
  }
  return relShowers;
}

/**
 * Return vector containing unique PDG codes of primary particles in MCTrack or MCShower, for ALL events
 */
std::vector<int> UniquePDGs(std::vector<int> codes, std::vector<sim::MCTrack> fRelTracks,std::vector<sim::MCShower> fRelShowers, int option){
  if(option==1){
    for(size_t s=0; s<fRelShowers.size();s++){
      auto const& mcshower = fRelShowers.at(s);
      if(mcshower.Process()=="primary"){
        if(codes.empty()){
          codes.push_back(mcshower.PdgCode());
        }
        else{
          if(find(codes.begin(),codes.end(),mcshower.PdgCode())==codes.end()){
            codes.push_back(mcshower.PdgCode());
          }
        }
      }
    }
  }
  //if unspecified, get relevant track codes
  else{
    for(size_t s=0; s<fRelTracks.size();s++){
      auto const& mctrack = fRelTracks.at(s);
      if(mctrack.Process()=="primary"){
        if(codes.empty()){
          codes.push_back(mctrack.PdgCode());
        }
        else{
          if(find(codes.begin(),codes.end(),mctrack.PdgCode())==codes.end()){
            codes.push_back(mctrack.PdgCode());
          }
        }
      }
    }
  }
  return codes;
}

/**
* reconstructed neutrino energy based on cross-section for CCQE interactions w/ heavy nuclei
*/
double FindRecoEnergy_nue(sim::MCShower mcshower){
  double m_n = 939.565; //neutron mass, in MeV
  double m_p = 938.272; //proton mass, in MeV
  double E_b = 340; //binding energy to bind with an Argon nucleus, in MeV
  double m_l = 0.510; //electron mass, in MeV
  auto E_l = mcshower.Start().E(); //lepton energy, in MeV
  auto p_l = 0.001*mcshower.Start().Momentum().Vect().Mag(); //lep momentum, MeV
  auto theta_l = mcshower.Start().Position().Vect().Theta(); //azimuthal angle b/w prod lepton + neutrino (z-axis)
  auto reco_energy = 0.5*(TMath::Power(m_p,2)-TMath::Power(m_n-E_b,2)-TMath::Power(m_l,2)+2*(m_n-E_b)*E_l)/(m_n-E_l+p_l*TMath::Cos(theta_l));
  reco_energy /= 1000; //convert from MeV to GeV
  return reco_energy;
}

double FindRecoEnergy_nue(sim::MCTrack mctrack){
  double m_n = 939.565; //neutron mass, in MeV
  double m_p = 938.272; //proton mass, in MeV
  double E_b = 340; //binding energy to bind with an Argon nucleus, in MeV
  double m_l = 0.510; //electron mass, in MeV
  auto E_l = mctrack.Start().E(); //lepton energy, in MeV
  auto p_l = 0.001*mctrack.Start().Momentum().Vect().Mag(); //lep momentum, MeV
  auto theta_l = mctrack.Start().Position().Vect().Theta(); //azimuthal angle b/w prod lepton + neutrino (z-axis)
  auto reco_energy = 0.5*(TMath::Power(m_p,2)-TMath::Power(m_n-E_b,2)-TMath::Power(m_l,2)+2*(m_n-E_b)*E_l)/(m_n-E_l+p_l*TMath::Cos(theta_l));
  reco_energy /= 1000; //convert from MeV to GeV
  return reco_energy;
}

double FindRecoEnergy_numu(sim::MCTrack mctrack){
  double m_n = 939.565; //neutron mass, in MeV
  double m_p = 938.272; //proton mass, in MeV
  double E_b = 340; //binding energy to bind with an Argon nucleus, in MeV
  auto m_l = 105.658; //muon mass, in MeV
  auto E_l = mctrack.Start().E(); //lepton energy, in MeV
  auto p_l = 0.001*mctrack.Start().Momentum().Vect().Mag(); //lep momentum, MeV
  auto theta_l = mctrack.Start().Position().Vect().Theta(); //azimuthal angle b/w prod lepton + neutrino (z-axis)
  auto reco_energy = 0.5*(TMath::Power(m_p,2)-TMath::Power(m_n-E_b,2)-TMath::Power(m_l,2)+2*(m_n-E_b)*E_l)/(m_n-E_l+p_l*TMath::Cos(theta_l));
  reco_energy /= 1000;
  return reco_energy;
}
  

/**
 * nu_e CC CUT 1 IMPLEMENTATION: nu_e CC signal, 80% efficiency
 */
void Cut1nue(simb::MCTruth mctruth, std::vector<sim::MCTrack> fRelTracks,std::vector<sim::MCShower> fRelShowers, std::vector<TH1D*> fCut1, std::mt19937 rng){
  std::uniform_int_distribution<std::mt19937::result_type> dist100(0,99); //distribution in range [0, 99]
  auto CCNC = mctruth.GetNeutrino().CCNC();
  auto energy = mctruth.GetNeutrino().Nu().E();
  if(CCNC==0){
    for(size_t t=0; t<fRelShowers.size();t++){
      auto mcshower = fRelShowers.at(t);
      E_shower_nu->Fill(energy, mcshower.Start().E()/1000, 1);
      //if shower created by electron
      if(mcshower.PdgCode()==11 && mcshower.Process()=="primary"){
        if(mcshower.Start().E()/1000>0.2){
          if(dist100(rng)<80) fCut1[1]->Fill(energy, 1);
          else fCut1[0]->Fill(energy, 1);
        }
        else fCut1[0]->Fill(energy, 1);
      }
    }
  }
}

int LinkShower_HadronParent(simb::MCFlux mcflux, simb::MCTruth mctruth, sim::MCShower mcshower){
  int parent = 0;
  if(mctruth.GetNeutrino().Lepton().PdgCode()==mcshower.PdgCode()){
    if(mctruth.GetNeutrino().Lepton().Position().Vect()==mcshower.Start().Position().Vect()){
      parent = mcflux.fptype;
    }
  }
  return parent;
}

void Cut1nue_reco(simb::MCFlux mcflux, simb::MCTruth mctruth, std::vector<sim::MCShower> mcshowers, std::vector<TH1D*> fCut1, std::vector<TH1D*> fig11, std::mt19937 rng){
  std::uniform_int_distribution<std::mt19937::result_type> dist100(0,99); //distribution in range [0, 99]
  for(size_t i=0;i<mcshowers.size();i++){
    auto mcshower = mcshowers.at(i);
    if(mcshower.PdgCode()==11){
      auto energy = FindRecoEnergy_nue(mcshower);
      if(mcshower.Start().E() > 200){
        if(dist100(rng)<80){
          //muon parent
          if(LinkShower_HadronParent(mcflux,mctruth,mcshower)==13) fig11[0]->Fill(energy, 1);
          //K0 parent
          if(LinkShower_HadronParent(mcflux,mctruth,mcshower)==130 || LinkShower_HadronParent(mcflux,mctruth,mcshower)==310 || LinkShower_HadronParent(mcflux,mctruth,mcshower)==311) fig11[1]->Fill(energy, 1);
          //K+ parent
          if(LinkShower_HadronParent(mcflux,mctruth,mcshower)==321){
            fig11[2]->Fill(energy, 1);
            familyConfusion << "hadron parent: " << LinkShower_HadronParent(mcflux,mctruth,mcshower) << "\treco energy: " << energy <<std::endl;
          }

          fCut1[1]->Fill(energy, 1);
        }
        else fCut1[0]->Fill(energy, 1);
      }
      else fCut1[0]->Fill(energy, 1);
    }
  }
}

/**
 * nu_e CC CUT 2 IMPLEMENTATION: NC photon production
 */
void Cut2nue(simb::MCTruth mctruth, std::vector<sim::MCTrack> fRelTracks, std::vector<sim::MCShower> fRelShowers, std::vector<TH1D*> hists, std::vector<TH1D*> fig11, std::mt19937 rng){
  std::uniform_int_distribution<std::mt19937::result_type> dist100(0,99); //distribution in range [0, 99]
  auto nuE = mctruth.GetNeutrino().Nu().E();
  auto nuVertex = mctruth.GetNeutrino().Nu().EndPosition(); 
  int count = 0; //count the number of photon showers produced by daughter pi0 decay
  double hadronE = 0.0; //total charged hadronic activity produced at nu vertex
  bool visible = false; //neutrino has "visible" interaction vertex if hadronE > 50 MeV at vertex
  bool reject = false; //does the neutrino interaction event pass the cuts?
  //keep events with photon showers above 200 MeV
  int nKids=0; //number of kids
  int momID;
  for(size_t n=0;n<fRelShowers.size();n++){
    auto const& mcshower = fRelShowers.at(n);
    dEdx->Fill(mcshower.dEdx(),1);
    
    if(mcshower.PdgCode()==22){
      E_gamma[0]->Fill(mcshower.Start().E(),1); //all photon showers
      pPass[0]->Fill(mcshower.Start().E(),1);
      dEdx_2[0]->Fill(mcshower.dEdx(),1);
    }
    if(mcshower.PdgCode()==11) dEdx_2[1]->Fill(mcshower.dEdx(),1);

    //leading photon shower produced by pi0 must have E_gamma > 200 MeV; 2nd photon shower must have E_gamma > 100 MeV
    if(mcshower.PdgCode()==22 && mcshower.MotherPdgCode()==111 && mcshower.MotherProcess()=="primary"){
      if(nKids==0){
        E_gamma[1]->Fill(mcshower.Start().E(),1); //# leading photon showers
        momID = mcshower.MotherTrackID(); //parent track ID of the first photon child
        if(mcshower.Start().E() > 200){
          E_gamma[2]->Fill(mcshower.Start().E(), 1); //# primary g-showers that pass the energy cut
          count++;
        }
        nKids++;
      }      
      else if(nKids==1){
        if(mcshower.MotherTrackID()==momID){
          E_gamma[3]->Fill(mcshower.Start().E(),1); //# 2nd photon showers
          if(count==1 && mcshower.Start().E() > 100){
            E_gamma[4]->Fill(mcshower.Start().E(), 1); //#2nd g-kids that pass the energy cut
            count++;
            //if the second photon converts (showers) within the active TPC, reject the event!
            if(AnybodyHome_SBND(mcshower.Start().Position())==true){
              reject=true;
              break; 
            } 
          }
        }
      }
    }
  }

  for(size_t r=0;r<fRelTracks.size();r++){
    auto const& mctrack = fRelTracks.at(r);
    //nu_e CC SELECTION CUT 2B
    //Look for charged hadron primary particles, i.e. NOT leptons, protons, pi0, or K0. Find energy total (>50 MeV to be visible)
    if(mctrack.PdgCode()!=11 && mctrack.PdgCode()!=13 && mctrack.PdgCode()!=15 && mctrack.PdgCode()!=111 && mctrack.PdgCode()!=130 && mctrack.PdgCode()!= 310 && mctrack.PdgCode()!=311 && mctrack.Process()=="primary"){
      double anarG = mctrack.Start().E(); //in MeV
      hadronE += anarG;
    }
  }
  E_hadron->Fill(nuE, hadronE/1000, 1);  
  hPass[0]->Fill(nuE, 1); //total count of interactions
  if(hadronE >0) hPass[1]->Fill(nuE, 1); //count of number of interactions that produced any charged hadron activity
  if(hadronE > 50){
    visible = true;
    hPass[2]->Fill(nuE, 1); //number of interactions that produced charged hadron activity abv the 50 MeV threshold
  }
  int nPhotonSh = 0; //total num of photon showers in the nu interaction
  int nfailed = 0; //num of photon showers that converted >3cm from vertex
  if(visible==true){
    for(size_t sh=0; sh<fRelShowers.size();sh++){
      auto const& mcshower = fRelShowers.at(sh);
      auto dist_from_vertex = FindDistance(mcshower.Start().Position(), nuVertex);
      if(mcshower.PdgCode()==22){
        nPhotonSh++;
        pPass[1]->Fill(mcshower.Start().E(), 1); //all photon showers produced by an interaction with a visible vertex
        if(dist_from_vertex > 3){
          nfailed++; 
          pPass[2]->Fill(mcshower.Start().E(), 1); //total photon showers that convert >3 cm from vertex
        } 
      } 
    }
    if(nfailed == nPhotonSh){
      reject=true; 
    }
  } 
  //96% photon rejection rate from dE/dx cut - have to loop back thru all relevant showers
  if(reject==false){
    for(size_t sh=0; sh<fRelShowers.size();sh++){
      auto const& mcshower = fRelShowers.at(sh);
      if(mcshower.PdgCode()==22){
        if(mcshower.dEdx()<2.3 && mcshower.dEdx()>1.7){
          hists[1]->Fill(nuE,1);
          fig11[3]->Fill(FindRecoEnergy_nue(mcshower),1);
        }
        else hists[0]->Fill(nuE,1);
        //if(dist100(rng)>96) hists[1]->Fill(nuE, 1);
        //else hists[0]->Fill(nuE, 1);
      }
    }  
  }
  //if rejecting event, add all photons to "failed" hist
  else{
    for(size_t sh=0; sh<fRelShowers.size();sh++){
      auto const& mcshower = fRelShowers.at(sh);
      if(mcshower.PdgCode()==22) hists[0]->Fill(nuE, 1);  
    }
  } 
} 

/**
 * nu_e CC CUT 3 IMPLEMENTATION
 */
void Cut3nue(simb::MCTruth mctruth, std::vector<sim::MCTrack> fRelTracks, std::vector<sim::MCShower> fRelShowers, std::vector<TH1D*> hists, std::vector<TH1D*> fig11, std::mt19937 rng){
  if(mctruth.GetNeutrino().Nu().PdgCode()==14 && mctruth.GetNeutrino().CCNC()==0){
    for(size_t t=0; t<fRelTracks.size();t++){
      auto const& mctrack = fRelTracks.at(t);
      if(mctrack.PdgCode()==13){
        //if visible track length > 1 m or track exits (cannot tell what total track length is), then cut
        if(ActiveTrackLength(mctrack)>=100 || AnybodyHome_SBND(mctrack.End().Position())==false){
          hists[0]->Fill(mctruth.GetNeutrino().Nu().E(),1);  
        }
        else{
          int ct=0; //counts how many primary showers are associated w/ vertex
          for(size_t s=0; s<fRelShowers.size();s++){
            auto const& mcshower = fRelShowers.at(s);
            ct++;
          }
          //if only one shower associated with the CC event vertex 
          if(ct==1){
            //run Cut 2 photon selection criteria, events that are not rejected are retained as BG for nu_e CC sample.
            //I guess only cuts 2+3 (dEdx and produced charged hadron activity) are relevant?
            bool reject = false; 
            bool visible = false;
            double hadronE = 0.0;
            for(size_t r=0;r<fRelTracks.size();r++){
              auto const& track = fRelTracks.at(r);
              if(track.PdgCode()!=11 && track.PdgCode()!=13 && track.PdgCode()!=15 && track.PdgCode()!=111 && track.PdgCode()!=130 && track.PdgCode()!= 310 && track.PdgCode()!=311 && track.Process()=="primary"){
                double anarG = track.Start().E(); //in MeV
                hadronE += anarG;
              }
            }
            if(hadronE > 50){
              visible = true;
            }
            int nPhotonSh = 0; //total num of photon showers in the nu interaction
            int nfailed = 0; //num of photon showers that converted >3cm from vertex
            if(visible==true){
              for(size_t sh=0; sh<fRelShowers.size();sh++){
              auto const& mcshower = fRelShowers.at(sh);
              auto dist_from_vertex = FindDistance(mcshower.Start().Position(), mctruth.GetNeutrino().Nu().EndPosition());
                if(mcshower.PdgCode()==22){
                  nPhotonSh++;
                  if(dist_from_vertex > 3){
                    nfailed++; 
                  }  
                } 
              }
              if(nfailed == nPhotonSh){
              reject=true; 
              }
            }        
            if(reject==false){
              for(size_t sh=0; sh<fRelShowers.size();sh++){
                auto const& mcshower = fRelShowers.at(sh);
                if(mcshower.PdgCode()==22){
                  if(mcshower.dEdx()<2.5 && mcshower.dEdx()>1.5){
                    hists[1]->Fill(mctruth.GetNeutrino().Nu().E(),1);
                    fig11[4]->Fill(FindRecoEnergy_nue(mcshower),1);
                  }
                  else hists[0]->Fill(mctruth.GetNeutrino().Nu().E(),1);
                }
              }  
            }
            //if rejecting event, add all photons to "failed" hist
            else{
              for(size_t sh=0; sh<fRelShowers.size();sh++){
              auto const& mcshower = fRelShowers.at(sh);
              if(mcshower.PdgCode()==22) hists[0]->Fill(mctruth.GetNeutrino().Nu().E(), 1);  
              }
            } 
    
          }
        }
      }
    } 
  }
}

void Cut3nue_reco(simb::MCTruth mctruth, std::vector<sim::MCTrack> fRelTracks, std::vector<sim::MCShower> fRelShowers, std::vector<TH1D*> hists, std::mt19937 rng){  
  for(size_t t=0; t<fRelTracks.size();t++){
    auto const& mctrack = fRelTracks.at(t);
    if(mctrack.PdgCode()==13){
      //if visible track length > 1 m or track exits (cannot tell what total track length is), then cut
      if(ActiveTrackLength(mctrack)>=100 || AnybodyHome_SBND(mctrack.End().Position())==false){
        hists[0]->Fill(FindRecoEnergy_nue(mctrack),1);  
      }
      else{
        int ct=0; //counts how many primary showers are associated w/ vertex
        for(size_t s=0; s<fRelShowers.size();s++){
          auto const& mcshower = fRelShowers.at(s);
          //muon start position as sub for nu vertex pos
          if(FindDistance(mcshower.Start().Position(),mctrack.Start().Position()) < 5) ct++; 
        }
        //if only one shower associated with the CC event vertex 
        //if(ct==1)
      }
    }
  }  
}

/**
 * Nu_mu candidate event cuts
 */
void numuCCCut(simb::MCTruth mctruth, std::vector<sim::MCTrack> fRelTracks, std::vector<TH1D*> hists, std::vector<TH1D*> reco_compare, std::vector<TH1D*> fmumu, std::vector<TH1D*> fpie){
  auto energy = mctruth.GetNeutrino().Nu().E();
  auto nuType = mctruth.GetNeutrino().Nu().PdgCode();
  auto CCNC = mctruth.GetNeutrino().CCNC(); //0 for CC interaction
  for(size_t t=0; t<fRelTracks.size();t++){
    auto const& mctrack = fRelTracks.at(t);
    //auto tracklength = ActiveTrackLength(mctrack);
    auto tracklength = FindDistance(mctrack.Start().Position(), mctrack.End().Position());

    if(mctrack.PdgCode()==13 || mctrack.PdgCode()==-13){
      fmumu[0]->Fill(energy, 1); 
      nummu++;
    }
    else if(mctrack.PdgCode()==211 || mctrack.PdgCode()==-211){
      fpie[0]->Fill(energy, 1); 
      numpi++;
    }
    
    if(AnybodyHome_SBND(mctrack.Start().Position())==true && AnybodyHome_SBND(mctrack.End().Position())==true){
      if(tracklength > 50){
        if(mctrack.PdgCode()==13){ 
          fmumu[1]->Fill(energy, 1); 
          reco_compare[0]->Fill(energy, 1); //passed muons
          nummupass++;
        } 
        else if(mctrack.PdgCode()==-13){
          fmumu[2]->Fill(energy, 1);
        }
        else if(mctrack.PdgCode()==211){ 
          fpie[1]->Fill(energy, 1); 
          reco_compare[1]->Fill(energy,1);//passed pions
          numpipass++;
        }
        else if(mctrack.PdgCode()==-211){
          fpie[2]->Fill(energy, 1); 
        }
        if(mctrack.PdgCode()==13 && mctrack.Process()=="primary") hists[1]->Fill(energy, 1);
      }
      else{
        if(mctrack.PdgCode()==13 && mctrack.Process()=="primary") hists[0]->Fill(energy, 1);  
      }
    }   
    
    //if exiting muon track, length of track in detector active volume must be >= 1 m to pass cut.    
    if(AnybodyHome_SBND(mctrack.Start().Position())==true && AnybodyHome_SBND(mctrack.End().Position())==false){
      if(tracklength >= 100){
        if(mctrack.PdgCode()==13){ 
          reco_compare[0]->Fill(energy,1); 
          fmumu[1]->Fill(energy, 1); 
          nummupass++;
        } 
        else if(mctrack.PdgCode()==-13){
          fmumu[2]->Fill(energy, 1); 
        }
        else if(mctrack.PdgCode()==211){ 
          reco_compare[1]->Fill(energy,1);
          fpie[1]->Fill(energy, 1); 
          numpipass++;
        }
        else if(mctrack.PdgCode()==-211){
          fpie[2]->Fill(energy, 1); 
        }
        if(mctrack.PdgCode()==13 && mctrack.Process()=="primary") hists[1]->Fill(energy, 1);
      }
      else{
        if(mctrack.PdgCode()==13 && mctrack.Process()=="primary") hists[0]->Fill(energy, 1);  
      }
    }      
  } 
}

void numu_reco(std::vector<sim::MCTrack> mctracks, TH1D* fnumu_reco1, std::vector<TH1D*> fnumu_reco2){
  for(size_t t=0; t<mctracks.size(); t++){
    auto mctrack = mctracks.at(t);
    auto tracklength = FindDistance(mctrack.Start().Position(), mctrack.End().Position());
    auto energy = FindRecoEnergy_numu(mctrack);
    if(AnybodyHome_SBND(mctrack.Start().Position())==true && AnybodyHome_SBND(mctrack.End().Position())==true){
      if(tracklength > 50){
        fnumu_reco1->Fill(energy, 1);
        if(mctrack.PdgCode()==13){ 
          fnumu_reco2[0]->Fill(energy, 1); 
        } 
        else if(mctrack.PdgCode()==211){ 
          fnumu_reco2[1]->Fill(energy, 1); 
        }
      }
    }   

    //if exiting muon track, length of track in detector active volume must be >= 1 m to pass cut.    
    if(AnybodyHome_SBND(mctrack.Start().Position())==true && AnybodyHome_SBND(mctrack.End().Position())==false){
      if(tracklength >= 100){
        fnumu_reco1->Fill(energy, 1);
        if(mctrack.PdgCode()==13){ 
          fnumu_reco2[0]->Fill(energy, 1); 
        } 
        else if(mctrack.PdgCode()==211){ 
          fnumu_reco2[1]->Fill(energy, 1); 
        }
      }
    }      
  }
}

void ExampleSelection::Finalize() {
  // Output our histograms to the ROOT file
  familyConfusion.close();
  fOutputFile->cd();
  fNuVertexXZHist->Write();  
  nuFlux->Write();
  E_shower_nu->Write();
  WriteHists(E_gamma, gammaStack);
  dEdx->Write();
  WriteHists(dEdx_2, dEdx_2_stack);
  E_hadron->Write();
  WriteHists(fCut1, fnueCC);
  WriteHists(fCut1_reco, fnueCC_reco);
  fig11[0]->SetFillColor(kSpring-7);
  fig11[1]->SetFillColor(kGreen+2);
  fig11[2]->SetFillColor(kGreen-7);
  fig11[3]->SetFillColor(kOrange+6);
  fig11[4]->SetFillColor(kBlue-9);
  WriteHists(fig11, fig11Stack);
  WriteHists(fCut2, fGammaNC);
  WriteHists(hPass, hadStack);
  WriteHists(pPass, phoStack);
  WriteHists(fCut3, fnumuCC);
  WriteHists(fnumuCut, fnumuCandidates);
  WriteHists(reco_compare, reco_compare_stack);
  WriteHists(CombineHists(fmumu, fpie, rng), fDiscriminate_numu_mupi);    
  fnumu_reco1->Write();
  WriteHists(fnumu_reco2, fnumu_reco_stack);  

  std::cout << "MCTrack primary particle codes: " << std::endl;
  for(size_t i=0;i<fTrackCodes.size();i++){
    std::cout << fTrackCodes[i] << std::endl;
  }
  std::cout <<"MCShower primary codes: " << std::endl;
  for(size_t j=0;j<fShowerCodes.size();j++){
    std::cout << fShowerCodes[j] << std::endl;
  }
}

bool ExampleSelection::ProcessEvent(gallery::Event& ev) {
  srand(time(NULL)); 

  if (fEventCounter % 1000 == 0) {
    std::cout << "ExampleSelection: Processing event " << fEventCounter << std::endl;
  }

  // Grab a data product from the event - gives a vector containing the neutrinos from each event
  // MCTruth refers to the truth information about each neutrino interaction 
  auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth>>(fTruthTag);
  auto const& mcfluxes = *ev.getValidHandle<std::vector<simb::MCFlux>>(fTruthTag);
  assert(mctruths.size() == mcfluxes.size());
  auto const& mcparticles = *ev.getValidHandle<std::vector<simb::MCParticle>>(fParticleTag);

  auto const& evtwts = *ev.getValidHandle<std::vector<evwgh::MCEventWeight>>(fGenieTag); //event weights
  auto const& fluxwts = *ev.getValidHandle<std::vector<evwgh::MCEventWeight>>(fFluxTag);

  auto const& mctracks = *ev.getValidHandle<std::vector<sim::MCTrack>>(fTrackTag);
  auto const& mcshowers = *ev.getValidHandle<std::vector<sim::MCShower>>(fTrackTag);

  //Define handles (like pointers but memory address can be changed? I don't know why this works)
  auto const& mctruths_handle = ev.getValidHandle<std::vector<simb::MCTruth>>(fTruthTag);
  auto const& mcparticles_handle = ev.getValidHandle<std::vector<simb::MCParticle>>(fParticleTag);

  // Fill in the custom branches
  fNuCount = mctruths.size();  // Number of neutrinos in this event
  fMyVar = fMyParam;

  fInteraction.clear();
  fPParticle.clear();
  fEnergy.clear();
  fEventWt.clear();
  fTotFluxWt.clear();
 
  // Iterate through the neutrinos
  for (size_t i=0; i<mctruths.size(); i++) {
    auto const& mctruth = mctruths.at(i);
    auto const& mcflux = mcfluxes.at(i);
    auto const& evtwt = evtwts.at(i);
    auto const& fluxwt = fluxwts.at(i); //flux weights per neutrino

    auto energy = mctruth.GetNeutrino().Nu().E();
    fEnergy.push_back(energy);
    nuFlux->Fill(energy, 1);  

    auto mode = mctruth.GetNeutrino().Mode();
    fInteraction.push_back(mode);
    //Get nu PDG code: mctruth.GetNeutrino().Nu().PdgCode()
    auto nuPDG = mcflux.fntype;

    //Find total flux weights
    FindFluxWeights(fluxwt, fTotFluxWt, fWts);

    // Fill neutrino vertex position histogram
    fNuVertexXZHist->Fill(mctruth.GetNeutrino().Nu().Vx(), mctruth.GetNeutrino().Nu().Vz());

    fRelTracks = FindRelevantTracks(mctruth, mctracks);
    fRelShowers = FindRelevantShowers(mctruth, mcshowers);
    fTrackCodes = UniquePDGs(fTrackCodes, fRelTracks, fRelShowers, 0);
    fShowerCodes = UniquePDGs(fShowerCodes, fRelTracks, fRelShowers, 1); 

    Cut1nue(mctruth, fRelTracks, fRelShowers, fCut1, rng);
    Cut1nue_reco(mcflux, mctruth, mcshowers, fCut1_reco, fig11, rng);
    Cut2nue(mctruth, fRelTracks, fRelShowers, fCut2, fig11, rng);
    Cut3nue(mctruth, fRelTracks, fRelShowers, fCut3, fig11, rng);

/*    
    numuCCCut(mctruth, fRelTracks, fnumuCut, reco_compare, fmumu, fpie);
    numu_reco(mctracks, fnumu_reco1, fnumu_reco2);
*/
  }  
 
  fEventCounter++;
  
  return true;
}

  }  // namespace ExampleAnalysis
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::ExampleAnalysis::ExampleSelection)
