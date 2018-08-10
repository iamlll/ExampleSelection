#ifndef __sbnanalysis_ana_ExampleAnalysis_ExampleSelection__
#define __sbnanalysis_ana_ExampleAnalysis_ExampleSelection__

/**
 * \file ExampleSelection.h
 *
 * An example event selection processor.
 *
 * This is an implementation of the core::SelectionBase class. We define
 * the methods called for initialization, finalization, and event-by-event
 * processing.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>
 */

// Includes
#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include <fstream>
#include "lardataobj/MCBase/MCShower.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"
#include <random>

// Forward declarations
class TH2D;
class TH1D;
class THStack;
class TLorentzVector;

/** All analysis code is defined in namespace "ana" */
namespace ana {

  /** Code specific to the ExampleAnalysis. */
  namespace ExampleAnalysis {

/**
 * \class ExampleSelection
 * \brief An example selection analysis
 *
 * This selection analysis doesn't actually select events, it just
 * demonstrates the framework!
 */
class ExampleSelection : public core::SelectionBase {
public:
  /** Constructor. */
  ExampleSelection();

  /**
   * Initialization.
   *
   * Here we load configuration parameters, set up histograms for output, and
   * add our own branches to the output tree.
   *
   * \param config A configuration, as a JSON object
   */
  void Initialize(Json::Value* config=NULL);

  /** Finalize and write objects to the output file. */
  void Finalize();

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \return True to keep event
   */
  bool ProcessEvent(gallery::Event& ev);

protected:
  unsigned fEventCounter;  //!< Count processed events

  /** Configuration parameters */
  std::mt19937 rng;

  art::InputTag fTruthTag;  //!< art tag for MCTruth information
  art::InputTag fGenieTag; //art tag for event weight info (MCEventWeight)
  art::InputTag fFluxTag; //art tag for flux weight info (MCEventWeight)
  art::InputTag fParticleTag; //MCParticle
  art::InputTag fTrackTag; //MCShower + MCTrack
  int fMyParam;  //!< A parameter from the configuration file

  /** Custom data branches */
  int fNuCount;  //!< Number of neutrino interactions in the event
  int fMyVar;  //!< Another variable of interest
  std::vector<int> fInteraction; //interaction mode
  std::vector<int> fPParticle; //hadron parent particles
  std::vector<double> fEnergy; //neutrino interaction energy

  std::vector<sim::MCTrack> fRelTracks; //MCTracks within 5 cm of neutrino interaction vertex
  std::vector<sim::MCShower> fRelShowers; //MCShowers within 5 cm of neutrino interaction vertex
  std::vector<int>  fTrackCodes; 
  std::vector<int> fShowerCodes;

  std::vector<std::vector<double>> fEventWt; 
  std::vector<std::vector<double>> fTotFluxWt;
  std::vector<double> fWts;

  /** Histograms */
  TH2D* fNuVertexXZHist;  //!< Neutrino vertex XZ projection
  std::vector<TH1D*> fCut1; //cuts for selection criteria nu_e CC electron with E_e > 200 MeV
  THStack* fnueCC;
  std::vector<TH1D*> fCut1_reco; //Cut 1 with reconstructed neutrino energy
  THStack* fnueCC_reco;
  std::vector<TH1D*> fig11;
  THStack* fig11Stack;

  std::vector<TH1D*> fCut2;
  THStack* fGammaNC;
  std::vector<TH1D*> fCut3;
  THStack* fnumuCC;

  std::vector<TH1D*> fnumuCut;
  THStack* fnumuCandidates; 
  std::vector<TH1D*> fmumu; //0 is total muons; 1 is muons passing geometric constraints (begin inside detector, track length sufficiently long)
  std::vector<TH1D*> fpie; //0 is total #pi0s; 1 is pi0s that passed the length cut, i.e. mistaken for muon CC events only based on length cut; 2 is total num antiparticles that passed
};

  int nummu;
  int numpi;
  int nummupass;
  int numpipass;

  //ALL of my variable cut histograms
  std::ofstream familyConfusion;
  TH1D* nuFlux; //neutrino count vs energy
  std::vector<TH1D*> E_gamma; //ct vs photon shower energy
  THStack* gammaStack;
  TH2D* E_shower_nu; //shower energy vs neutrino energy
  TH1D* E_shower; //ct vs shower energy
  TH1D* dEdx;
  std::vector<TH1D*> dEdx_2; //differentiates b/w gamma and e- showers
  THStack* dEdx_2_stack;
  TH1D* dist_from_vertex; //ct vs distance from neutrino vertex
  TH2D* E_hadron; //count vs hadron energies 
  std::vector<TH1D*> hPass; //hadrons that pass the gamma (energy) cuts
  THStack* hadStack;
  std::vector<TH1D*> pPass; //interactions that have a visible vertex  showers that have     
  THStack* phoStack;

  TH1D* fnumu_reco1; 
  std::vector<TH1D*> fnumu_reco2; //distinguish b/w PDG codes
  THStack* fnumu_reco_stack;
  std::vector<TH1D*> reco_compare;
  THStack* reco_compare_stack;
  THStack* fDiscriminate_numu_mupi; 
  }  // namespace ExampleAnalysis
}  // namespace ana

#endif  // __sbnanalysis_ana_ExampleAnalysis_ExampleSelection__



