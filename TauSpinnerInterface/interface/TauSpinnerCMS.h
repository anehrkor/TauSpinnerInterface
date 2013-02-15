// -*- C++ -*-
//
// Package:    TauSpinnerInterface
// Class:      TauSpinnerCMS
// 
/**\class TauSpinnerCMS TauSpinnerCMS.cc 

*/
//
// Original Author:  Ian Nugent  
//         Created:  Fri Feb 15 2013

#ifndef TauSpinnerCMS_h
#define TauSpinnerCMS_h


#include <iostream>

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

// essentials !!!
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimpleParticle.h"


using namespace edm;
using namespace std;
using namespace TauSpinner;

class TauSpinnerCMS : public edm::EDProducer
{
  
 public:
  
  //
  explicit TauSpinnerCMS( const edm::ParameterSet& ) ;
  virtual ~TauSpinnerCMS() {} // no need to delete ROOT stuff
  // as it'll be deleted upon closing TFile
  
  virtual void produce( edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endRun( const edm::Run&, const edm::EventSetup& ) ;
  virtual void endJob() ;
  
 private:
  bool isReco_;
  string LHAPDFname_; 
  edm::InputTag gensrc_;
  int readParticlesfromReco(edm::Event& e,SimpleParticle &X,SimpleParticle &tau,SimpleParticle &tau2,
			    std::vector<SimpleParticle> &tau_daughters,std::vector<SimpleParticle> &tau2_daughters);
  void GetRecoDaughters(const reco::GenParticle *Particle,std::vector<SimpleParticle> &daughters,int parentpdgid);

}; 
#endif
