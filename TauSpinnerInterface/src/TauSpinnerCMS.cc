#include "TauSpinnerInterface/TauSpinnerInterface/interface/TauSpinnerCMS.h"

//MC-TESTER header files
#include "Tauola.h"
#include "tau_reweight_lib.h"
#include "Tauola_wrapper.h"
#include "TauSpinnerInterface/TauSpinnerInterface/interface/read_particles_from_HepMC.h"

bool TauSpinnerCMS::isTauSpinnerConfigure=false;

TauSpinnerCMS::TauSpinnerCMS( const ParameterSet& pset ) :
  isReco_(pset.getParameter<bool>("isReco"))
  ,isTauolaConfigured_(pset.getParameter<bool>("isTauolaConfigured" ))
  ,isLHPDFConfigured_(pset.getParameter<bool>("isLHPDFConfigured" ))
  ,LHAPDFname_(pset.getUntrackedParameter("LHAPDFname",(string)("MSTW2008nnlo90cl.LHgrid")))
  ,CMSEnergy_(pset.getParameter<double>("CMSEnergy"))//GeV
  ,gensrc_(pset.getParameter<edm::InputTag>("gensrc"))
  ,MotherPDGID_(pset.getUntrackedParameter("MotherPDGID",(int)(-1)))
{
  produces<double>("TauSpinnerWT").setBranchAlias("TauSpinnerWT");
  produces<double>("TauSpinnerWTFlip").setBranchAlias("TauSpinnerWTFlip");
  produces<double>("TauSpinnerWThplus").setBranchAlias("TauSpinnerWThplus");
  produces<double>("TauSpinnerWThminus").setBranchAlias("TauSpinnerWThminus");
}

void TauSpinnerCMS::beginJob()
{
  if(!isTauolaConfigured_){
    Tauolapp::Tauola::initialize();
  }
  if(!isLHPDFConfigured_){
    LHAPDF::initPDFSetByName(LHAPDFname_);   
  }
  if(!isTauSpinnerConfigure){
    isTauSpinnerConfigure=true;
    bool Ipp = true;  // for pp collisions 
    // Initialize TauSpinner
    //Ipol - polarization of input sample
    //nonSM2 - nonstandard model calculations
    //nonSMN
    int Ipol=0,nonSM2=0,nonSMN=0;
    TauSpinner::initialize_spinner(Ipp,Ipol,nonSM2,nonSMN,CMSEnergy_);
  }
}

void TauSpinnerCMS::produce( edm::Event& e, const edm::EventSetup& iSetup){
  double WT=1.0;
  double WTFlip=1.0;
  double polSM=-999; //range [-1,1]
  if(!e.isRealData()){
    double WT = 1.0; 
    SimpleParticle X, tau, tau2;
    vector<SimpleParticle> tau_daughters, tau_daughters2;
    int stat(0);
    if(isReco_){
      stat=readParticlesfromReco(e,X,tau,tau2,tau_daughters,tau_daughters2);
    }
    else{
      Handle< HepMCProduct > EvtHandle ;
      e.getByLabel( "generator", EvtHandle ) ;
      const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
      stat=readParticlesFromHepMC(Evt,X,tau,tau2,tau_daughters,tau_daughters2);
    }  
    if(MotherPDGID_<0 || abs(X.pdgid())==MotherPDGID_){
      if(stat!=1){
	// Determine the weight      
	if( abs(X.pdgid())==24 ||  abs(X.pdgid())==37 ){
	  WT = TauSpinner::calculateWeightFromParticlesWorHpn(X, tau, tau2, tau_daughters); // note that tau2 is tau neutrino
	  polSM=getTauSpin();
	  WTFlip=(2.0-WT)/WT;
	}
	else if( X.pdgid()==25 || X.pdgid()==36 || X.pdgid()==22 || X.pdgid()==23 ){
	  WT = TauSpinner::calculateWeightFromParticlesH(X, tau, tau2, tau_daughters,tau_daughters2);
	  polSM=getTauSpin();
	  if(X.pdgid()==25 || X.pdgid()==22 || X.pdgid()==23 ){
	    if(X.pdgid()==25) X.setPdgid(23);
	    if( X.pdgid()==22 || X.pdgid()==23 ) X.setPdgid(25);
	    double WTother=TauSpinner::calculateWeightFromParticlesH(X, tau, tau2, tau_daughters,tau_daughters2);
	    WTFlip=WTother/WT;
	  }
	}
	else{
	  cout<<"TauSpinner: WARNING: Unexpected PDG for tau mother: "<<X.pdgid()<<endl;
	}
      }
      
    }
  }
  // regular weight
  std::auto_ptr<double> TauSpinnerWeight(new double);
  *TauSpinnerWeight =WT;    
  e.put(TauSpinnerWeight,"TauSpinnerWT");  
  
  // flipped weight (ie Z->H or H->Z)
  std::auto_ptr<double> TauSpinnerWeightFlip(new double);
  *TauSpinnerWeightFlip =WTFlip;
  e.put(TauSpinnerWeightFlip,"TauSpinnerWTFlip");
  
  // h+ polarization
  double WThplus=WT;
  if(polSM<0.0 && polSM!=-999) WT=0; 
  std::auto_ptr<double> TauSpinnerWeighthplus(new double);
  *TauSpinnerWeighthplus = WThplus;
  e.put(TauSpinnerWeighthplus,"TauSpinnerWThplus");

  // h- polarization
  double WThminus=WT;
  if(polSM>0.0&& polSM!=-999) WT=0;
  std::auto_ptr<double> TauSpinnerWeighthminus(new double);
  *TauSpinnerWeighthminus = WThminus;
  e.put(TauSpinnerWeighthminus,"TauSpinnerWThminus");
  
  return ;
}  

void TauSpinnerCMS::endRun( const edm::Run& r, const edm::EventSetup& ){}

void TauSpinnerCMS::endJob(){}

int TauSpinnerCMS::readParticlesfromReco(edm::Event& e,SimpleParticle &X,SimpleParticle &tau,SimpleParticle &tau2, 
					 std::vector<SimpleParticle> &tau_daughters,std::vector<SimpleParticle> &tau2_daughters){
  edm::Handle<reco::GenParticleCollection> genParticles;
  e.getByLabel(gensrc_, genParticles);
  for(reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr!= genParticles->end(); ++itr){
    int pdgid=abs(itr->pdgId());
    if(pdgid==24 || pdgid==37 || pdgid ==25 || pdgid==36 || pdgid==22 || pdgid==23 ){
      const reco::GenParticle *recotau1=NULL;
      const reco::GenParticle *recotau2=NULL;
      unsigned int ntau(0),ntauornu(0);
      for(unsigned int i=0; i<itr->numberOfDaughters(); i++){
	const reco::Candidate *dau=itr->daughter(i);
	if(abs(dau->pdgId())!=pdgid){
	  if(abs(dau->pdgId())==15 || abs(dau->pdgId())==16){
	    if(ntau==0 && abs(dau->pdgId())==15){
	      recotau1=static_cast<const reco::GenParticle*>(dau);
	      ntau++;
	    }
	    else if((ntau==1 && abs(dau->pdgId())==15) || abs(dau->pdgId())==16){
	      recotau2=static_cast<const reco::GenParticle*>(dau);
	      if(abs(dau->pdgId())==15) ntau++;
	    }
	    ntauornu++;
	  }
	}
      }
      if((ntau==2 && ntauornu==2) || (ntau==1 && ntauornu==2)){
	X.setPx(itr->p4().Px());
	X.setPy(itr->p4().Py());
	X.setPz(itr->p4().Pz());
	X.setE (itr->p4().E());
	X.setPdgid(itr->pdgId());
	tau.setPx(recotau1->p4().Px());
	tau.setPy(recotau1->p4().Py());
	tau.setPz(recotau1->p4().Pz());
	tau.setE (recotau1->p4().E());
        tau.setPdgid(recotau1->pdgId());
	GetRecoDaughters(recotau1,tau_daughters,recotau1->pdgId());
	tau2.setPx(recotau2->p4().Px());
        tau2.setPy(recotau2->p4().Py());
        tau2.setPz(recotau2->p4().Pz());
        tau2.setE (recotau2->p4().E());
        tau2.setPdgid(recotau2->pdgId());
	if(ntau==2)GetRecoDaughters(recotau2,tau2_daughters,recotau2->pdgId());
	return 0;
      }
    }
  }
  return 1;
}


void TauSpinnerCMS::GetRecoDaughters(const reco::GenParticle *Particle,std::vector<SimpleParticle> &daughters, int parentpdgid){
  if(Particle->pdgId()!=parentpdgid){
    SimpleParticle tp(Particle->p4().Px(), Particle->p4().Py(), Particle->p4().Pz(), Particle->p4().E(), Particle->pdgId());
    daughters.push_back(tp);
  }
  for (unsigned int i=0; i< Particle->numberOfDaughters(); i++){
    const reco::Candidate *dau=Particle->daughter(i);
    if(abs(Particle->pdgId())!=111)GetRecoDaughters(static_cast<const reco::GenParticle*>(dau),daughters,Particle->pdgId());
  }
}

DEFINE_FWK_MODULE(TauSpinnerCMS);
