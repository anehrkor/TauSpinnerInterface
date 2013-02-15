#include "TauSpinnerInterface/TauSpinnerInterface/interface/TauSpinnerCMS.h"

//MC-TESTER header files
#include "Tauola.h"
#include "tau_reweight_lib.h"
#include "Tauola_wrapper.h"
#include "TauSpinnerInterface/TauSpinnerInterface/interface/read_particles_from_HepMC.h"

TauSpinnerCMS::TauSpinnerCMS( const ParameterSet& pset ) :
  isReco_(pset.getUntrackedParameter("isReco",(bool)(false)))
  ,LHAPDFname_(pset.getUntrackedParameter("LHAPDFname",(string)("MSTW2008nnlo90cl.LHgrid")))
  ,gensrc_(pset.getParameter<edm::InputTag>( "gensrc" ))
{
  produces<double>("TauSpinerWT").setBranchAlias("TauSpinerWT");
}

void TauSpinnerCMS::beginJob()
{
    
  if(isReco_){
    Tauolapp::Tauola::initialize();
    string name="MSTW2008nnlo90cl.LHgrid";
    LHAPDF::initPDFSetByName(name);   
  }
  double CMSENE = 8000.0; // center of mass system energy.
                          // used in PDF calculation. For pp collisions only
  bool Ipp = true;  // for pp collisions 
  // Initialize TauSpinner
  //Ipol - polarization of input sample
  //nonSM2 - nonstandard model calculations
  //nonSMN
  int Ipol=0,nonSM2=0,nonSMN=0;
  TauSpinner::initialize_spinner(Ipp,Ipol,nonSM2,nonSMN,CMSENE);

  return ;
  
}

void TauSpinnerCMS::produce( edm::Event& e, const edm::EventSetup& iSetup)
{
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

    if(stat!=1){
      // Determine the weight      

      if( abs(X.pdgid())==24 ||  abs(X.pdgid())==37 ){
	WT = TauSpinner::calculateWeightFromParticlesWorHpn(X, tau, tau2, tau_daughters); // note that tau2 is tau neutrino
      }
      else if( X.pdgid()==25 || X.pdgid()==36 || X.pdgid()==22 || X.pdgid()==23 ){
	WT = TauSpinner::calculateWeightFromParticlesH(X, tau, tau2, tau_daughters,tau_daughters2);
      }
      else{
	cout<<"TauSpiner: WARNING: Unexpected PDG for tau mother: "<<X.pdgid()<<endl;
      }
    }
    
    cout<<"Stat: "<<stat<<" WT: "<<WT<<endl;
    
    std::auto_ptr<double> TauSpinerWeight(new double);
    *TauSpinerWeight =WT;    
    
    e.put(TauSpinerWeight,"TauSpinerWT");  
  }
  return ;
}  

void TauSpinnerCMS::endRun( const edm::Run& r, const edm::EventSetup& )
{

   return;

}


void TauSpinnerCMS::endJob()
{

   return ;
}



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
