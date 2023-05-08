//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================

// Framework include files
#include "DualCrystalCalorimeterHit.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"
#include <cmath>








namespace CalVision {

  G4double fromEvToNm(G4double energy)
  {
    // there is some bug somewhere.  shouldn't need this factor
    return 1239.84187 / energy*1000.;
  }


  int SCECOUNT=0;

  class DualCrystalCalorimeterSD {
  public:
    typedef DualCrystalCalorimeterHit Hit;
    // If we need special data to personalize the action, be put it here
    //int mumDeposits = 0;
    //double integratedDeposit = 0;
  };
}

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {
  /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
  namespace sim   {

    using namespace CalVision;
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //               Geant4SensitiveAction<MyTrackerSD>
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /** \addtogroup Geant4SDActionPlugin
     *
     * @{
     * \package DualCrystalCalorimeterSDAction
     *
     * @}
     */

    /// Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<DualCrystalCalorimeterSD>::defineCollections()    {
      m_collectionID = declareReadoutFilteredCollection<CalVision::DualCrystalCalorimeterSD::Hit>();
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <> bool Geant4SensitiveAction<DualCrystalCalorimeterSD>::process(const G4Step* step,G4TouchableHistory* /*hist*/ ) {


      bool SCEPRINT=(SCECOUNT<10000);
      //if(SCEPRINT) std::cout<<"scecount is "<<SCECOUNT<<" print is "<<SCEPRINT<<std::endl;



      G4StepPoint *thePrePoint = step->GetPreStepPoint();
      G4StepPoint *thePostPoint = step->GetPostStepPoint();
      //      const G4ThreeVector &thePrePosition = thePrePoint->GetPosition();
      //const G4ThreeVector &thePostPosition = thePostPoint->GetPosition();
      G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
      G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();
      G4String thePrePVName = "";
      if (thePrePV)
	thePrePVName = thePrePV->GetName();
      G4String thePostPVName = "";
      if (thePostPV)
	thePostPVName = thePostPV->GetName();
      //G4Track *theTrack = step->GetTrack();
      //G4int TrPDGid = theTrack->GetDefinition()->GetPDGEncoding();

      //      if(thePrePVName.contains("slice")==0) {
      //std::cout<<"entering DualCrystalAction"<<std::endl;
      //  std::cout<<" prevolume is "<<thePrePVName<<std::endl;
      //  std::cout<<" postvolume is "<<thePostPVName<<std::endl;
      //  std::cout<<" pid is "<<TrPDGid<<std::endl;
	  //}


      Geant4StepHandler h(step);
      HitContribution contrib = DualCrystalCalorimeterHit::extractContribution(step);

      Geant4HitCollection*  coll    = collection(m_collectionID);
      VolumeID cell = 0;

      try {
        cell = cellID(step);
      } catch(std::runtime_error &e) {
	std::stringstream out;
        out << std::setprecision(20) << std::scientific;
        out << "ERROR: " << e.what()  << std::endl;
        out << "Position: "
            << "Pre (" << std::setw(24) << step->GetPreStepPoint()->GetPosition() << ") "
            << "Post (" << std::setw(24) << step->GetPostStepPoint()->GetPosition() << ") "
            << std::endl;
        out << "Momentum: "
            << " Pre (" <<std::setw(24) << step->GetPreStepPoint() ->GetMomentum()  << ") "
            << " Post (" <<std::setw(24) << step->GetPostStepPoint()->GetMomentum() << ") "
            << std::endl;

	std::cout << out.str();

        return true;
      }


      DualCrystalCalorimeterHit* hit = coll->findByKey<DualCrystalCalorimeterHit>(cell);
      if ( !hit ) {
        Geant4TouchableHandler handler(step);
	DDSegmentation::Vector3D pos = m_segmentation.position(cell);
        Position global = h.localToGlobal(pos);
        hit = new DualCrystalCalorimeterHit(global);
        hit->cellID = cell;
        coll->add(cell, hit);
        printM2("CREATE hit with deposit:%e MeV  Pos:%8.2f %8.2f %8.2f  %s",
                contrib.deposit,pos.X,pos.Y,pos.Z,handler.path().c_str());
	std::cout<<"DRcalo deposit "<<contrib.deposit<<" position ("<<pos.X<<","<<pos.Y<<","<<pos.Z<<") string "<<handler.path().c_str()<<" volume id "<<cell<<std::endl;

        if ( 0 == hit->cellID )  { // for debugging only!
          hit->cellID = cellID(step);
          except("+++ Invalid CELL ID for hit!");
        }
      } else {
	//	std::cout<<"updating old hit"<<std::endl;
      }


      G4Track * track =  step->GetTrack();

      if(SCEPRINT) std::cout<< (track->GetDefinition())->GetParticleName()<<std::endl;

      //photons
      if( track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() )  {
	if(SCEPRINT) std::cout<<"     in volume ID "<<cell<<std::endl;

	SCECOUNT+=1;
	//	if(SCEPRINT) std::cout<<"optical photon"<<std::endl;

	bool OptN = (track->GetCreatorProcess()->G4VProcess::GetProcessName() == "CerenkovPhys")||(track->GetCreatorProcess()->G4VProcess::GetProcessName() == "ScintillationPhys");

	//if(track->GetParentID()!=1) SCEPRINT=1;
	if( (track->GetCreatorProcess()->G4VProcess::GetProcessName() != "CerenkovPhys")&&(track->GetCreatorProcess()->G4VProcess::GetProcessName() != "ScintillationPhys") ) SCEPRINT=1;

	float wavelength=fromEvToNm(track->GetTotalEnergy()/eV);
	int ibin=-1;
	float binsize = (hit->wavelenmax-hit->wavelenmin)/hit->nbin;
	ibin = (wavelength-hit->wavelenmin)/binsize;

  /*float tof=step->GetPreStepPoint()->GetGlobalTime();
  int ibin_t=-1;
  float binsize_t = (hit->timemax-hit->timemin)/hit->nbin_t;
  ibin_t = (tof-hit->timemin)/binsize_t;*/

  float poi_x,poi_y,poi_z; //Pre step if obtained in the filter will be for the previous silicone grease gap
  poi_x = step->GetPreStepPoint()->GetMomentum().x();
  poi_y = step->GetPreStepPoint()->GetMomentum().y();
  poi_z = step->GetPreStepPoint()->GetMomentum().z();
  
  float nrm_x_long = 0., nrm_y_long = 1., nrm_z_long = 0.; //y direction where the silicone_gap - filter interface is
  float nrm_x_right = 0.,  nrm_y_right = 0., nrm_z_right = 1.; //z direction i.e. the right side crystal-silicone gap interface
  float aoi_long, aoi_right;
	aoi_long=180*acos((poi_x*nrm_x_long + poi_y*nrm_y_long + poi_z*nrm_z_long)/(sqrt(pow(poi_x,2) + pow(poi_y,2) + pow(poi_z,2))*sqrt(pow(nrm_x_long,2) + pow(nrm_y_long,2) + pow(nrm_z_long,2))))/(atan(1)*4);
 	aoi_right=180*acos((poi_x*nrm_x_right + poi_y*nrm_y_right + poi_z*nrm_z_right)/(sqrt(pow(poi_x,2) + pow(poi_y,2) + pow(poi_z,2))*sqrt(pow(nrm_x_right,2) + pow(nrm_y_right,2) + pow(nrm_z_right,2))))/(atan(1)*4);

	int ibin_a=-1;
	float binsize_a = (hit->anglemax-hit->anglemin)/hit->nbin_a;
   
  int phstep = track->GetCurrentStepNumber();
	

	if ( track->GetCreatorProcess()->G4VProcess::GetProcessName() == "CerenkovPhys")
  {
	  if(SCEPRINT) std::cout<<" found cerenkov photon"<<std::endl;
	  std::string amedia = ((track->GetMaterial())->GetName());
    if(amedia.find("kill")!=std::string::npos)     
	  //if(((track->GetMaterial())->GetName())=="killMedia") 
	    { 
	      if(SCEPRINT) std::cout<<"killing photon"<<std::endl;
        if(phstep>1)// don't count photons created in kill media i.e. with step number 1 because they are unphysical
       {
	      hit->ncerenkov+=1;
	      if(ibin>-1&&ibin<hit->nbin) ((hit->ncerwave).at(ibin))+=1;
        //if(ibin_t>-1&&ibin_t<hit->nbin_t) ((hit->ncertime).at(ibin_t))+=1; //Time of arrival at killMedia for Cerenkov
       }
	      track->SetTrackStatus(fStopAndKill); //kill ALL killMedia photons after the first step, regardless of whether they were created there or not (to avoid overcounting)
     }
	  
    else {
	    //if( (track->GetParentID()==1)&&(track->GetCurrentStepNumber()==1)  ) hit->ncerenkov+=1;
         	    if( (phstep==1)  )
              {
               hit->ncerenkov+=1;
               if(ibin>-1&&ibin<hit->nbin) ((hit->ncerwave).at(ibin))+=1;
               //if(ibin_t>-1&&ibin_t<hit->nbin_t) ((hit->ncertime).at(ibin_t))+=1; //Time of production (mainly in crystal) for Cerenkov
              }
	  }
           //if(((track->GetMaterial())->GetName())=="Air" && step->GetPreStepPoint()->GetPosition() != step->GetPostStepPoint()->GetPosition())
          if(/*((track->GetMaterial())->GetName())=="Air"*/amedia.find("Air")!=std::string::npos && (fabs(step->GetPreStepPoint()->GetPosition().x() - step->GetPostStepPoint()->GetPosition().x())>0.1 || fabs(step->GetPreStepPoint()->GetPosition().y() - step->GetPostStepPoint()->GetPosition().y())>0.1 || fabs(step->GetPreStepPoint()->GetPosition().z() - step->GetPostStepPoint()->GetPosition().z())>0.1))
          {
           hit->ncerenkov+=1;
           //track->SetTrackStatus(fStopAndKill); //So all the photons go to air including the ones that suffer at least one TIR and are counted in Air and aborted
           //Only count, no abort
           //std::cout<<"escaped cerenkov photon"<<std::endl;
          }
         //return false;
         
         if(amedia.find("OD2")!=std::string::npos && step->GetPreStepPoint()->GetPosition().y() == 12.5/*13.5*/) //interference filter OD2_600
         {
        	ibin_a = (aoi_long-hit->anglemin)/binsize_a;
          if(ibin_a>-1&&ibin_a<hit->nbin_a) ((hit->ncerangle).at(ibin_a))+=1; //Angle of incidence (at filter) for Cerenkov
         }
         
         else if(amedia.find("Silicone_gap_right")!=std::string::npos && step->GetPreStepPoint()->GetPosition().z() == 30.)
         {
        	ibin_a = (aoi_right-hit->anglemin)/binsize_a;
          if(ibin_a>-1&&ibin_a<hit->nbin_a) ((hit->ncerangle).at(ibin_a))+=1; //Angle of incidence (at right silicone interface) for Cerenkov         
         }
        }
         
	else if (  track->GetCreatorProcess()->G4VProcess::GetProcessName() == "ScintillationPhys"  )
  {
     if(SCEPRINT) std::cout<<" found scintillation photon"<<std::endl;
     std::string amedia = ((track->GetMaterial())->GetName());
     if(amedia.find("kill")!=std::string::npos)
     //if(((track->GetMaterial())->GetName())=="killMedia") 
	    {
	      if(SCEPRINT) std::cout<<"killing photon"<<std::endl;
	      if(phstep>1) // don't count photons created in kill media i.e. with step number 1 because they are unphysical      
        {
	       hit->nscintillator+=1;
	       if((ibin>-1)&&(ibin<hit->nbin)) ((hit->nscintwave).at(ibin))+=1;
         //if(ibin_t>-1&&ibin_t<hit->nbin_t) ((hit->nscinttime).at(ibin_t))+=1; //Time of arrival at killMedia for scintillation 
        }
       track->SetTrackStatus(fStopAndKill); //kill ALL killMedia photons after the first step, regardless of whether they were created there or not (to avoid overcounting)
       }
        
	  else {
	    //if( (track->GetParentID()==1)&&(track->GetCurrentStepNumber()==1) ) hit->nscintillator+=1; 
         	    if( (phstep==1) ) 
              {
               hit->nscintillator+=1;
               if((ibin>-1)&&(ibin<hit->nbin)) ((hit->nscintwave).at(ibin))+=1;
               //if(ibin_t>-1&&ibin_t<hit->nbin_t) ((hit->nscinttime).at(ibin_t))+=1; //Time of production (mainly in crystal) for scintillation
              }
	  }
          //if(((track->GetMaterial())->GetName())=="Air" && step->GetPreStepPoint()->GetPosition() != step->GetPostStepPoint()->GetPosition())
          if(/*((track->GetMaterial())->GetName())=="Air"*/amedia.find("Air")!=std::string::npos && (fabs(step->GetPreStepPoint()->GetPosition().x() - step->GetPostStepPoint()->GetPosition().x())>0.1 || fabs(step->GetPreStepPoint()->GetPosition().y() - step->GetPostStepPoint()->GetPosition().y())>0.1 || fabs(step->GetPreStepPoint()->GetPosition().z() - step->GetPostStepPoint()->GetPosition().z())>0.1))

          {
           hit->nscintillator+=1;
           //std::cout<<"escaped scintillator photon"<<std::endl;
          }
          //return false;
          if(amedia.find("OD2")!=std::string::npos && step->GetPreStepPoint()->GetPosition().y() == 12.5/*13.5*/) //interference filter OD2_600
         {
        	ibin_a = (aoi_long-hit->anglemin)/binsize_a;         
          if(ibin_a>-1&&ibin_a<hit->nbin_a) ((hit->nscintangle).at(ibin_a))+=1; //Angle of incidence (at filter) for Scintillation
         }
         else if(amedia.find("Silicone_gap_right")!=std::string::npos && step->GetPreStepPoint()->GetPosition().z() == 30.)
         {
        	ibin_a = (aoi_right-hit->anglemin)/binsize_a;
          if(ibin_a>-1&&ibin_a<hit->nbin_a) ((hit->nscintangle).at(ibin_a))+=1; //Angle of incidence (at right silicone interface) for Scintillation         
         }
        }
	else {
          //if(SCEPRINT) std::cout<<"      other photon"<<std::endl;
          //track->SetTrackStatus(fStopAndKill);
          //return false;
	}
 
 	if(SCEPRINT) {
	  std::cout<<"     SCECOUNT="<<SCECOUNT<<std::endl;
	
	  std::cout<<"     will robinson have photon "<<track->GetCreatorProcess()->G4VProcess::GetProcessName() <<std::endl;
	  std::cout<<"     photon mother is "<<track->GetParentID()<<std::endl;
	  std::cout<<"     photon material is "<<(track->GetMaterial())->GetName()<<std::endl;
	  std::cout<<"     photon creator process is "<<(track->GetCreatorProcess())->GetProcessName()<<std::endl;
	  std::cout<<"     photon  process  type is "<<(track->GetCreatorProcess())->GetProcessType()<<std::endl;
	  std::cout<<"     photon sub process is "<<(track->GetCreatorProcess())->GetProcessSubType()<<std::endl;
	  std::cout<<"     photon current step number is "<<track->GetCurrentStepNumber()<<std::endl;
    std::cout<<"     photon pre position is: x coordinate- "<<step->GetPreStepPoint()->GetPosition().x() <<" y coordinate- "<<step->GetPreStepPoint()->GetPosition().y()<<" z coordinate- "<<step->GetPreStepPoint()->GetPosition().z()<<std::endl;
    std::cout<<"     photon post position is: x coordinate- "<<step->GetPostStepPoint()->GetPosition().x() <<" y coordinate- "<<step->GetPostStepPoint()->GetPosition().y() <<" z coordinate- "<<step->GetPostStepPoint()->GetPosition().z()<<std::endl;
    
    std::cout<<"     photon pre momentum in eV is: x coordinate- "<<(step->GetPreStepPoint()->GetMomentum().x())*pow(10.0,6) <<" y coordinate- "<<(step->GetPreStepPoint()->GetMomentum().y())*pow(10.0,6) <<" z coordinate- "<< (step->GetPreStepPoint()->GetMomentum().z())*pow(10.0,6) <<std::endl;
    std::cout<<"     photon post momentum in eV is: x coordinate- "<<(step->GetPostStepPoint()->GetMomentum().x())*pow(10.0,6) <<" y coordinate- "<<(step->GetPostStepPoint()->GetMomentum().y())*pow(10.0,6) <<" z coordinate- "<<(step->GetPostStepPoint()->GetMomentum().z())*pow(10.0,6) <<std::endl;
    
    
    std::cout<<"     photon angle of incidence at long face is "<<180*acos((step->GetPreStepPoint()->GetMomentum().x()*nrm_x_long + step->GetPreStepPoint()->GetMomentum().y()*nrm_y_long + step->GetPreStepPoint()->GetMomentum().z()*nrm_z_long)/(sqrt(pow(step->GetPreStepPoint()->GetMomentum().x(),2) + pow(step->GetPreStepPoint()->GetMomentum().y(),2) + pow(step->GetPreStepPoint()->GetMomentum().z(),2))*sqrt(pow(nrm_x_long,2) + pow(nrm_y_long,2) + pow(nrm_z_long,2))))/(atan(1)*4)<<std::endl;
    
    std::cout<<"     photon angle of incidence at right face is "<<180*acos((step->GetPreStepPoint()->GetMomentum().x()*nrm_x_right + step->GetPreStepPoint()->GetMomentum().y()*nrm_y_right + step->GetPreStepPoint()->GetMomentum().z()*nrm_z_right)/(sqrt(pow(step->GetPreStepPoint()->GetMomentum().x(),2) + pow(step->GetPreStepPoint()->GetMomentum().y(),2) + pow(step->GetPreStepPoint()->GetMomentum().z(),2))*sqrt(pow(nrm_x_right,2) + pow(nrm_y_right,2) + pow(nrm_z_right,2))))/(atan(1)*4)<<std::endl;  
    
 	  std::cout<<"     the pre volume name is "<<thePrePVName<<std::endl;
	  std::cout<<"     the post volume name is "<<thePostPVName<<std::endl;
	//(track->GetCreatorProcess())->DumpInfo();
	  std::cout<<"     photon energy is "<<track->GetTotalEnergy()/eV<<std::endl;
	  std::cout<<"     photon wavelength is "<<fromEvToNm(track->GetTotalEnergy()/eV)<<std::endl;
	  std::cout<<"     number of cherenkov is "<<hit->ncerenkov<<std::endl;
	  std::cout<<"     number of scintillation is "<<hit->nscintillator<<std::endl;
	}


      }

    else {   // particles other than optical photons
	
      if(SCEPRINT) std::cout<<"NOT optical photon"<<std::endl;

        hit->energyDeposit += contrib.deposit;
        hit->truth.emplace_back(contrib);

        //return true;
      }

        mark(h.track);
	
      return true;

    }
	
  }
} // end namespace calvision







namespace dd4hep { namespace sim {

    using namespace CalVision;

    struct WavelengthMinimumCut : public dd4hep::sim::Geant4Filter  {
  /// Energy cut value
      double m_wavelengthCut;
    public:
  /// Constructor.
      WavelengthMinimumCut(dd4hep::sim::Geant4Context* c, const std::string& n);
  /// Standard destructor
      virtual ~WavelengthMinimumCut();
  /// Filter action. Return true if hits should be processed
      virtual bool operator()(const G4Step* step) const  override  final  {
	bool test=true;
	G4Track *theTrack = step->GetTrack();
	if(theTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) {
	  float energy=theTrack->GetTotalEnergy()/eV;
	  float wave=fromEvToNm(energy);
	  if(wave < m_wavelengthCut) {
	    theTrack->SetTrackStatus(fStopAndKill);
	    test=false;}
	}
	return test;
      }
      virtual bool operator()(const Geant4FastSimSpot* spot) const  override  final  {
	return true;
      }
    };

  /// Constructor.
    WavelengthMinimumCut::WavelengthMinimumCut(Geant4Context* c, const std::string& n)
      : Geant4Filter(c,n) {
      InstanceCount::increment(this);
      declareProperty("Cut",m_wavelengthCut=0.0);
    }

  /// Standard destructor
    WavelengthMinimumCut::~WavelengthMinimumCut() {
      InstanceCount::decrement(this);
    }



    struct WavelengthnmwindCut : public dd4hep::sim::Geant4Filter  {
  /// Energy cut value
      double m_wavelengthstart;
    public:
  /// Constructor.
      WavelengthnmwindCut(dd4hep::sim::Geant4Context* c, const std::string& n);
  /// Standard destructor
      virtual ~WavelengthnmwindCut();
  /// Filter action. Return true if hits should be processed
      virtual bool operator()(const G4Step* step) const  override  final  {
	bool test=true;
	G4Track *theTrack = step->GetTrack();
	if(theTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) {
	  float energy=theTrack->GetTotalEnergy()/eV;
	  float wave=fromEvToNm(energy);
	  if((wave<m_wavelengthstart) || (wave > m_wavelengthstart+0.5) ) {
	    theTrack->SetTrackStatus(fStopAndKill);
	    test=false;}
	}
	return test;
      }
      virtual bool operator()(const Geant4FastSimSpot* spot) const  override  final  {
	return true;
      }
    };

  /// Constructor.
    WavelengthnmwindCut::WavelengthnmwindCut(Geant4Context* c, const std::string& n)
      : Geant4Filter(c,n) {
      InstanceCount::increment(this);
      declareProperty("Cut",m_wavelengthstart=0.0);
    }

  /// Standard destructor
    WavelengthnmwindCut::~WavelengthnmwindCut() {
      InstanceCount::decrement(this);
    }


  }}  // end using namespace





//--- Factory declaration
namespace dd4hep { namespace sim {
    typedef Geant4SensitiveAction<DualCrystalCalorimeterSD> DualCrystalCalorimeterSDAction;
  }}
DECLARE_GEANT4SENSITIVE(DualCrystalCalorimeterSDAction)
DECLARE_GEANT4ACTION(WavelengthMinimumCut)
DECLARE_GEANT4ACTION(WavelengthnmwindCut)
