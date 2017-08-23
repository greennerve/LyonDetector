#include "G4VProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

#include "LyonSensitiveDetector.hh"

LyonSensitiveDetector::LyonSensitiveDetector(G4String name)
  : G4VSensitiveDetector(name),
    collectionID(-1)
{
  collectionName.insert("LyonTrackHits");
  efficiency = 0.95;
  precision = 5*CLHEP::mm/sqrt(12);
  fRandomGauss = new CLHEP::RandGauss( CLHEP::HepRandom::getTheEngine() );
  stdDev = 1.7*cm/sqrt(12);
}

LyonSensitiveDetector::~LyonSensitiveDetector()
{
  //  delete fRandomGauss;
}
void LyonSensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
  if (collectionID < 0) 
    collectionID = GetCollectionID(0);
    
  hitsCollection = new LyonTrackHitsCollection(
					       SensitiveDetectorName,
					       collectionName[0]);
    
  HCE->AddHitsCollection(
			 collectionID,
			 hitsCollection);

  NbOfTracks = -1;
  CalID.clear();
  RPCID.clear();
}


G4bool LyonSensitiveDetector::ProcessHits(G4Step* currentStep,G4TouchableHistory*)
{ 

  if(NbOfTracks<0)
    {
      event = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
      NbOfTracks = event->GetNumberOfPrimaryVertex();
      
      CalID.resize(NbOfTracks);
      RPCID.resize(NbOfTracks);
      
      for(G4int i = 0; i<NbOfTracks; i++) 
	{
	  CalID.push_back(-1);
	  RPCID.push_back(-1);
	}
    }

  G4int TID = currentStep->GetTrack()->GetTrackID();    

  if(TID > NbOfTracks)
    return false;
  
  if( G4UniformRand() < efficiency)
    {
      
      G4TouchableHistory * theTouchable = 
	(G4TouchableHistory *) currentStep->GetPreStepPoint()->GetTouchable();

      G4int fChamber = theTouchable->GetReplicaNumber(1);
      G4int fCalo = theTouchable->GetReplicaNumber(2);
          
      G4int TrackNb = TID - 1; //TrackID starts from 1, track number starts from 0;


      if(CalID[TrackNb] != fCalo || RPCID[TrackNb] != fChamber)//To avoid collecting multiple hits in the same layer
	{
	  LyonTrackHit* aHit= new LyonTrackHit();
	  G4ThreeVector ExitPoint = currentStep->GetPostStepPoint()->GetPosition();
	  G4int fB = theTouchable->GetReplicaNumber(0);
	  if(fCalo == 0 || fCalo == 1)
	    {
	      ExitPoint += G4ThreeVector(fRandomGauss->fire(0,stdDev),
					 fRandomGauss->fire(0,stdDev),
					 0.0);
	    }
	  else
	    {
	      if(fCalo == 2 || fCalo == 3)
		{
		  ExitPoint += G4ThreeVector(0.0,
					     fRandomGauss->fire(0,stdDev),
					     fRandomGauss->fire(0,stdDev));
		}
	      else
		{
		  ExitPoint += G4ThreeVector(fRandomGauss->fire(0,stdDev),
					     0.0,					     
					     fRandomGauss->fire(0,stdDev));
		}

	    }
	  aHit->SetExitPoint(ExitPoint);
	  aHit->SetChamber(fChamber);
	  aHit->SetCalo(fCalo);
	  aHit->SetBasicVolumeID(fB);
	  aHit->SetTrackID(TID);
      
	  hitsCollection->insert(aHit);
	  CalID[TrackNb] = fCalo;
	  RPCID[TrackNb] = fChamber;

	  G4cout<<"TrackID: "<<TID
		<<"; CalID: "<<fCalo
		<<"; RPCID: "<<fChamber<<G4endl;
	  return true;
	}
      else
	return false;
    }
  
  else
    return false;
}

/*
  if(currentStep->GetTrack()->GetTrackID() > NbOfTracks)
  return true;
  
  if( G4UniformRand() < efficiency)
  {
  G4ThreeVector precisionRandomVector(
  0.0,
  0.0,
  0.0);
  LyonTrackHit* aHit= new LyonTrackHit(currentStep, precisionRandomVector);
  G4int TrackNb = aHit->GetTrackID() - 1; //TrackID starts from 1, track number starts from 0;
  G4int Calid = aHit->GetCalo();
  G4int RPCid = aHit->GetChamber();
  if(CalID[TrackNb] != Calid || RPCID[TrackNb] != RPCid)//To avoid collecting multiple hits in the same layer
  {
  hitsCollection->insert(aHit);
  CalID[TrackNb] = Calid;
  RPCID[TrackNb] = RPCid;

  G4cout<<"TrackID: "<<aHit->GetTrackID()
  <<"; CalID: "<<Calid
  <<"; RPCID: "<<RPCid<<G4endl;
  }
  else
  {
  delete aHit;
  }
  }
  
  return true;
*/

//....ooOOO0OOOoo........ooOOO0OOOoo........ooOOO0OOOoo........ooOOO0OOOoo....
//a temperary SD

/*
  LyonAirSensitiveDetector::LyonAirSensitiveDetector(G4String name)
  : G4VSensitiveDetector(name),
  collectionID(-1)
  {
  collectionName.insert("LyonAirHits");
  }

  void LyonAirSensitiveDetector::Initialize(G4HCofThisEvent* HCE)
  {
  if (collectionID < 0) 
  collectionID = GetCollectionID(0);
    
  AirhitsCollection = new LyonAirHitsCollection(
  SensitiveDetectorName,
  collectionName[0]);
    
  HCE->AddHitsCollection(
  collectionID,
  AirhitsCollection);
  }


  G4bool LyonAirSensitiveDetector::ProcessHits(G4Step* currentStep,G4TouchableHistory*)
  { 
  //  G4cout<<"particle name: "<<currentStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()<<G4endl;
  LyonAirHit* aHit= new LyonAirHit(currentStep);

  G4cout<<"TrackID: "<<aHit->GetTrackID()
  <<";  AirLayerID: "<<aHit->GetBasicVolumeID()<<G4endl;
  AirhitsCollection->insert(aHit);
  return true;
  }

*/


// As a a first simple appracoh we only look at muons
//((currentStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()== "mu-")
//    || (currentStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()== "mu+"))
//   &&

// As a a first simple appracoh we only look at muons
//((currentStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()== "mu-")
//    || (currentStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()== "mu+"))
//   &&
