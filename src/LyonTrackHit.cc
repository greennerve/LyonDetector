#include "LyonTrackHit.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4VProcess.hh"
#include "G4UIcommand.hh"
#include "g4root.hh"
#include "LyonDetectorConstruction.hh"

LyonTrackHit::LyonTrackHit():G4VHit(),fExitPoint(0,0,0),
   fChamber(0), fCalo(0),fBasicVolume(0), TrackID(0){}


LyonTrackHit::LyonTrackHit(const G4Step* aStep, G4ThreeVector& precisionRandomVector)
{

  fExitPoint=aStep->GetPostStepPoint()->GetPosition()+ precisionRandomVector;  
  G4TouchableHistory * theTouchable = 
    (G4TouchableHistory *) aStep->GetPreStepPoint()->GetTouchable();
  fBasicVolume = theTouchable->GetReplicaNumber(0);
  fChamber = theTouchable->GetReplicaNumber(1);
  fCalo = theTouchable->GetReplicaNumber(2);
  TrackID = aStep->GetTrack()->GetTrackID();
}
/*
LyonAirHit::LyonAirHit(const G4Step* aStep)
{
  fExitPoint=aStep->GetPostStepPoint()->GetPosition(); 
  G4TouchableHistory * theTouchable = 
    (G4TouchableHistory *) aStep->GetPreStepPoint()->GetTouchable();
  fBasicVolume = theTouchable->GetReplicaNumber(0);
  TrackID = aStep->GetTrack()->GetTrackID();
}
*/


  //  fEnergyDeposited=aStep->GetTotalEnergyDeposit();  
  // following used to compute x/y distribution
  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // analysisManager->FillH1(5, fExitPoint.getX()/CLHEP::cm);
  // analysisManager->FillH1(6, fExitPoint.getY()/CLHEP::cm);
    
  //  fpdgID=aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();
  //  fTime=aStep->GetPostStepPoint()->GetGlobalTime();
