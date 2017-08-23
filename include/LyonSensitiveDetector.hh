#ifndef LyonSensitiveDetector_h
#define LyonSensitiveDetector_h

#include "G4VSensitiveDetector.hh"
#include "LyonTrackHit.hh"
#include "Randomize.hh"

class G4Step;
class G4HCofThisEvent;
class G4Event;

class LyonSensitiveDetector : public G4VSensitiveDetector
{
public:
  LyonSensitiveDetector(G4String name);
  virtual ~LyonSensitiveDetector();
  virtual void Initialize(G4HCofThisEvent*);
  virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent*){;}

private:
  LyonTrackHitsCollection * hitsCollection;
  G4int collectionID;
  G4double efficiency;
  G4double precision;
  G4double stdDev;
  G4int NbOfTracks;
  std::vector<G4int> CalID;
  std::vector<G4int> RPCID;
  G4Event* event;
  CLHEP::RandGauss* fRandomGauss;
  
};
/*
class LyonAirSensitiveDetector : public G4VSensitiveDetector
{
public:
  LyonAirSensitiveDetector(G4String name);
  virtual ~LyonAirSensitiveDetector() {;}
  virtual void Initialize(G4HCofThisEvent*);
  virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent*){;}

private:
  LyonAirHitsCollection * AirhitsCollection;
  G4int collectionID;
};
*/
#endif
