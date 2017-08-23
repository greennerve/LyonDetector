#ifndef LyonTrackHit_h
#define LyonTrackHit_h

#include "G4VHit.hh"
#include "G4ThreeVector.hh"

class G4Step;

class LyonTrackHit : public G4VHit
{
public:
  LyonTrackHit();
  LyonTrackHit(const G4Step* aStep, G4ThreeVector& precisionRandomVector);
  virtual ~LyonTrackHit() {;};
  // virtual void Draw();
  // virtual void Print();
public:

  void SetExitPoint(G4ThreeVector xyz){fExitPoint = xyz;};
  void SetChamber(G4int fc) { fChamber = fc;};
  void SetCalo(G4int fca) { fCalo = fca;};
  void SetBasicVolumeID(G4int fB) { fBasicVolume = fB;};
  void SetTrackID(G4int TI) { TrackID = TI;};
  
  const G4ThreeVector& GetExitPoint() const {return fExitPoint;}
  inline G4int GetChamber() const {return fChamber;}
  inline G4int GetCalo() const {return fCalo;}
  inline G4int GetBasicVolumeID() const {return fBasicVolume;}
  inline G4int GetTrackID () const {return TrackID;}

  
private:
  //  G4double fEnergyDeposited;
  G4ThreeVector fExitPoint;
  G4int fChamber;
  G4int fCalo;
  G4int fBasicVolume;
  //  G4int fpdgID;
  //  G4double fTime;
  G4int TrackID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/*
class LyonAirHit : public G4VHit
{
public:
  LyonAirHit(const G4Step* aStep);
  virtual ~LyonAirHit() {;}

private:
  G4ThreeVector fExitPoint;
  G4int fBasicVolume;
  G4int TrackID;

public:
  const G4ThreeVector& GetExitPoint() const {return fExitPoint;}
  inline G4int GetBasicVolumeID() const {return fBasicVolume;}
  inline G4int GetTrackID () const {return TrackID;}
};

*/

#include "G4THitsCollection.hh"
typedef G4THitsCollection<LyonTrackHit> LyonTrackHitsCollection;
//typedef G4THitsCollection<LyonAirHit> LyonAirHitsCollection;

#endif
