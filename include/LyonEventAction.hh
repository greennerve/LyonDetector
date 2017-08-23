#ifndef LyonEventAction_h
#define LyonEventAction_h 1

#include "G4UserEventAction.hh"
#include "LyonTrackHit.hh"
#include "globals.hh"

class LyonEventActionMessenger;

class LyonEventAction : public G4UserEventAction
{
public:
  LyonEventAction();
  virtual ~LyonEventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);

  
  void SetRunDuration(G4double val){RunDuration = val;}
  inline G4double GetRunDuration(){return RunDuration;}

  std::vector<double>& GetDir () {return TrajDir;};
  std::vector<double>& GetPoint () {return TrajPoint;};
                     
    
private:
  //  LyonTrackHitsCollection* GetHitsCollection(const G4String& hcName,
  //                                          const G4Event* event) const;

  G4double RunDuration;
  LyonEventActionMessenger* MSN;

  G4int HCID;
  G4int AirHCID;
  G4int NbOfPV;
  LyonTrackHitsCollection* hitsCollection;
  std::vector<G4double> TrajDir;
  std::vector<G4double> TrajPoint;
  G4double RunTime;
  //  LyonAirHitsCollection* AirhitsCollection;
};

#endif

    
