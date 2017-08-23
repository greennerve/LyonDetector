#ifndef LyonRunAction_h
#define LyonRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class LyonEventAction;

class LyonRunAction : public G4UserRunAction
{

public:
  LyonRunAction(LyonEventAction* EvA);
  ~LyonRunAction();
  
public:
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

private:
  LyonEventAction* EvtA;
};


#endif
