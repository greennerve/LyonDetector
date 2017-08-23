//******************************************************************************
// PrimaryGeneratorAction.hh
//
// This class is a class derived from G4VUserPrimaryGeneratorAction for 
// constructing the process used to generate incident particles.
//
// 1.00 JMV, LLNL, JAN-2007:  First version.
//******************************************************************************
// 
#ifndef LyonPrimaryGeneratorAction_h
#define LyonPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4DataVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "CRYSetup.h"
#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYUtils.h"
#include "vector"
#include "RNGWrapper.hh"
#include "LyonPrimaryGeneratorMessenger.hh"
#include <vector>

class G4Event;

class LyonPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  LyonPrimaryGeneratorAction(const char * filename);
  ~LyonPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);
  void InputCRY();
  void UpdateCRY(std::string* MessInput);
  void CRYFromFile(G4String newValue);
  static inline G4double GetDuration(){return TimeSimulated;};
  G4double GetAverageEnergy(){return energy/genParticles;};
  G4ThreeVector Momen;
  G4ThreeVector Posit;
  inline G4int GetNbOfPrimaryParticles() const {return NbOfPrimaryParticles;};

private:
  std::vector<CRYParticle*> *vect; // vector of generated particles
  G4ParticleTable* particleTable;
  G4ParticleGun* particleGun;
  CRYGenerator* gen;
  G4int InputState;
  LyonPrimaryGeneratorMessenger* gunMessenger;
  G4double energy;
  G4int genParticles;
  G4double CalHalfLength;
  G4double DetectorHalfSizeZ;
  G4double CalHalfDistance;
  G4double WorldHalfSizeZ;
  G4int NbOfPrimaryParticles;
  static G4double TimeSimulated;
};

#endif
