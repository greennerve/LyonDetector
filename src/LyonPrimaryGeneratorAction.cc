//******************************************************************************
// PrimaryGeneratorAction.cc
//
// 1.00 JMV, LLNL, Jan-2007:  First version.
//******************************************************************************
//

#define Rad2Deg 57.2957795131

#include "g4root.hh"
#include "G4RunManager.hh"

#include <iomanip>

#include "LyonPrimaryGeneratorAction.hh"
#include "LyonConstants.hh"

#include "G4SystemOfUnits.hh"

#include "G4Event.hh"

//----------------------------------------------------------------------------//

LyonPrimaryGeneratorAction::LyonPrimaryGeneratorAction(const char *inputfile)
{
  // define a particle gun
  particleGun = new G4ParticleGun();

  energy = 0;
  genParticles = 0;

  // Read the cry input file
  std::ifstream inputFile;
  inputFile.open(inputfile, std::ios::in);
  char buffer[1000];

  if (inputFile.fail())
  {
    if (*inputfile != 0) //....only complain if a filename was given
      G4cout << "PrimaryGeneratorAction: Failed to open CRY input file= " << inputfile << G4endl;
    InputState = -1;
  }
  else
  {
    std::string setupString("");
    while (!inputFile.getline(buffer, 1000).eof())
    {
      setupString.append(buffer);
      setupString.append(" ");
    }

    CRYSetup *setup = new CRYSetup(setupString, "../data");

    gen = new CRYGenerator(setup);

    // set random number generator
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState = 0;
  }
  // create a vector to store the CRY particle properties
  vect = new std::vector<CRYParticle *>;

  // Create the table containing all particle names
  particleTable = G4ParticleTable::GetParticleTable();

  // Create the messenger file
  gunMessenger = new LyonPrimaryGeneratorMessenger(this);

  CalHalfLength = RPC_BoxLength / 2;
  DetectorHalfSizeZ = DetectorSizeZ / 2;
  WorldHalfSizeZ = DetectorHalfSizeZ + 5 * cm; // minus 0.5mm to make sure it's in the world
  CalHalfDistance = CalDistance / 2;
}

LyonPrimaryGeneratorAction::~LyonPrimaryGeneratorAction()
{
  delete gunMessenger;
  delete particleGun;
  delete gen;
  delete vect;
  //  G4cout<<"Average muon flux density: "<<genParticles/(gen->timeSimulated()*gen->boxSizeUsed()*gen->boxSizeUsed())<<G4endl;
}

// The folowing functions, is provided by CRY package no info give here
void LyonPrimaryGeneratorAction::InputCRY()
{
  InputState = 1;
}
// The folowing functions, is provided by CRY package no info give here
void LyonPrimaryGeneratorAction::UpdateCRY(std::string *MessInput)
{
  CRYSetup *setup = new CRYSetup(*MessInput, "../data");
  gen = new CRYGenerator(setup);
  // set random number generator
  RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
  setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
  InputState = 0;
}
// The folowing functions, is provided by CRY package no info give here
void LyonPrimaryGeneratorAction::CRYFromFile(G4String newValue)
{
  // Read the cry input file
  std::ifstream inputFile;
  inputFile.open(newValue, std::ios::in);
  char buffer[1000];

  if (inputFile.fail())
  {
    G4cout << "Failed to open input file " << newValue << G4endl;
    G4cout << "Make sure to define the cry library on the command line" << G4endl;
    InputState = -1;
  }
  else
  {
    std::string setupString("");
    while (!inputFile.getline(buffer, 1000).eof())
    {
      setupString.append(buffer);
      setupString.append(" ");
    }

    CRYSetup *setup = new CRYSetup(setupString, "../data");

    gen = new CRYGenerator(setup);

    // set random number generator
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState = 0;
  }
}

G4double LyonPrimaryGeneratorAction::TimeSimulated = 0;

void LyonPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  // This function is called at the very begining of an event

  //  G4cout<<"Time: "<<gen->timeSimulated()<<G4endl;

  if (InputState != 0)
  {
    G4String *str = new G4String("CRY library was not successfully initialized");
    //G4Exception(*str);
    G4Exception("PrimaryGeneratorAction", "1", RunMustBeAborted, *str);
    delete str;
  }

  // Clear current vect of generated particles

  // Generate a new event from CRY package
  bool ReadyToLaunch = false;

  // NB : this loop can taake a WHILE, especially if fired area is small
  //  G4RootAnalysisManager* analysisManager = G4RootAnalysisManager::Instance();
  std::vector<CRYParticle *> GoodMuons;
  GoodMuons.clear();

  do
  {
    // While the particle is not considered as interesting we generate another one
    for (unsigned i = 0; i < vect->size(); i++)
      delete (*vect)[i];
    vect->clear();
    gen->genEvent(vect);
    for (unsigned j = 0; j < vect->size(); j++)
    {
      if ((*vect)[j]->w() >= 0 && (*vect)[j]->ke() <= 10.0)
        continue; // (Mev) Avoid generating particles that can't be detected

      bool IsGoodMuon = false;
      G4ThreeVector pos((*vect)[j]->x() * m, (*vect)[j]->y() * m, (*vect)[j]->z() * m + DetectorHalfSizeZ);
      G4ThreeVector mom((*vect)[j]->u(), (*vect)[j]->v(), (*vect)[j]->w());

      if (abs(pos[0]) < CalHalfLength && abs(pos[1]) < CalHalfLength)
      {
        G4double ft = (-CalHalfDistance - pos[2]) / mom[2];

        IsGoodMuon = fabs(mom[0] * ft + pos[0]) < CalHalfLength;
        IsGoodMuon = IsGoodMuon && (fabs(mom[1] * ft + pos[1]) < CalHalfLength);
        if (IsGoodMuon)
        {
          GoodMuons.push_back((*vect)[j]);
          ReadyToLaunch = true;
          goto ForLoopEnd;
        }
      }

    /*
      if (mom[0] != 0)
      {
        G4double var[2];
        var[0] = (-CalHalfLength - pos[0]) / mom[0];
        var[1] = (CalHalfLength - pos[0]) / mom[0];

        for (G4int k = 0; k < 2; k++)
        {
          IsGoodMuon = fabs(mom[1] * var[k] + pos[1]) < CalHalfLength;
          IsGoodMuon = IsGoodMuon && (fabs(mom[2] * var[k] + pos[2]) < CalHalfDistance);
          if (IsGoodMuon)
          {
            GoodMuons.push_back((*vect)[j]);
            ReadyToLaunch = true;
            goto ForLoopEnd;
          }
        }
      }
      if (mom[1] != 0)
      {
        G4double var[2];
        var[0] = (-CalHalfLength - pos[1]) / mom[1];
        var[1] = (CalHalfLength - pos[1]) / mom[1];

        for (G4int k = 0; k < 2; k++)
        {
          IsGoodMuon = fabs(mom[0] * var[k] + pos[0]) < CalHalfLength;
          IsGoodMuon = IsGoodMuon && (fabs(mom[2] * var[k] + pos[2]) < CalHalfDistance);
          if (IsGoodMuon)
          {
            GoodMuons.push_back((*vect)[j]);
            ReadyToLaunch = true;
            goto ForLoopEnd;
          }
        }
      }
      */
    ForLoopEnd:;
    }
  } while (!ReadyToLaunch);

  NbOfPrimaryParticles = GoodMuons.size();

  for (G4int j = 0; j < NbOfPrimaryParticles; j++)
  {
    particleGun->GeneratePrimaryVertex(anEvent);
    particleGun->SetParticleDefinition(particleTable->FindParticle(GoodMuons[j]->PDGid()));
    particleGun->SetParticleEnergy(GoodMuons[j]->ke() * MeV);
    particleGun->SetParticlePosition(G4ThreeVector(GoodMuons[j]->x() * m, GoodMuons[j]->y() * m, GoodMuons[j]->z() * m + DetectorHalfSizeZ));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(GoodMuons[j]->u(), GoodMuons[j]->v(), GoodMuons[j]->w()));
    particleGun->SetParticleTime(GoodMuons[j]->t() * s);
  }

  for (unsigned i = 0; i < vect->size(); i++)
    delete (*vect)[i];

  vect->clear();

  TimeSimulated = gen->timeSimulated();
}

////////////////////////////////////////////////////////////////////////

/** This function is used during primary generation of particles, it
 * checks if particle is coming on this calorimeter, by looking at
 * the position of the point in the z Plane coresponding to Top or 
 * !Top = Bot of the calorimeter **/
/*
  bool Trigger(G4ThreeVector gunPos, G4ThreeVector gunMom)
  {
 
	 
  bool check(1);
	
  for(unsigned int i(0); i < LyonDetectorConstruction::fGeometryData.Scintillators.size(); ++i)
  {
  G4double x(0);
  G4double y(0);
  G4double z(LyonDetectorConstruction::fGeometryData.Scintillators[i][0].getZ()/CLHEP::m);
  G4double t(0);
		
  t = (z - gunPos.getZ()/CLHEP::m)/gunMom.getZ();
  x = gunPos.getX()/CLHEP::m + t*gunMom.getX();
  y = gunPos.getY()/CLHEP::m + t*gunMom.getY();
	
  check = check && fabs(x-LyonDetectorConstruction::fGeometryData.Scintillators[i][0].getX()/CLHEP::m)<0.5*LyonDetectorConstruction::fGeometryData.Scintillators[i][1].getX()/CLHEP::m;
  check = check && fabs(y-LyonDetectorConstruction::fGeometryData.Scintillators[i][0].getY()/CLHEP::m)<0.5*LyonDetectorConstruction::fGeometryData.Scintillators[i][1].getY()/CLHEP::m;
  }
	
  return check;
  }
*/

//	  G4cout<<"Initial muon Z coordinates: "<<(*vect)[j]->z()<<G4endl;
//	  bool checkIn = LyonDetectorConstruction::fCalorimeters[0].IsParticleComing(pos,mom,1); // Fires first layer
//	  bool checkOut = LyonDetectorConstruction::fCalorimeters[1].IsParticleComing(pos,mom,0); // Fires last layer
//	  bool checkSci = Trigger(pos,mom); //Fires scintillators

//	      check = (*vect)[j]->ke()>10.0; // (Mev) Avoid do generate particles that we can't detect
//	  check = checkIn && checkOut && checkEnergy; // Global check
//	  check =  check && (fabs(pos[0]) < CalHalfLength);
//	  check =  check && (fabs(pos[1]) < CalHalfLength);
//	  check =  check && (fabs(mom[0]/fabs(mom[2])*DetectorSizeZ + pos[0]) < CalHalfLength);
//	  check =  check && (fabs(mom[1]/fabs(mom[2])*DetectorSizeZ + pos[1]) < CalHalfLength);

//	  ArrivalCheck[1] = (fabs(mom[0]/fabs(mom[2])*DetectorSizeZ + pos[0]) < CalHalfLength);

//	      energy += (*vect)[j]->ke(); // For energy average
//	      genParticles ++; // For energy average

//	  G4cout<<"Theta = "<<180.0-mom.getTheta()*Rad2Deg<<G4endl;

//	  analysisManager->FillH1(1, -atan(mom.x()/mom.z())*Rad2Deg); // X = Long axis
// analysisManager->FillH1(2, -atan(mom.y()/mom.z())*Rad2Deg); // Y = Small axis
// analysisManager->FillH1(3, 180-mom.getTheta()*Rad2Deg); // Theta
// analysisManager->FillH1(4, mom.getPhi()*Rad2Deg); // Phi

//  Momen = G4ThreeVector((*vect)[0]->u(), (*vect)[0]->v(), (*vect)[0]->w());
//  Posit = G4ThreeVector((*vect)[0]->x()*m, (*vect)[0]->y()*m, (*vect)[0]->z()*m+DetectorSizeZ/2);
/*
  if(vect->size() > 1)
  {
  G4String* str = new G4String("More than one muon generated");
  G4Exception("LyonPrimaryGenerationAction", "1", RunMustBeAborted, *str);
  }
*/
// Once we found a good one (or a set of one), we simply fire it/them
