#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UIExecutive.hh"

#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
//#include "G4ScoringManager.hh"

#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"


#include "LyonEventAction.hh"
#include "LyonDetectorConstruction.hh"
#include "LyonPrimaryGeneratorAction.hh"
#include "LyonRunAction.hh"

#include "Randomize.hh" //modif pour random seed



int main(int argc,char** argv)
{
  //modif pour random seed
  G4int seed=-1;
  if(argc > 2) seed=atoi(argv[2]);
  //if(seed<0)   seed=time(0);
  if(seed<0)   seed=1;
  G4cout << "Chosen seed is=" << seed << G4endl;
  CLHEP::HepRandom::setTheSeed(seed);

  // Construct the default run manager
  //  LyonRunManager* runManager = new LyonRunManager;
  G4RunManager* runManager = new G4RunManager;
  //Activate command-based scorer 
  //should be called right afetr G4RunManager instantition
//  G4ScoringManager::GetScoringManager();

  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  //G4VUserPhysicsList* physics;
  auto physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // set mandatory initialization classes
  G4VUserDetectorConstruction* detector = new LyonDetectorConstruction;
  runManager->SetUserInitialization(detector);

 
  
  // set mandatory user action class
  LyonPrimaryGeneratorAction* priAct = new LyonPrimaryGeneratorAction("setup.file"); // Sale histoire pour récupérer le nombre de particules générés mais non gardées
  runManager->SetUserAction(priAct);

  LyonEventAction* EvA = new LyonEventAction;
  
  runManager->SetUserAction(new LyonRunAction(EvA));

  runManager->SetUserAction(EvA);
 
  // Initialize G4 kernel
  runManager->Initialize();

  // Get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
 
  
  if(argc==1)  // Define (G)UI terminal for interactive mode
    {
      G4UIExecutive* UIExe=new G4UIExecutive(argc,argv);
      UI->ApplyCommand("/control/execute mac/vis.mac");
      //session->SessionStart();
      UIExe->SessionStart();
      //delete session;
      delete UIExe;
    }
  else   // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  
  delete visManager;
  delete runManager;

  return 0;
}
