#include "G4Run.hh"
#include "G4RunManager.hh"
//#include "G4VUserPrimaryGeneratorAction.hh"

#include "LyonRunAction.hh"
//#include "LyonRunManager.hh"
#include "LyonEventAction.hh"
#include "LyonDetectorConstruction.hh"
#include "LyonPrimaryGeneratorAction.hh"
#include "g4root.hh"
#include "G4UIcommand.hh"
#include "G4SystemOfUnits.hh"

//LyonRunManager* LyonRunManager::lyonRunMan;

LyonRunAction::LyonRunAction(LyonEventAction *EvA) : EvtA(EvA) {}

LyonRunAction::~LyonRunAction() { ; }

void LyonRunAction::BeginOfRunAction(const G4Run *run)
{
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType()
         << " analysis manager" << G4endl;

  // Open an output file
  G4String fileName = "data";
  analysisManager->OpenFile(fileName);

  analysisManager->CreateNtuple("N1", "X Y Z DTHETA");
  analysisManager->CreateNtupleDColumn("X");
  analysisManager->CreateNtupleDColumn("Y");
  analysisManager->CreateNtupleDColumn("Z");
  analysisManager->CreateNtupleDColumn("DTHETA");
  analysisManager->CreateNtupleDColumn("Time");
  analysisManager->FinishNtuple();
  /*
  analysisManager->CreateNtuple("N2", "Absorbed muon tracks");
  analysisManager->CreateNtupleDColumn("DirVec", EvtA->GetDir());
  analysisManager->CreateNtupleDColumn("Point", EvtA->GetPoint());
  analysisManager->CreateNtupleDColumn("Time");
  analysisManager->FinishNtuple();
  */
  /*
  analysisManager->CreateNtuple("N3", "2 Tracks");
  analysisManager->CreateNtupleDColumn("px1"); //0
  analysisManager->CreateNtupleDColumn("py1"); //1
  analysisManager->CreateNtupleDColumn("pz1"); //2
  analysisManager->CreateNtupleDColumn("mx1"); //3
  analysisManager->CreateNtupleDColumn("my1"); //4
  analysisManager->CreateNtupleDColumn("mz1"); //5
  analysisManager->CreateNtupleDColumn("px2"); //6
  analysisManager->CreateNtupleDColumn("py2"); //7
  analysisManager->CreateNtupleDColumn("pz2"); //8
  analysisManager->CreateNtupleDColumn("mx2"); //9
  analysisManager->CreateNtupleDColumn("my2"); //10
  analysisManager->CreateNtupleDColumn("mz2"); //11
  analysisManager->CreateNtupleDColumn("DTHETA"); //12
  analysisManager->FinishNtuple();
  */
}

void LyonRunAction::EndOfRunAction(const G4Run *aRun)
{
  G4int nofEvents = aRun->GetNumberOfEvent();
  if (nofEvents == 0)
    return;

  G4cout << "Time =" << LyonPrimaryGeneratorAction::GetDuration() << " seconds" << G4endl;

  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

  // complete cleanup
  delete G4RootAnalysisManager::Instance();
}

//  analysisManager->SetHistoDirectoryName("histograms");
//  analysisManager->SetFirstHistoId(0);

// Creating histograms
//  analysisManager->CreateH2("0","Muon flux distribution at Z = 1.4m ", 50, -2.5, 2.5, 50, -2.5, 2.5);
//  analysisManager->CreateH2("1","Muon flux distribution at Z = 0 ", 50, -2.5, 2.5, 50, -2.5, 2.5);
//  analysisManager->CreateH2("2","Muon flux distribution at Z = -1.4m", 50, -2.5, 2.5, 50, -2.5, 2.5);

/* 
  analysisManager->CreateH1("1","Angular distribution (X)", 45, -90.0, 90.0); 
  analysisManager->CreateH1("2","Angular distribution (Y)", 45, -90.0, 90.0);
  analysisManager->CreateH1("3","Angular distribution (Theta)", 45, 0.0, 90.0);
  analysisManager->CreateH1("4","Angular distribution (Phi)", 15, -180.0, 180.0);
  */
/*
     analysisManager->CreateNtuple("N1", "DTHETA_X DTHETA_Y DTHETA Phi");
     analysisManager->CreateNtupleDColumn("DTHETA_X");
     analysisManager->CreateNtupleDColumn("DTHETA_Y");
     analysisManager->CreateNtupleDColumn("DTHETA");
     analysisManager->CreateNtupleDColumn("Phi");
     analysisManager->FinishNtuple();
  */

//     analysisManager->CreateH1("5","Angular deviation (X)", 75, -15.0, 15.0);
//     analysisManager->CreateH1("6","Angular deviation (Y)", 75, -15.0, 15.0);
//      analysisManager->CreateH1("7","Angular deviation (Theta)", 75, -10.0, 10.0);
//      analysisManager->CreateH1("8","Angular deviation (Phi)", 75, -15.0, 15.0);

//     analysisManager->CreateH1("9","Angular deviation (X)", 75, 0.0, 10.0);
//     analysisManager->CreateH1("10","Angular deviation (X)", 75, 0.0, 10.0);
//     analysisManager->CreateH1("11","Angular deviation (X)", 75, 0.0, 10.0);
//     analysisManager->CreateH1("12","Angular deviation (X)", 75, 0.0, 15.0);
/*  
     analysisManager->CreateH1("13","Angularzerzn (X)", 50, -300.0, 300.0); 
     analysisManager->CreateH1("14","Angulza deviation (X)", 50, -300.0, 300.0); 
     analysisManager->CreateH1("15","Anguleazdeviation (X)", 50, -300.0, 300.0);
     //analysisManager->CreateH3("16","Anguleazdeviation (X)", 50, -300.0, 300.0, 50, -300.0, 300.0, 50, -300.0, 300.0);
     */

//  G4double eff = 0.5;

//  LyonRunManager* runManager = LyonRunManager::GetLyonRunMan();
//  G4RunManager* runManager = G4RunManager::GetRunManager();
//  G4cout << "Time =" << (runManager->GetLyonPrimaryGeneratorAction()->GetDuration()/60.0)/eff << " min" << G4endl;
//  G4cout << "Average energy =" << runManager->GetLyonPrimaryGeneratorAction()->GetAverageEnergy() << " Mev" << G4endl;
//  LyonPrimaryGeneratorAction* _LyonPrimaryGeneratorAction =  (LyonPrimaryGeneratorAction*)(runManager->GetUserPrimaryGeneratorAction());
//  LyonPrimaryGeneratorAction* _LyonPrimaryGeneratorAction = (LyonPrimaryGeneratorAction*) runManager->GetUserPrimaryGeneratorAction();

/*
     analysisManager->CreateNtuple("N2", "X Y Z DTHETA");
     analysisManager->CreateNtupleDColumn("X");
     analysisManager->CreateNtupleDColumn("Y");
     analysisManager->CreateNtupleDColumn("Z");
     analysisManager->CreateNtupleDColumn("DTHETA");
     analysisManager->FinishNtuple();

     analysisManager->CreateNtuple("N3", "X Y Z DTHETA");
     analysisManager->CreateNtupleDColumn("X");
     analysisManager->CreateNtupleDColumn("Y");
     analysisManager->CreateNtupleDColumn("Z");
     analysisManager->CreateNtupleDColumn("DTHETA");
     analysisManager->FinishNtuple();
     */

//  G4cout << "Average energy =" << _LyonPrimaryGeneratorAction->GetAverageEnergy() << " Mev" << G4endl;

// print histogram statistics

//  G4cout<<"Run aborted but EndOfRunAction evoqued"<<G4endl;
