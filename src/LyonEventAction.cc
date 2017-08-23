#include "LyonEventAction.hh"
#include "LyonEventActionMessenger.hh"
#include "LyonRunAction.hh"
#include "LyonSensitiveDetector.hh"
#include "LyonTrackHit.hh"
#include "LyonAnalysisToolkit.hh"
#include "LyonPrimaryGeneratorAction.hh"

#include "g4root.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4GenericMessenger.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <iomanip>

LyonEventAction::LyonEventAction() : G4UserEventAction(),
                                     HCID(-1),
                                     TrajDir(3, 0),
                                     TrajPoint(3, 0)
{
  MSN = new LyonEventActionMessenger(this);
  RunDuration = 10.;
  RunTime = 0.;
}

LyonEventAction::~LyonEventAction()
{
  delete MSN;
}

void LyonEventAction::BeginOfEventAction(const G4Event *event)
{
  G4int eventID = event->GetEventID();
  G4cout << "\n---> Begin of event: " << eventID << G4endl;
  if (HCID < 0)
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID("CalorimeterSD/LyonTrackHits");
  hitsCollection = static_cast<LyonTrackHitsCollection *>(event->GetHCofThisEvent()->GetHC(HCID));
  if (!hitsCollection)
  {
    G4String *str = new G4String("Cannot access hitsCollection ");
    G4Exception("LyonEventAction::EndOfEventAction", "1", RunMustBeAborted, *str);
  }

  NbOfPV = event->GetNumberOfPrimaryVertex();
  G4cout << "Number of primary particles: " << NbOfPV << G4endl;

  RunTime = LyonPrimaryGeneratorAction::GetDuration();
  G4cout << "Time: " << RunTime << " seconds" << G4endl;

  if (RunTime > RunDuration)
  {
    G4RunManager *runManager = G4RunManager::GetRunManager();
    runManager->AbortRun();
    G4cout << "\n\n....ooOO00OOoo........ooOO00OOoo........ooOO00OOoo...." << G4endl;
    G4cout << "time out" << G4endl;
    G4cout << "....ooOO00OOoo........ooOO00OOoo........ooOO00OOoo....\n\n"
           << G4endl;
  }
}

void LyonEventAction::EndOfEventAction(const G4Event *event)
{
  LyonAnalysisToolkit *analysis = new LyonAnalysisToolkit(hitsCollection, NbOfPV);
  G4RootAnalysisManager *analysisManager = G4RootAnalysisManager::Instance();

  for (G4int i = 0; i < NbOfPV; i++)
  {

    G4int state = analysis->Analyze(i);

    if (state == 2)
    {
      if (analysis->GetFit(i))
      {
        G4double dtheta = analysis->GetDAngle(i);
        G4ThreeVector Point(analysis->GetPoint(i));

        analysisManager->FillNtupleDColumn(0, 0, Point.getX());
        analysisManager->FillNtupleDColumn(0, 1, Point.getY());
        analysisManager->FillNtupleDColumn(0, 2, Point.getZ());
        analysisManager->FillNtupleDColumn(0, 3, dtheta);
        analysisManager->FillNtupleDColumn(0, 4, RunTime);
        analysisManager->AddNtupleRow(0);

        G4cout << "X = " << Point.getX()
               << "; Y = " << Point.getY()
               << "; Z = " << Point.getZ()
               << "; dTheta = " << dtheta << G4endl;
      }
    }
    else
      G4cout << ">>>>>>>>>Doesn't meet hits requirement<<<<<<<<<<" << G4endl;
    /*
    else
    {
      if (state == 1)
      {
        analysis->GetSingleFit(i);

        G4ThreeVector Vec = analysis->GetSingleDir(i);
        TrajDir[0] = Vec.getX();
        TrajDir[1] = Vec.getY();
        TrajDir[2] = Vec.getZ();

        Vec = analysis->GetPoint(i);
        TrajPoint[0] = Vec.getX();
        TrajPoint[1] = Vec.getY();
        TrajPoint[2] = Vec.getZ();
        analysisManager->FillNtupleDColumn(1, 2, RunTime);
        analysisManager->AddNtupleRow(1);

        G4cout << ">>>>>>>>>Muon absorbed<<<<<<<<<<" << G4endl;
      }
      else
        G4cout << ">>>>>>>>>Doesn't meet hits requirement<<<<<<<<<<" << G4endl;
    }
  }
  */
  }
    hitsCollection = NULL;
    //  AirhitsCollection = NULL;
    delete analysis;
}
/*
        std::vector<G4double> Trks(analysis->GetTracks(i));

        for(G4int j = 0; j<13; j++)
          analysisManager->FillNtupleDColumn(2, j, Trks[j]);
        analysisManager->AddNtupleRow(2);
*/
//    if(state ==0 || state ==1 || state ==2)
//      {

//	analysisManager->FillH1(1, analysis.GetTheta(0,1)); // X = Long axis
//	analysisManager->FillH1(2, analysis.GetTheta(0,2)); // Y = Small axis
//	analysisManager->FillH1(3, analysis.GetTheta(0,0)); // Theta
//	analysisManager->FillH1(4, analysis.GetPhi(0)); // Phi

//      G4cout << "[0] Tot part !" << G4endl;

//	analysisManager->FillNtupleDColumn(0, analysis.GetTheta(0,1));
//	analysisManager->FillNtupleDColumn(1, analysis.GetTheta(0,2));
//	analysisManager->FillNtupleDColumn(2, analysis.GetTheta(0,0));
//	analysisManager->FillNtupleDColumn(3, analysis.GetPhi(0));
//	analysisManager->AddNtupleRow();

//	G4double dphi= analysis.GetPhi(0)-analysis.GetPhi(1);
//	    G4double dthetaX= analysis.GetTheta(0,1)-analysis.GetTheta(1,1);
//	    G4double dthetaY= analysis.GetTheta(0,2)-analysis.GetTheta(1,2);

//	    analysisManager->FillH1(5, dthetaX); // X = Long axis
//	    analysisManager->FillH1(6, dthetaY); // Y = Small axis
//	  analysisManager->FillH1(1, dtheta); // Theta
//	analysisManager->FillH1(2, dphi); // Phi

//	    analysisManager->FillH1(9, fabs(dthetaX)); // X = Long axis
//	    analysisManager->FillH1(10, fabs(dthetaY)); // Y = Small axis
//	    analysisManager->FillH1(11, fabs(dtheta)); // Theta
//	    analysisManager->FillH1(12, fabs(dphi)); // Phi

//	  analysisManager->FillNtupleDColumn(0, 0, Point.getX());
//	  analysisManager->FillNtupleDColumn(0, 1, Point.getY());
//	  analysisManager->FillNtupleDColumn(0, 2, Point.getZ());
//	  analysisManager->FillNtupleDColumn(0, 3, dtheta);
//	  analysisManager->AddNtupleRow(0);

/*
  if(dtheta > 2)
  {
  analysisManager->FillNtupleDColumn(1, 0, Point.getX());
  analysisManager->FillNtupleDColumn(1, 1, Point.getY());
  analysisManager->FillNtupleDColumn(1, 2, Point.getZ());
  analysisManager->FillNtupleDColumn(1, 3, dtheta);
  analysisManager->AddNtupleRow(1);
  }

  if(dtheta > 5)
  {
  analysisManager->FillNtupleDColumn(2, 0, Point.getX());
  analysisManager->FillNtupleDColumn(2, 1, Point.getY());
  analysisManager->FillNtupleDColumn(2, 2, Point.getZ());
  analysisManager->FillNtupleDColumn(2, 3, dtheta);
  analysisManager->AddNtupleRow(2);
  }
*/

/*
  if (AirHCID < 0) 
  AirHCID = G4SDManager::GetSDMpointer()->GetCollectionID("AirSD/LyonAirHits");
  AirhitsCollection = static_cast<LyonAirHitsCollection*>(event->GetHCofThisEvent()->GetHC(AirHCID));
  if ( ! AirhitsCollection )
  {
  G4String* str = new G4String("Cannot access AirhitsCollection ");
  G4Exception("LyonEventAction::EndOfEventAction", "1", RunMustBeAborted, *str);
  }*/

/*	  
	  for(G4int j = 0; j < AirhitsCollection->entries(); j++)
	  {
	  LyonAirHit* AHit = (*AirhitsCollection)[j];
	  if(AHit->GetTrackID() == (i+1))
	  {
	  G4ThreeVector Po(AHit->GetExitPoint());
	  G4int BasicVolumeID = AHit->GetBasicVolumeID();
	  if( BasicVolumeID == 0 || BasicVolumeID == 1 ||  BasicVolumeID == 2)
	  {
	  analysisManager->FillH2(BasicVolumeID,
	  Po.getX()/1000,
	  Po.getY()/1000);
	  G4cout<<"X = "<<Po.getX()<<"  Y = "<<Po.getY()<<G4endl;
	  }
	  else
	  {
	  G4String* str = new G4String("Air layers hit volume id error");
	  G4Exception("LyonEventAction::EndOfEventAction", "1", RunMustBeAborted, *str);
	  }
	  }
	  }
*/
