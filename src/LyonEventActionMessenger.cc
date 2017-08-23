//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: ExN03EventActionMessenger.cc,v 1.12 2006/10/26 14:28:00 allison Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "LyonEventActionMessenger.hh"
#include "LyonEventAction.hh"
//#include "G4UIdirectory.hh"
//#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LyonEventActionMessenger::LyonEventActionMessenger(LyonEventAction* eve)
  :EventAct(eve)
{
  //  LyonDetectordetDir = new G4UIdirectory("/field/");
  //  LyonDetectordetDir->SetGuidance("LyonDetector field tracking control.");
  EventCmd = new G4UIcmdWithADoubleAndUnit("/run/SetRunDuration",this);  
  EventCmd->SetGuidance("Set the duration of a run");
  EventCmd->SetParameterName("RunDuration",false,false);
  EventCmd->SetDefaultValue(10.);
  EventCmd->SetUnitCategory("Time");
  EventCmd->SetDefaultUnit("s");
  EventCmd->SetRange("RunDuration > 0");
}

LyonEventActionMessenger::~LyonEventActionMessenger()
{
  delete EventCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LyonEventActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == EventCmd)
    {
      EventAct->SetRunDuration(EventCmd->GetNewDoubleValue(newValue)/second);
      G4cout<<"Run duration is set to "<<newValue<<G4endl;
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
