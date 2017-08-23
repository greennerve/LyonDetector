#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4UIcommand.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4SDParticleFilter.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "LyonDetectorConstruction.hh"
#include "LyonSensitiveDetector.hh"
#include "LyonConstants.hh"

struct VolumeInfo
{
    G4String VolumeName;
    G4String Material;
    G4double Thickness;

    VolumeInfo(G4String name, G4String Mate, G4double Thic) : VolumeName(name),
                                                              Material(Mate),
                                                              Thickness(Thic){};
};

LyonDetectorConstruction::LyonDetectorConstruction() : fVolumeInfo(14, nullptr),
                                                       CalLV(3, nullptr),
                                                       LogicalVolumes(_NbOfVolumes, nullptr),
                                                       RotCollector(4, nullptr),
                                                       VisCollector(2, nullptr)
{

    MuonFilter = new G4SDParticleFilter("MuonFilter");
    MuonFilter->add("mu+");
    MuonFilter->add("mu-");

    SDetector = new LyonSensitiveDetector("CalorimeterSD");
    SDetector->SetFilter(MuonFilter);
    G4SDManager::GetSDMpointer()->AddNewDetector(SDetector);

    VisCollector[0] = new G4VisAttributes(G4Colour::Red());
    VisCollector[1] = new G4VisAttributes(G4Colour::Blue());

    for (auto visAt : VisCollector)
    {
        visAt->SetForceSolid(true);
    };
}

LyonDetectorConstruction::~LyonDetectorConstruction()
{
    for (auto VisAt : VisCollector)
    {
        delete VisAt;
    }

    for (auto Rot : RotCollector)
    {
        delete Rot;
    }
    for (auto vInfo : fVolumeInfo)
    {
        delete vInfo;
    }
    delete MuonFilter;
}

void LyonDetectorConstruction::DefineMaterials()
{

    /* ------------------------------------------------------------------ */
    /*                              MATERIALS                             */
    /* ------------------------------------------------------------------ */

    //Get standard stuff from NIST database :
    G4NistManager *nistManager = G4NistManager::Instance();

    //G4Material *Air=nistManager->FindOrBuildMaterial("G4_AIR");
    G4Material *Fer = nistManager->FindOrBuildMaterial("G4_Fe");
    G4Material *Cr = nistManager->FindOrBuildMaterial("G4_Cr");
    G4Material *Ni = nistManager->FindOrBuildMaterial("G4_Ni");
    //  G4Material *Al=nistManager->FindOrBuildMaterial("G4_Al");

    G4double steeldensity = 7.87 * g / cm3;
    G4double fractionMassFe = 0.70611;
    G4double fractionMassCr = 0.18784;
    G4double fractionMassNi = 0.10605;

    G4Material *Steel = new G4Material("Steel", steeldensity, 3);
    Steel->AddMaterial(Fer, fractionMassFe);
    Steel->AddMaterial(Cr, fractionMassCr);
    Steel->AddMaterial(Ni, fractionMassNi);

    //Definition du g10
    double a = 1.01 * g / mole;
    G4Element *elH = new G4Element("Hydrogen", "H", 1, a);

    a = 12.01 * g / mole;
    G4Element *elC = new G4Element("Carbon", "C", 6, a);

    a = 16.00 * g / mole;
    G4Element *elO = new G4Element("Oxygen", "O", 8, a);

    a = 32.06 * g / mole;
    G4Element *elS = new G4Element("Sulfur", "S", 16, a);

    a = 19.00 * g / mole;
    G4Element *elF = new G4Element("Fluor", "F", 9, a);

    G4double epoxydensity = 1.3 * g / cm3;
    G4Material *epoxy = new G4Material("epoxy", epoxydensity, 3);
    epoxy->AddElement(elC, 15);
    epoxy->AddElement(elO, 7);
    epoxy->AddElement(elH, 44);

    G4Material *Si02 = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4Material *Cl = nistManager->FindOrBuildMaterial("G4_Cl");

    G4double g10density = 1.7 * g / cm3;
    G4Material *g10 = new G4Material("g10", g10density, 3);
    g10->AddMaterial(epoxy, 14.7 * perCent);
    g10->AddMaterial(Cl, 8 * perCent);
    g10->AddMaterial(Si02, 77.3 * perCent);

    /* = Gaz definition  ================================================ */

    //SF6
    G4double SF6density = 6.27e-3 * g / cm3;
    G4Material *SF6 = new G4Material("SF6", SF6density, 2);
    SF6->AddElement(elS, 1);
    SF6->AddElement(elF, 6);

    //Isobutane
    G4double Isobutanedensity = 2.51e-3 * g / cm3;
    G4Material *Isobutane = new G4Material("Isobutane", Isobutanedensity, 2);
    Isobutane->AddElement(elC, 4);
    Isobutane->AddElement(elH, 10);

    //CO2
    G4double CO2density = 1.799e-3 * g / cm3;
    G4Material *CO2 = new G4Material("CO2", CO2density, 2);
    CO2->AddElement(elC, 1);
    CO2->AddElement(elO, 2);

    //TFE
    G4double TFEdensity = 4.25e-3 * g / cm3;
    G4Material *TFE = new G4Material("TFE", TFEdensity, 3);
    TFE->AddElement(elC, 2);
    TFE->AddElement(elH, 2);
    TFE->AddElement(elF, 4);

    //RPCGaz
    G4double RPCGazdensity = 4.13e-3 * g / cm3;
    G4Material *RPCGaz = new G4Material("RPCGaz", RPCGazdensity, 3);
    RPCGaz->AddMaterial(TFE, 93 * perCent);
    RPCGaz->AddMaterial(Isobutane, 5 * perCent);
    RPCGaz->AddMaterial(SF6, 2 * perCent);
}

void LyonDetectorConstruction::DefineVolumeInfo()
{
    //Differents volume thickness
    fVolumeInfo[_RPC_Absorber] = new VolumeInfo("RPC_Absorber",
                                                "G4_Al",
                                                RPC_Absorber_Thickness); //0,Absorber Layer

    fVolumeInfo[_RPC_Electronics] = new VolumeInfo("RPC_Electronics",
                                                   "g10",
                                                   RPC_ChipPackage_Thickness); //1, electronics (chip + epoxy)

    fVolumeInfo[_RPC_PCB] = new VolumeInfo("RPC_PCB",
                                           "g10",
                                           RPC_PCB_Thickness); // 2,PCB Layer

    fVolumeInfo[_RPC_mylarAnode] = new VolumeInfo("RPC_mylarAnode",
                                                  "G4_MYLAR",
                                                  RPC_mylarAnode_Thickness); //3,mylar (anode)

    fVolumeInfo[_RPC_GraphiteAnode] = new VolumeInfo("RPC_GraphiteAnode",
                                                     "G4_GRAPHITE",
                                                     RPC_GraphiteAnode_Thickness); //4,graphite

    fVolumeInfo[_RPC_ThinGlass] = new VolumeInfo("RPC_ThinGlass",
                                                 "G4_Pyrex_Glass",
                                                 RPC_ThinGlass_Thickness); //5,glass (thin)

    fVolumeInfo[_RPC_GasChamber] = new VolumeInfo("RPC_GasChamber",
                                                  "RPCGaz",
                                                  RPC_GasChamber_Thickness); //6,Gas gap

    fVolumeInfo[_RPC_ThickGlass] = new VolumeInfo("RPC_ThickGlass",
                                                  "G4_Pyrex_Glass",
                                                  RPC_ThickGlass_Thickness); //7,glass (thick)

    fVolumeInfo[_RPC_GraphiteCathode] = new VolumeInfo("RPC_GraphiteCathode",
                                                       "G4_GRAPHITE",
                                                       RPC_GraphiteCathode_Thickness); //8,graphite (cathode)

    fVolumeInfo[_RPC_mylarCathode] = new VolumeInfo("RPC_mylarCathode",
                                                    "G4_MYLAR",
                                                    RPC_mylarCathode_Thickness); //9,mylar (cathode)

    fVolumeInfo[_RPC_Free] = new VolumeInfo("RPC_Free",
                                            "G4_AIR",
                                            RPC_Free_Thickness); //10, left space

    fVolumeInfo[_RPC_AdditionalAlu] = new VolumeInfo("RPC_AdditionalAlu",
                                                     "G4_Al",
                                                     RPC_AdditionalAlu_Thickness); //11, additional Alu

    fVolumeInfo[_RPC] = new VolumeInfo("RPC",
                                       "G4_AIR",
                                       RPC_Thickness); //12,RPC
    fVolumeInfo[_Calorimeter] = new VolumeInfo("Calorimeter",
                                               "G4_AIR",
                                               Cal_Thickness); //13,
}

G4LogicalVolume *LyonDetectorConstruction::ConstructCalorimeter(G4double XHalfLength, G4double YHalfLength)
{
    G4NistManager *nistManager = G4NistManager::Instance();

    G4Box *SolidBox = new G4Box("PhysicalRPC", XHalfLength, YHalfLength, fVolumeInfo[_RPC]->Thickness / 2);
    LogicalVolumes[_RPC] = new G4LogicalVolume(SolidBox,
                                               nistManager->FindOrBuildMaterial(fVolumeInfo[_RPC]->Material),
                                               "Logical" + fVolumeInfo[_RPC]->VolumeName);

    G4double Position_Z = -fVolumeInfo[_RPC]->Thickness / 2; //Z axis coordinates

    bool fCheckOverlaps = 1;
    G4int NbOfLayers = _NbOfVolumes - 2;

    for (G4int i = 0; i < NbOfLayers; i++)
    {
        SolidBox = new G4Box("Physical" + fVolumeInfo[i]->VolumeName, XHalfLength, YHalfLength, fVolumeInfo[i]->Thickness / 2);
        LogicalVolumes[i] = new G4LogicalVolume(SolidBox,
                                                nistManager->FindOrBuildMaterial(fVolumeInfo[i]->Material),
                                                "Logical" + fVolumeInfo[i]->VolumeName);
        Position_Z += fVolumeInfo[i]->Thickness / 2;
        new G4PVPlacement(0,                               // no rotation
                          G4ThreeVector(0, 0, Position_Z), // at (x,y,z)
                          LogicalVolumes[i],               // its logical volume
                          fVolumeInfo[i]->VolumeName,      // its name
                          LogicalVolumes[_RPC],            // its mother  volume
                          false,                           // no boolean operations
                          0,                               // copy number
                          fCheckOverlaps);                 // checking overlaps
        Position_Z += fVolumeInfo[i]->Thickness / 2;
    }

    //==================== = Calorimeter itselfs =================================

    SolidBox = new G4Box("PhysicalCalorimeter", XHalfLength, YHalfLength, fVolumeInfo[_Calorimeter]->Thickness / 2);

    LogicalVolumes[_Calorimeter] = new G4LogicalVolume(SolidBox,
                                                       nistManager->FindOrBuildMaterial(fVolumeInfo[_Calorimeter]->Material),
                                                       "Logical" + fVolumeInfo[_Calorimeter]->VolumeName);

    Position_Z = (-Cal_Thickness + RPC_Thickness) / 2;
    G4double deltaZ = RPC_Thickness + RPC_Distance;
    for (G4int i = 0; i < NbOfRPC; i++) // Layers placement
    {
        new G4PVPlacement(
            0, // 180° rotation due to old orientation of this code
            G4ThreeVector(0, 0, Position_Z),
            LogicalVolumes[_RPC],                    // Logical volume
            "RPC" + G4UIcommand::ConvertToString(i), // Name
            LogicalVolumes[_Calorimeter],            // Mother's volume
            false,
            i,               // Copy number
            fCheckOverlaps); // checking overlaps
        Position_Z += deltaZ;
    }

    //  Cal_Thickness +=6*mm; //Avoid overlaps
    LogicalVolumes[_RPC_GasChamber]->SetSensitiveDetector(SDetector); // In order to make possible analyze

    /* = Visualization attributes =================================== */

    LogicalVolumes[_Calorimeter]->SetVisAttributes(G4VisAttributes::Invisible);
    LogicalVolumes[_RPC]->SetVisAttributes(G4VisAttributes::Invisible);

    LogicalVolumes[_RPC_AdditionalAlu]->SetVisAttributes(VisCollector[0]); //for debug

    LogicalVolumes[_RPC_Absorber]->SetVisAttributes(VisCollector[1]);

    return LogicalVolumes[_Calorimeter];
}

G4VPhysicalVolume *LyonDetectorConstruction::Construct()
{
    DefineMaterials();
    DefineVolumeInfo();

    CalLV[0] = ConstructCalorimeter(RPC_BoxLength / 2 + RPC_Distance, RPC_BoxLength / 2 + RPC_Distance);
    //  CalLV[1] = ConstructCalorimeter(RPC_BoxLength/2, RPC_BoxLength/2);
    //CalLV[2] = ConstructCalorimeter(RPC_BoxLength / 2, CalDistance / 2);
    //  CalLV[3] = ConstructCalorimeter(RPC_BoxLength/2+RPC_Distance*6, CalDistance/2);

    G4NistManager *man = G4NistManager::Instance();
    G4Material *Air = man->FindOrBuildMaterial("G4_AIR");
    G4Material *Pb = man->FindOrBuildMaterial("G4_Pb");
    G4Material *Steel = man->FindOrBuildMaterial("Steel");
    G4Material *Al = man->FindOrBuildMaterial("G4_Al");

    //G4Material *Water = man->FindOrBuildMaterial("G4_WATER");

    G4Isotope *isoU235 = new G4Isotope("U235", 92, 235, 235.0439242 * g / mole);
    G4Isotope *isoU238 = new G4Isotope("U238", 92, 238, 238.0507847 * g / mole);
    G4Element *elenrichedU = new G4Element("enriched U", "U", 2);
    elenrichedU->AddIsotope(isoU235, 15. * perCent);
    elenrichedU->AddIsotope(isoU238, 85. * perCent);
    G4Material *matenrichedU = new G4Material("U for nuclear power generation",
                                              19.050 * g / cm3, 1, kStateSolid);
    matenrichedU->AddElement(elenrichedU, 1.);

    /* = World construction ============================================= */

    G4Box *solidWorld = new G4Box(
        "World",
        WorldSizeX / 2,
        WorldSizeY / 2,
        WorldSizeZ / 2);

    G4LogicalVolume *logicWorld = new G4LogicalVolume(
        solidWorld, //its solid
        Air,        //its material
        "World");   //its name

    G4VPhysicalVolume *physiWorld = new G4PVPlacement(
        0,               //no rotation
        G4ThreeVector(), //at (0,0,0)
        logicWorld,      //its logical volume
        "World",         //its name
        0,               //its mother  volume
        false,           //no boolean operation
        0,               //copy number
        true);           //check overlaps

    new G4PVPlacement(
        0,                                                      //no rotation
        G4ThreeVector(0, 0, (Cal_Thickness + CalDistance) / 2), //
        CalLV[0],                                               //its logical volume
        "Cal0",                                                 //its name
        logicWorld,                                             //its mother  volume
        false,                                                  //no boolean operation
        0,                                                      //copy number
        true);                                                  //check overlaps

    new G4PVPlacement(
        0,                                                       //no rotation
        G4ThreeVector(0, 0, -(Cal_Thickness + CalDistance) / 2), //
        CalLV[0],                                                //its logical volume
        "Cal1",                                                  //its name
        logicWorld,                                              //its mother  volume
        false,                                                   //no boolean operation
        1,                                                       //copy number
        true);                                                   //check overlaps
                                                                 /*
    G4Box *IronBoxOut = new G4Box(
        "IronBoxOut", //its name
        2 * m,        //size
        1. * m,
        1. * m);

    G4LogicalVolume *IronBoxOutLV = new G4LogicalVolume(
        IronBoxOut,    //its solid
        Steel,         //its material
        "IronBoxOut"); //its name
    IronBoxOutLV->SetVisAttributes(VisCollector[1]);
    new G4PVPlacement(
        0, //no rotation
        G4ThreeVector(),
        IronBoxOutLV, //its logical volume
        "IronBoxOut", //its name
        logicWorld,   //its mother  volume
        false,        //no boolean operation
        1,            //copy number
        true);
        
    
  G4Box *IronBoxIn = new G4Box(
      "IronBoxOut", //its name
      1.9 * m,      //size
      0.9 * m,
      0.9 * m);

  G4LogicalVolume *IronBoxInLV = new G4LogicalVolume(
      IronBoxIn,    //its solid
      Air,          //its material
      "IronBoxIn"); //its name
 
  IronBoxInLV->SetVisAttributes(VisCollector[1]);
  new G4PVPlacement(
      0, //no rotation
      G4ThreeVector(),
      IronBoxInLV,  //its logical volume
      "IronBoxIn",  //its name
      IronBoxOutLV, //its mother  volume
      false,        //no boolean operation
      1,            //copy number
      true);

                                                                 
    G4double TarL = 15 * cm;

    G4Box *solidTarget2 = new G4Box(
        "Target2", //its name
        TarL,      //size
        TarL,
        TarL);

    G4LogicalVolume *logicTarget2 = new G4LogicalVolume(
        solidTarget2, //its solid
        matenrichedU, //its material
        "Target2");   //its name
    logicTarget2->SetVisAttributes(VisCollector[0]);
    new G4PVPlacement(
        0, //no rotation
        G4ThreeVector(1.5 * m, 0.75 * m, 0.75 * m),
        logicTarget2, //its logical volume
        "U1",         //its name
        IronBoxOutLV, //its mother  volume
        false,        //no boolean operation
        0,            //copy number
        true);

    //G4double WaterL = 1 * m;
    
G4Tubs* EyeSoli = new G4Tubs(
                          "EyeSoli",
                          0,
                          30*cm,
                          15*cm,
                          0,
                          CLHEP::twopi
                        );

  G4LogicalVolume *EyeLV = new G4LogicalVolume(
      EyeSoli, //its solid
      matenrichedU, //its material
      "EyeLV");   //its name

  new G4PVPlacement(
      0, //no rotation
      G4ThreeVector(1.3 * m, 1.3 * m, 0),
      EyeLV, //its logical volume
      "Eye1",         //its name
      logicWorld,   //its mother  volume
      false,        //no boolean operation
      1,            //copy number
      true);

  new G4PVPlacement(
      0, //no rotation
      G4ThreeVector(-1.3 * m, 1.3 * m, 0),
      EyeLV, //its logical volume
      "Eye2",         //its name
      logicWorld,   //its mother  volume
      false,        //no boolean operation
      2,            //copy number
      true);

  G4Tubs *MouthSoli = new G4Tubs(
      "MouthSoli",
      1.2*m,
      1.4 * m,
      15 * cm,
      CLHEP::pi,
      CLHEP::pi);

  G4LogicalVolume *MouthLV = new G4LogicalVolume(
      MouthSoli, //its solid
      matenrichedU, //its material
      "MouthLV");   //its name

  new G4PVPlacement(
      0, //no rotation
      G4ThreeVector(0, 0, 0),
      MouthLV, //its logical volume
      "Mouth",         //its name
      logicWorld,   //its mother  volume
      false,        //no boolean operation
      1,            //copy number
      true);
      */
    /*  
  G4double TarL = 15 * cm;
  G4Box *solidTarget2 = new G4Box(
      "Target2", //its name
      TarL,      //size
      TarL,
      TarL);
  G4LogicalVolume *logicTarget2 = new G4LogicalVolume(
      solidTarget2, //its solid
      matenrichedU, //its material
      "Target2");   //its name
  logicTarget2->SetVisAttributes(VisCollector[0]);
  new G4PVPlacement(
      0, //no rotation
      G4ThreeVector(-2.2 * m, 2.2 * m, 1.2 * m),
      logicTarget2, //its logical volume
      "U1",         //its name
      logicWorld,   //its mother  volume
      false,        //no boolean operation
      1,            //copy number
      true);

  new G4PVPlacement(
      0, //no rotation
      G4ThreeVector(2.2 * m, 2.2 * m, 0.4 * m),
      logicTarget2, //its logical volume
      "U2",         //its name
      logicWorld,   //its mother  volume
      false,        //no boolean operation
      2,            //copy number
      true);

  new G4PVPlacement(
      0, //no rotation
      G4ThreeVector(2.2 * m, -2.2 * m, -0.4 * m),
      logicTarget2, //its logical volume
      "U3",         //its name
      logicWorld,   //its mother  volume
      false,        //no boolean operation
      3,            //copy number
      true);

  new G4PVPlacement(
      0, //no rotation
      G4ThreeVector(-2.2 * m, -2.2 * m, -1.2 * m),
      logicTarget2, //its logical volume
      "U4",         //its name
      logicWorld,   //its mother  volume
      false,        //no boolean operation
      4,            //copy number
      true);
*/

    G4double TarL = 15 * cm;

    G4Box *solidTarget = new G4Box(
        "Target", //its name
        TarL,     //size
        TarL,
        TarL);
    G4LogicalVolume *logicTarget = new G4LogicalVolume(
        solidTarget, //its solid
        Pb,          //its material
        "Target");   //its name
    logicTarget->SetVisAttributes(VisCollector[0]);
    new G4PVPlacement(
        0, //no rotation
        G4ThreeVector(1.5 * m, 0., 0.),
        logicTarget, //its logical volume
        "Target",    //its name
        logicWorld,  //its mother  volume
        false,       //no boolean operation
        0,           //copy number
        true);

    G4Box *solidTarget2 = new G4Box(
        "Target2", //its name
        TarL,      //size
        TarL,
        TarL);
    G4LogicalVolume *logicTarget2 = new G4LogicalVolume(
        solidTarget2, //its solid
        matenrichedU, //its material
        "Target2");   //its name
    logicTarget2->SetVisAttributes(VisCollector[0]);
    new G4PVPlacement(
        0, //no rotation
        G4ThreeVector(-1.5 * m, 0., 0.),
        logicTarget2, //its logical volume
        "Target2",    //its name
        logicWorld,   //its mother  volume
        false,        //no boolean operation
        0,            //copy number
        true);

    G4Box *solidTarget3 = new G4Box(
        "Target3", //its name
        TarL,      //size
        TarL,
        TarL);
    G4LogicalVolume *logicTarget3 = new G4LogicalVolume(
        solidTarget3, //its solid
        Al,           //its material
        "Target3");   //its name
    logicTarget3->SetVisAttributes(VisCollector[0]);
    new G4PVPlacement(
        0, //no rotation
        G4ThreeVector(0., 1.5 * m, 0.),
        logicTarget3, //its logical volume
        "Target3",    //its name
        logicWorld,   //its mother  volume
        false,        //no boolean operation
        0,            //copy number
        true);

    G4Box *solidTarget4 = new G4Box(
        "Target4", //its name
        TarL,      //size
        TarL,
        TarL);
    G4LogicalVolume *logicTarget4 = new G4LogicalVolume(
        solidTarget4, //its solid
        Steel,        //its material
        "Target4");   //its name
    logicTarget4->SetVisAttributes(VisCollector[0]);
    new G4PVPlacement(
        0, //no rotation
        G4ThreeVector(0., -1.5 * m, 0.),
        logicTarget4, //its logical volume
        "Target4",    //its name
        logicWorld,   //its mother  volume
        false,        //no boolean operation
        0,            //copy number
        true);

    //....ooOO00OOoo........ooOO00OOoo........ooOO00OOoo........ooOO00OOoo........ooOO00OOoo....
    /*
    G4RotationMatrix *RotCal2 = new G4RotationMatrix();
    RotCal2->rotateZ(-90 * deg);
    RotCal2->rotateX(90 * deg);
    RotCollector[0] = RotCal2;

    new G4PVPlacement(
        RotCal2,                                                   // rotation
        G4ThreeVector(-(Cal_Thickness + RPC_BoxLength) / 2, 0, 0), //
        CalLV[2],                                                  //its logical volume
        "Cal2",                                                    //its name
        logicWorld,                                                //its mother  volume
        false,                                                     //no boolean operation
        2,                                                         //copy number
        true);                                                     //check overlaps

    G4RotationMatrix *RotCal3 = new G4RotationMatrix();
    RotCal3->rotateZ(90 * deg);
    RotCal3->rotateX(90 * deg);
    RotCollector[1] = RotCal3;

    new G4PVPlacement(
        RotCal3,                                                  // rotation
        G4ThreeVector((Cal_Thickness + RPC_BoxLength) / 2, 0, 0), //
        CalLV[2],                                                 //its logical volume
        "Cal3",                                                   //its name
        logicWorld,                                               //its mother  volume
        false,                                                    //no boolean operation
        3,                                                        //copy number
        true);                                                    //check overlaps

    G4RotationMatrix *RotCal4 = new G4RotationMatrix();
    RotCal4->rotateX(270 * deg);
    RotCollector[2] = RotCal4;

    new G4PVPlacement(
        RotCal4,                                                   // rotation
        G4ThreeVector(0, -(Cal_Thickness + RPC_BoxLength) / 2, 0), //
        CalLV[2],                                                  //its logical volume
        "Cal4",                                                    //its name
        logicWorld,                                                //its mother  volume
        false,                                                     //no boolean operation
        4,                                                         //copy number
        true);                                                     //check overlaps

    G4RotationMatrix *RotCal5 = new G4RotationMatrix();
    // RotCal5->rotateZ(90*deg);
    RotCal5->rotateX(90 * deg);
    RotCollector[3] = RotCal5;

    new G4PVPlacement(
        RotCal5,                                                  // rotation
        G4ThreeVector(0, (Cal_Thickness + RPC_BoxLength) / 2, 0), //
        CalLV[2],                                                 //its logical volume
        "Cal5",                                                   //its name
        logicWorld,                                               //its mother  volume
        false,                                                    //no boolean operation
        5,                                                        //copy number
        true);                                                    //check overlaps
*/
    logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

    return physiWorld; //always return the physical World
}

//G4VisAttributes* mylarAttributes = new G4VisAttributes(G4Colour::Green());
//mylarAttributes->SetForceSolid(true);
//LogicalVolumes[_RPC_mylarCathode]->SetVisAttributes(mylarAttributes);
//LogicalVolumes[_RPC_mylarAnode]->SetVisAttributes(mylarAttributes);

//G4VisAttributes* electronicAttributes = new G4VisAttributes(G4Colour::Yellow());
//electronicAttributes->SetForceSolid(true);
//LogicalVolumes[_RPC_Electronics]->SetVisAttributes(electronicAttributes);
//LogicalVolumes[_RPC_PCB]->SetVisAttributes(electronicAttributes);

// LogicalVolumes[_RPC_AdditionalAlu]->SetVisAttributes(absorberAttributes);
// LogicalVolumes[_RPC_GraphiteCathode]->SetVisAttributes(graphiteAttributes);
// LogicalVolumes[_RPC_GraphiteAnode]->SetVisAttributes(graphiteAttributes);

//  AirSDetector = new LyonAirSensitiveDetector("AirSD");
//  AirSDetector->SetFilter(MuonFilter);
//  G4SDManager::GetSDMpointer()->AddNewDetector(AirSDetector);
//  NbOfCal = 6;
//  NbOfRPC = 6;
//  RPC_Distance = 6.0*cm;
//  RPC_BoxLength = 5*m;
//  CalDistance = 3*m;
//  WorldFreeSpace = 50*cm;
//  CalThickness = fCalorimeter.GetThickness();
//  DetectorSizeZ = 2*CalThickness + CalDistance;
//  WorldSizeY = RPC_BoxLength + CalThickness * 2+ WorldFreeSpace;
//  WorldSizeX = WorldSizeY ;
//  WorldSizeZ = DetectorSizeZ + WorldFreeSpace;

//#include "G4VProcess.hh"
//#include "G4VSDFilter.hh"
//#include "G4SubtractionSolid.hh"
//#include "G4TransportationManager.hh"
//#include "G4FieldManager.hh"
//#include "G4UniformMagField.hh"
//#include "CalParameterisation.hh"
//#include "G4PVParameterised.hh"
//#include "G4Tubs.hh"
//#include "G4PVReplica.hh"
//#include "LyonDetectorFieldSetup.hh"

/*
  bool target = 1; //Whether or not there is a target
  bool targetIn = 1; //If it's a screed

  //The target size
  G4double TargetSizeX = 30.0*CLHEP::cm;
  G4double TargetSizeY = 30.0*CLHEP::cm;
  G4double TargetSizeZ = 30*CLHEP::cm; 
  G4double TargetThickness = 5.0*CLHEP::cm;
  G4double TargetDistance = -1.0*CLHEP::m;
  */

// Calorimeter Parameters

/*	
  // Calorimeter Out
  G4int CaloOutNbLayer = 6;
  G4double CaloOutLayerDistance = 20.0*CLHEP::cm;
  G4double CaloOutRadius = 4*2*37*CLHEP::cm/2;
	
  G4double CaloDistance = 400.0*CLHEP::cm;
  */

/* = Calorimeters =================================================== */

/* The following need some explanations : we need some calorimeters
   * informations mostly geometrical in other parts of the code, so I
   * fill a static array */

/*
  fCalorimeters.push_back(LyonCalorimeter(
  0, //ID 
    CaloInNbLayer, //NbLayers
    CaloInLayerDistance, //Layers spacement
    SDetector,
    CaloInRadius));
  */
/*	fCalorimeters.push_back(LyonCalorimeter(
	1, //ID
	CaloOutNbLayer, //NbLayers
	CaloOutLayerDistance, //Layers spacement
	detector,
	CaloOutRadius));
  */

/*	
  fGeometryData.Area = (2*CaloInRadius/CLHEP::m)*(2*CaloInRadius/CLHEP::m);
  fGeometryData.WorldSizeZ = WorldSizeZ;
  fGeometryData.Scintillators.resize(3);
	
  for(unsigned int i(0); i<fGeometryData.Scintillators.size(); ++i)
    fGeometryData.Scintillators[i].resize(2);
	
  G4ThreeVector SizeTall = G4ThreeVector(0.15*CLHEP::m,0.25*CLHEP::m,0.0); //m
  G4ThreeVector SizeSmall = G4ThreeVector(0.07*CLHEP::m,0.25*CLHEP::m,0.0); //m
	
  fGeometryData.Scintillators[0][0]=G4ThreeVector(-3.5*CLHEP::cm,0.0,3*.140*CLHEP::m); //m
  fGeometryData.Scintillators[1][0]=G4ThreeVector(-1.5*CLHEP::cm,0.0,-0.140*CLHEP::m); //m
  fGeometryData.Scintillators[2][0]=G4ThreeVector(0.0,0.0,-3*0.140*CLHEP::m); //m
	
  fGeometryData.Scintillators[0][1]=SizeSmall;
  fGeometryData.Scintillators[1][1]=SizeTall;
  fGeometryData.Scintillators[2][1]=SizeTall;
  */

/*===========================Put Calorimeters in the world===========================*/
//  G4double PositionZ = (CalDistance + fCalorimeter->GetThickness())/2;
/*
  //The first one, on the bottom
  new G4PVPlacement(  
		    0,                                               //no rotation
		    G4ThreeVector(0, 0, -PositionZ),                 //at (0,0,0)
		    fCalorimeter->GetLogicalCalorimeter(),           //its logical volume   
		    "Cal0",                                          //its name
		    logicWorld,                                      //its mother  volume
		    false,                                           //no boolean operation
		    0,                                               //copy number
		    true);                                           //Checking overlaps
  //The Second one, in the air
  new G4PVPlacement(  
		    0,                                               //no rotation
		    G4ThreeVector(0, 0, PositionZ),                  //at (0,0,0)
		    fCalorimeter->GetLogicalCalorimeter(),           //its logical volume   
		    "Cal1",                                          //its name
		    logicWorld,                                      //its mother  volume
		    false,                                           //no boolean operation
		    1,                                               //copy number
		    true);                                           //Checking overlaps
  */
/* 
  G4VPVParameterisation* CalParam =
    new CalParameterisation(NbOfCal,   // NoChambers
			    CalDistance,  // Z of center of first
			    RPC_BoxLength, // Z spacing of centers
			    fCalorimeter->GetThickness());   // chamber width
			 
			   
  // dummy value : kZAxis -- modified by parameterised volume

  new G4PVParameterised("Calorimeter",       // their name
                        fCalorimeter->GetLogicalCalorimeter(),   // their logical volume
                        logicWorld,       // Mother logical volume
                        kZAxis,          // Are placed along this axis 
                        NbOfCal,    // Number of chambers
                        CalParam,    // The parametrisation
                        true); // checking overlaps 
  */
/* = Target ========================================================= */

/*	G4Box *solidTarget = new G4Box(
	"Target",     //its name
	TargetSizeX/2,    //size
	TargetSizeY/2,
	TargetSizeZ/2);
		
	G4LogicalVolume *logicTarget = new G4LogicalVolume(
	solidTarget,      //its solid
	Pb, //its material
	"Target");    //its name
	 
	G4Box *solidTargetIn = new G4Box(
	"Target",     //its name
	(TargetSizeX - 2*TargetThickness)/2 ,    //size
	(TargetSizeY - 2*TargetThickness)/2,
	(TargetSizeZ - 2*TargetThickness)/2);
		
	G4LogicalVolume *logicTargetIn = new G4LogicalVolume(
	solidTargetIn,      //its solid
	Air, //its material
	"TargetIn");    //its name 

	if(target)
	{
	new G4PVPlacement(
	0,          //no rotation
	G4ThreeVector(0.5*CLHEP::m,-0.75*CLHEP::m, TargetDistance), 
	logicTarget,  //its logical volume
	"Target",     //its name
	logicWorld,     //its mother  volume
	false,        //no boolean operation
	0,
	true);        //copy number
	}
	if(target && targetIn)
	{
	new G4PVPlacement(
	0,          //no rotation
	G4ThreeVector(0,0, 0), 
	logicTargetIn,  //its logical volume
	"Target",     //its name
	logicTarget,     //its mother  volume
	false,        //no boolean operation
	0,
	true);        //copy number
	}*/

// Create calorimeters and make them accesible because static !

//  fCalorimeters[0].PlaceIn(logicWorld, /*(CaloDistance+fCalorimeters[0].GetThickness()/2*/0.0);

//	fCalorimeters[1].PlaceIn(logicWorld, -(CaloDistance+fCalorimeters[1].GetThickness())/2);

/* = Visualization attributes ======================================= */

//	G4VisAttributes* targetAttributes = new G4VisAttributes(G4Colour::Grey());
//	targetAttributes->SetForceSolid(true);
//logicTargetIn->SetVisAttributes(targetAttributes);
//	logicTarget->SetVisAttributes(G4Colour::Grey());
//logicWorld->SetVisAttributes (G4Colour::Grey());

//#define PI 3.1415926

//GeometryData LyonDetectorConstruction::fGeometryData=GeometryData();

//std::vector<LyonCalorimeter> LyonDetectorConstruction::fCalorimeters;
/*
//helper print function
void printInterval(ostream& flux, const char* name, G4double middle, G4double halfsize)
{
  flux << name << " " << middle << " middle of [" << middle-halfsize << " , " << middle+halfsize <<"]" << endl;
}
*/

//#define MAGNETIC_FIELD

/*#ifdef MAGNETIC_FIELD
    G4MagneticField* magField = new G4UniformMagField(G4ThreeVector(0.0,3.0*tesla,0.0));
    G4FieldManager* globalFieldManager=G4TransportationManager::GetTransportationManager()->GetFieldManager();
    globalFieldManager->SetDetectorField(magField);
    globalFieldManager->CreateChordFinder(magField);
    #endif*/

//Get standard stuff from NIST database :
//fEmFieldSetup = new LyonDetectorFieldSetup();

/*
  G4Box* AirLayerSV = new G4Box(
				"AirSV",
				RPC_BoxLength/2,
				RPC_BoxLength/2,
			        0.01*mm);

  G4LogicalVolume *AirLayerLV = new G4LogicalVolume(
						    AirLayerSV,     //its solid
						    Air,  //its material
						    "AirLayerLV");       //its name

  new G4PVPlacement(  
		    0,                //no rotation
		    G4ThreeVector(0., 0., 1.4*m),  //at (0,0,0)
		    AirLayerLV,       //its logical volume   
		    "AirLayerPV",          //its name
		    logicWorld,                //its mother  volume
		    false,            //no boolean operation
		    0,                //copy number
		    true);            //check overlaps
  new G4PVPlacement(  
		    0,                //no rotation
		    G4ThreeVector(0., 0., 0.),  //at (0,0,0)
		    AirLayerLV,       //its logical volume   
		    "AirLayerPV",          //its name
		    logicWorld,                //its mother  volume
		    false,            //no boolean operation
		    1,                //copy number
		    true);            //check overlaps
  new G4PVPlacement(  
		    0,                //no rotation
		    G4ThreeVector(0., 0., -1.4*m),  //at (0,0,0)
		    AirLayerLV,       //its logical volume   
		    "AirLayerPV",          //its name
		    logicWorld,                //its mother  volume
		    false,            //no boolean operation
		    2,                //copy number
		    true);            //check overlaps

		    AirLayerLV->SetSensitiveDetector(AirSDetector);
*/
