//
////******************************************************************************
//what's wrong with this editor??


#ifndef LyonConstants_h
#define LyonConstants_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

constexpr G4int NbOfCal = 2;
constexpr G4int NbOfRPC = 3;
constexpr G4double  RPC_Distance = 25.0*cm;
constexpr G4double  RPC_BoxLength = 5*m;
constexpr G4double  CalDistance = 3*m;
constexpr G4double  WorldFreeSpace = 20*cm;
  
constexpr G4double RPC_Absorber_Thickness = 0.2*cm; 
constexpr G4double RPC_AdditionalAlu_Thickness = 0.3*cm; //For K7 cover on the thick glass side
constexpr G4double RPC_Inner_Thickness = 6.131*mm;//by rhan from 6.00 to 6.031
constexpr G4double RPC_PCB_Thickness = 1.200*mm;//by rhan from 0.8 to 1.2
constexpr G4double RPC_mylarAnode_Thickness = 0.050*mm;
constexpr G4double RPC_mylarCathode_Thickness = 0.180*mm;//by rhan from 0.2 to 0.18
constexpr G4double RPC_GraphiteAnode_Thickness = 0.050*mm;
constexpr G4double RPC_GraphiteCathode_Thickness = 0.050*mm;//by rhan from 0.1 to 0.05
constexpr G4double RPC_ThinGlass_Thickness = 0.700*mm;
constexpr G4double RPC_GasChamber_Thickness = 1.200*mm;
constexpr G4double RPC_ThickGlass_Thickness = 1.100*mm;
constexpr G4double RPC_ChipPackage_Thickness = 1.600*mm;//by rhan from 1.40 to 1.60
constexpr G4double RPC_Free_Thickness =   RPC_Inner_Thickness
                                        - RPC_ChipPackage_Thickness
                                        - RPC_PCB_Thickness
                                        - RPC_mylarAnode_Thickness
                                        - RPC_mylarCathode_Thickness
                                        - RPC_GraphiteAnode_Thickness
                                        - RPC_GraphiteCathode_Thickness
                                        - RPC_ThinGlass_Thickness
                                        - RPC_ThickGlass_Thickness
                                        - RPC_GasChamber_Thickness;

constexpr G4double RPC_Thickness =   RPC_Absorber_Thickness
                                   + RPC_Inner_Thickness 
                                   + RPC_AdditionalAlu_Thickness;//Total RPC thickness include the absorber and additional aluminum

constexpr G4double Cal_Thickness = NbOfRPC * (RPC_Thickness + RPC_Distance) - RPC_Distance;
constexpr G4double DetectorSizeZ = 2*Cal_Thickness + CalDistance;

//constexpr G4double WorldSizeY = RPC_BoxLength + Cal_Thickness * 2+ WorldFreeSpace;
constexpr G4double WorldSizeY = 10*m;
constexpr G4double WorldSizeX = WorldSizeY ;
constexpr G4double WorldSizeZ = DetectorSizeZ + WorldFreeSpace;

#endif
 
