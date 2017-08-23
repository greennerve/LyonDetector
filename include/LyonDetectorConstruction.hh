#ifndef LyonDetectorConstruction_H
#define LyonDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"

#include "globals.hh"



class G4VPhysicalVolume;
class G4LogicalVolume;

class G4SDParticleFilter;
class LyonSensitiveDetector;
struct VolumeInfo;


class LyonDetectorConstruction : public G4VUserDetectorConstruction
{
private:

  enum _Layers {
    _RPC_Absorber         = 0,
    _RPC_Electronics      = 1,
    _RPC_PCB              = 2,
    _RPC_mylarAnode       = 3,
    _RPC_GraphiteAnode    = 4,
    _RPC_ThinGlass        = 5,
    _RPC_GasChamber       = 6,
    _RPC_ThickGlass       = 7,
    _RPC_GraphiteCathode  = 8,
    _RPC_mylarCathode     = 9,
    _RPC_Free             = 10,
    _RPC_AdditionalAlu    = 11,
    _RPC                  = 12,
    _Calorimeter          = 13,
    _NbOfVolumes          = 14
  };

  std::vector< VolumeInfo* > fVolumeInfo;
  std::vector<G4LogicalVolume*> CalLV;
  std::vector<G4LogicalVolume*> LogicalVolumes;

  std::vector<G4RotationMatrix*> RotCollector;
  std::vector<G4VisAttributes*> VisCollector;

  G4SDParticleFilter* MuonFilter;
  LyonSensitiveDetector* SDetector;

public:

  LyonDetectorConstruction();
  ~LyonDetectorConstruction();
  
  G4VPhysicalVolume* Construct();

  void DefineMaterials();
  void DefineVolumeInfo();

  G4LogicalVolume* ConstructCalorimeter(G4double XHalfLength, G4double YhalfLength);

  //  inline G4double GetCalBoxSize() const { return RPC_BoxLength;}
  //  inline G4double GetDetectorSizeZ() const { return DetectorSizeZ;}
  //  inline G4double GetWorldSizeZ() const { return WorldSizeZ;}
  //  inline G4double GetCalDistance() const { return CalDistance;}
};

#endif

/*

  G4int NbOfCal;
  G4int NbOfRPC;
  G4double RPC_Distance;
  G4double RPC_BoxLength;
  G4double CalDistance;
  G4double CalThickness;
  G4double DetectorSizeZ;
  G4double WorldFreeSpace;
  G4double WorldSizeY;
  G4double WorldSizeX;
  G4double WorldSizeZ;  

*/
