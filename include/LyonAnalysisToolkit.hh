#ifndef LYONANALYSISTOOLKIT_HH
#define LYONANALYSISTOOLKIT_HH

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "G4ThreeVector.hh"
#include "LyonTrackHit.hh"

/*!
 * Class LinearFit3D.
 * Fit a 3D line in 3D space. The line equation is the following one :
 *       x = a*z + b
 *       y = c*z + d
 *       z = t
 */

class LinearFit3D 
{
protected :

  G4double _chi2;
  G4double _params[4];
  G4double _paramsError[4];
  void ComputeChi2();
  G4ThreeVector _direction;
  G4ThreeVector _point;

public :
  LinearFit3D( const std::vector<G4ThreeVector> & _positions);
  virtual ~LinearFit3D(){;}
  inline G4double getChi2(){return _chi2;}
  inline G4double* getFitParameters(){return _params;}

  inline G4double* getFitParError(){return _paramsError;}
  inline G4ThreeVector GetDirection(){return _direction;}
  inline G4ThreeVector GetPoint(){return _point;}
};
 
class LyonAnalysisToolkit
{
private :
  LyonTrackHitsCollection* hitsCollection;

  G4int NbOfTracks;
  G4int Nb_Cal;

  
  std::vector< std::vector< G4ThreeVector > > Dir;

  std::vector< G4ThreeVector > Point;
  std::vector<G4double> dAngle;

  std::vector<std::vector<std::vector<G4int> > >IDList;
  std::vector<std::vector<G4int> >CalList;

public :

  LyonAnalysisToolkit(LyonTrackHitsCollection* hc0, G4int NumberOfTracks);
    
  G4int Analyze(G4int TrackNb);
  bool GetFit(G4int TNb);
  void GetSingleFit(G4int TNb);

  inline const G4ThreeVector& GetPoint(G4int TrackNb) const {return Point[TrackNb];};
  inline G4double GetDAngle(G4int TNb) const {return dAngle[TNb];};
  inline const G4ThreeVector& GetSingleDir(G4int TNb) const { return Dir[TNb][0]; }

 

  inline G4int GetNbOfTracks() const {return NbOfTracks;}

};


#endif

// inline const std::vector<G4double> & GetTracks(G4int TNb) const { return Tracks[TNb]; }
  //std::vector< std::vector< G4double > >Tracks;
  //  std::vector<std::vector<std::vector<std::vector<LyonTrackHit*> > > >hits;
  //  std::vector<std::vector<int> > singleHits;
  //  std::vector< G4ThreeVector > DirOut;
  //  G4ThreeVector DirTot;
  //  LinearFit3D GetFit();
  //  G4ThreeVector GetDirOut(){return DirOut;}
  //  G4ThreeVector GetDirIn(){return DirIn;}
  //  G4ThreeVector GetDirTot(){return DirTot;}
  //  std::vector<G4int> TrackIDList;
  //  G4double GetPhi(G4int TrackNb, G4int inout);w
  //  inline G4int GetTrackID(G4int nb) const {return TrackIDList[nb];};
  //  std::vector<G4ThreeVector> _positions;
  //G4double getFitParameters(int i){return _params[i];}
