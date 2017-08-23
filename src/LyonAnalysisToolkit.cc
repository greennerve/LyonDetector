#include "LyonAnalysisToolkit.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
//#include "LyonDetectorConstruction.hh"

//#include "LyonPrimaryGeneratorAction.hh"

#define PI 3.1415926
#define Rad2Degree 57.2957795131
//#define Nb_Cal 2

LinearFit3D::LinearFit3D(const std::vector<G4ThreeVector> &_positions)
{
  _chi2 = _params[0] = _params[1] = _params[2] = _params[3] = 0;
  _paramsError[0] = _paramsError[1] = _paramsError[2] = _paramsError[3] = 0;
  //  _positions = pos;
  _params[0] = _params[1] = _params[2] = _params[3] = 0;

  G4double xsum = 0.0;
  G4double ysum = 0.0;
  G4double zsum = 0.0;
  G4double zzsum = 0.0;
  G4double xzsum = 0.0;
  G4double yzsum = 0.0;

  for (unsigned int i = 0; i < _positions.size(); i++)
  {
    //for equation 1
    zsum = zsum + _positions.at(i).z();
    xsum = xsum + _positions.at(i).x();
    zzsum = zzsum + (_positions.at(i).z() * _positions.at(i).z());
    xzsum = xzsum + _positions.at(i).x() * _positions.at(i).z();

    //for equation 2
    ysum = ysum + _positions.at(i).y();
    yzsum = yzsum + _positions.at(i).y() * _positions.at(i).z();
  }

  G4double A1 = zsum;
  G4double B1 = _positions.size();
  G4double C1 = xsum;
  G4double D1 = zzsum;
  G4double E1 = xzsum;

  G4double C2 = ysum;
  G4double E2 = yzsum;

  _params[0] = (D1 * C1 - E1 * A1) / (B1 * D1 - A1 * A1);
  _params[1] = (E1 * B1 - C1 * A1) / (B1 * D1 - A1 * A1);
  _params[2] = (D1 * C2 - E2 * A1) / (B1 * D1 - A1 * A1);
  _params[3] = (E2 * B1 - C2 * A1) / (B1 * D1 - A1 * A1);

  _paramsError[0] = sqrt(D1 / (B1 * D1 - A1 * A1));
  _paramsError[1] = sqrt(B1 / (B1 * D1 - A1 * A1));
  _paramsError[2] = sqrt(D1 / (B1 * D1 - A1 * A1));
  _paramsError[3] = sqrt(B1 / (B1 * D1 - A1 * A1));

  _direction = G4ThreeVector(_params[1], _params[3], 1);

  _point = G4ThreeVector(_params[0], _params[2], 0.0);

  //  ComputeChi2();
}

void LinearFit3D::ComputeChi2()
{

  _chi2 = 0;
}

LyonAnalysisToolkit::LyonAnalysisToolkit(LyonTrackHitsCollection *hc0, G4int NumberOfTracks)
    : hitsCollection(hc0),
      NbOfTracks(NumberOfTracks)
//    Tracks(NbOfTracks, std::vector<G4double>(13, 0))
{
  // 2 : Number of calorimeters -> if you want to change to N, interacts with static members of geometry
  IDList.clear();
  IDList.resize(NbOfTracks);

  CalList.clear();
  CalList.resize(NbOfTracks);

  for (G4int i = 0; i < NbOfTracks; i++)
  {
    IDList[i].clear();
    IDList[i].resize(2);
  }

  for (G4int i = 0; i < hitsCollection->entries(); ++i)
  {
    LyonTrackHit *hit = (*hitsCollection)[i];
    IDList[hit->GetTrackID() - 1][hit->GetCalo()].push_back(i);
  }

  Dir.clear(); // Direction "momentum" of entering particle
  Dir.resize(NbOfTracks);

  //  DirTot = G4ThreeVector(); // Used to build total direction if no target, ie no deviation
  Point.clear();
  Point.resize(NbOfTracks);
  dAngle.clear();
  dAngle.resize(NbOfTracks);
}

G4int LyonAnalysisToolkit::Analyze(G4int TNb)
{

  CalList[TNb].clear();

  G4int Touched = 0;

  for (G4int i = 0; i < 2; i++)
  {
    if (IDList[TNb][i].size() > 0)
    {
      Touched++;

      if (IDList[TNb][i].size() == 3)
        CalList[TNb].push_back(i);
    }
  }

  if ((CalList[TNb].size() == 2) && (Touched == 2))
    return 2;

  if ((CalList[TNb].size() == 1) && (Touched == 1))
  {
    if (CalList[TNb][0] != 1)
    {
      return 1;
    }
  }

  return 0;
}

bool LyonAnalysisToolkit::GetFit(G4int TNb)
{
  std::vector<std::vector<G4ThreeVector>> pos;
  pos.clear();
  pos.resize(2);

  Dir[TNb].clear();
  Dir[TNb].resize(2);

  for (G4int i = 0; i < 2; i++)
  {
    G4int CalID = CalList[TNb][i];
    for (unsigned j = 0; j < IDList[TNb][CalID].size(); j++)
    {
      G4int HitID = IDList[TNb][CalID][j];
      LyonTrackHit *hit = (*hitsCollection)[HitID];
      pos[i].push_back(hit->GetExitPoint());
    }
  }

  LinearFit3D *Fit0 = new LinearFit3D(pos[0]);
  LinearFit3D *Fit1 = new LinearFit3D(pos[1]);

  Dir[TNb][0] = Fit0->GetDirection();
  Dir[TNb][1] = Fit1->GetDirection();

  G4ThreeVector Pt0(Fit0->GetPoint());
  G4ThreeVector Pt1(Fit1->GetPoint());

  dAngle[TNb] = Dir[TNb][0].angle(Dir[TNb][1]) * Rad2Degree;

  G4ThreeVector u(Fit0->GetDirection());
  G4ThreeVector v(Fit1->GetDirection());
  G4ThreeVector w = Pt0 - Pt1;
  G4double a = u.dot(u); // always >= 0
  G4double b = u.dot(v);
  G4double c = v.dot(v); // always >= 0
  G4double d = u.dot(w);
  G4double e = v.dot(w);
  G4double D = a * c - b * b; // always >= 0
  G4double sc, tc;

  // compute the line parameters of the two closest points
  if (D < 0.00000001) // the lines are almost parallel
  {
    sc = 0.0;
    tc = (b > c ? d / b : e / c); // use the largest denominator
  }
  else
  {
    sc = (b * e - c * d) / D;
    tc = (a * e - b * d) / D;
  }

  u *= sc;
  Pt0 += u;
  v *= tc;
  Pt1 += v;

  if ((Pt0 - Pt1).mag() <= 100)
  {
    G4double x = (Pt0.getX() + Pt1.getX()) / 2;
    G4double y = (Pt0.getY() + Pt1.getY()) / 2;
    G4double z = (Pt0.getZ() + Pt1.getZ()) / 2;

    if (abs(x) < 2500 && abs(y) < 2500 && abs(z) < 1500)
    {
      Point[TNb] = G4ThreeVector(x, y, z);

      if (Fit0)
        delete Fit0;
      if (Fit1)
        delete Fit1;

      return true;
    }
  }

  if (Fit0)
    delete Fit0;
  if (Fit1)
    delete Fit1;

  return false;
}

void LyonAnalysisToolkit::GetSingleFit(G4int TNb)
{
  std::vector<G4ThreeVector> pos;
  pos.clear();

  Dir[TNb].clear();
  Dir[TNb].resize(1);

  G4int CalID = CalList[TNb][0];

  for (unsigned j = 0; j < IDList[TNb][CalID].size(); j++)
  {
    G4int HitID = IDList[TNb][CalID][j];
    LyonTrackHit *hit = (*hitsCollection)[HitID];
    pos.push_back(hit->GetExitPoint());
  }

  LinearFit3D *Fit0 = new LinearFit3D(pos);

  Dir[TNb][0] = Fit0->GetDirection();
  Point[TNb] = Fit0->GetPoint();

  if (Fit0)
    delete Fit0;
}

/*
      for (G4int i = 0; i < 3; i++)
      {
        Tracks[TNb][i] = Pt0[i];
      }

      for (G4int i = 0; i < 3; i++)
      {
        Tracks[TNb][i + 3] = Dir[TNb][0][i];
      }

      for (G4int i = 0; i < 3; i++)
      {
        Tracks[TNb][i + 6] = Pt1[i];
      }

      for (G4int i = 0; i < 3; i++)
      {
        Tracks[TNb][i + 9] = Dir[TNb][i][i];
      }

      Tracks[TNb][12] = dAngle[TNb];
      */

// This functions return a state as int :
// * -1 : enough hit to build a track acrross both detector (ie no target)
//  * 0 : not enough
// * 1 : enough to build in
//  * 2 : enough to build in/out
//  * Condition of building is simply to get a single hit per event per
//  * layer and at least a certains nb of hits per calo
// * This functions also realize the fit of tracks if it can, that
// * then you acess in EventAction knowing the state ;)

/*
G4double LyonAnalysisToolkit::GetPhi(G4int TrackNb, G4int inout)
{
  G4ThreeVector Dir = G4ThreeVector();
  G4ThreeVector Normal = G4ThreeVector(-1.0,0.0,0.0); 
	
  if(inout ==0)
    Dir = DirIn[TrackNb];
  else if(inout ==1)
    Dir = DirOut[TrackNb];
  //  else if(inout ==-1)
  //    Dir = DirTot;

	
  Dir.setZ(0.0);
  Dir = Dir.unit();

  return acos(Normal.dot(Dir))*(Dir.getY()/fabs(Dir.getY()))*Rad2Degree;
}

*/

/*
  
  for(G4int j=0; j < Nb_Cal; ++j) // 2 ! -> N calo
  for(unsigned i=0; i < hits[TrackNb][j].size(); ++i)
  if(hits[TrackNb][j][i].size()>0)	singleHits[TrackNb][j] ++;

  G4cout<<"Number of entering muon hits: "<<singleHits[TrackNb][0]<<G4endl;
  G4cout<<"Number of leaving muon hits: "<<singleHits[TrackNb][1]<<G4endl;

  // Has to be fixed more cleverly
  //  int minHitsTot = 3;
  int minHits = 4;
  int state =-1;
	
  LinearFit3D *fitIn = NULL;
  LinearFit3D *fitOut = NULL;
	
  if(singleHits[TrackNb][0]>=minHits && singleHits[TrackNb][1]>=minHits)
  {
  //      state = 1; // Entering particle !
  fitIn = GetFit(TrackNb,0);
  DirIn[TrackNb] = fitIn->GetDirection();		
 
  fitOut = GetFit(TrackNb,1);
  DirOut[TrackNb] = fitOut->GetDirection();
  state = 2;  
  }


  if(singleHits[TrackNb][1] == 0)
  {
  //      fitOut = GetFit(1);
  //     DirOut = fitOut->GetDirection();
  //      fitIn = GetFit(0);
  //      DirIn = fitIn->GetDirection();
  //           G4cout<<"MX "<< DirIn.getX()<<"     "
  //        	    <<"MY "<< DirIn.getY()<<"     "
  //  	    <<"MZ "<< DirIn.getZ()<<"     "
  // 	    <<G4endl;
  G4RunManager* runManager = G4RunManager::GetRunManager();
  LyonPrimaryGeneratorAction* _LyonPrimaryGeneratorAction = (LyonPrimaryGeneratorAction*) runManager->GetUserPrimaryGeneratorAction();
  G4ThreeVector P = _LyonPrimaryGeneratorAction->Posit;
  G4ThreeVector M = _LyonPrimaryGeneratorAction->Momen;
  G4cout<<"MX "<< M.getX()<<"     "
  <<"MY "<< M.getY()<<"     "
  <<"MZ "<< M.getZ()<<"     "
  <<G4endl;

  G4cout<<"PX "<< P.getX()<<"     "
  <<"PY "<< P.getY()<<"     "
  <<"PZ "<< P.getZ()<<"     "
  <<G4endl;
  //     delete fitOut;
  //      fitOut = NULL;
  }	

  if(state==2)
  {
  // Determine the "interaction point"

  //....ooOO0OOoo........ooOO0OOoo........ooOO0OOoo....
  //Maxime's simple approche


  */

/*
  G4double* t1 = fitIn->getFitParameters();
  G4double* t2 = fitOut->getFitParameters();
		
  G4double zx = - (t1[0]-t2[0])/(t1[1]-t2[1]);
  G4double x = t1[1]*zx + t1[0];
  G4double zy = - (t1[2]-t2[2])/(t1[3]-t2[3]);
  G4double y = t1[3]*zy + t1[2];      		
  Point[TrackNb] = G4ThreeVector(x,y,(zx+zy)*0.5);
*/
//....ooOO0OOoo........ooOO0OOoo........ooOO0OOoo....

//The closest point approch

/*
  
  G4ThreeVector u(fitIn->GetDirection());
  G4ThreeVector v(fitOut->GetDirection());
  G4ThreeVector w = fitIn->GetPoint() - fitOut->GetPoint();
  G4double    a = u.dot(u);         // always >= 0
  G4double    b = u.dot(v);
  G4double    c = v.dot(v);         // always >= 0
  G4double    d = u.dot(w);
  G4double    e = v.dot(w);
  G4double    D = a*c - b*b;        // always >= 0
  G4double    sc, tc;

  // compute the line parameters of the two closest points
  if (D < 0.00000001)  // the lines are almost parallel
  {         
  sc = 0.0;
  tc = (b>c ? d/b : e/c);    // use the largest denominator
  }
  else
  {
  sc = (b*e - c*d) / D;
  tc = (a*e - b*d) / D;
  }
      
  u *= sc;
  u += fitIn->GetPoint();
  v *= tc;
  v += fitOut->GetPoint();
      
  double x = (u.getX()+v.getX())/2;      
  double y = (u.getY()+v.getY())/2;
  double z = (u.getZ()+v.getZ())/2;
		
  Point[TrackNb] = G4ThreeVector(x,y,z);
      
  }
	
  if(fitOut) delete fitOut;
  if(fitIn) delete fitIn;

  return state;
*/
/*	
  	if(singleHits[0]+singleHits[1]>=minHitsTot)
	{
	state =0;
	DirTot= this->GetFit().GetDirection();
	}
*/
//}

/* Folowing functions are called, after we check there was enough hits to
 * build track, a more clever way to make it could be to build a function
 * bool GetFit(LinearFit3D, std::vector<std::vector<int> > ) returning 
 * fit state, and a array[calo ID][layers ID] specifying layer use to build
 * tracks */

/* Note that in more complex simulation user should maybe forget about 
 * Calorimeter class to move on Layer class, it maybe should be more handy 
 * to work with */

//  G4RunManager* runManager = G4RunManager::GetRunManager();
//  LyonPrimaryGeneratorAction* _LyonPrimaryGeneratorAction = (LyonPrimaryGeneratorAction*) runManager->GetUserPrimaryGeneratorAction();

//  NbOfTracks = _LyonPrimaryGeneratorAction->GetNbOfPrimaryParticles();
/*
  Nb_Cal = 6;
  TrackIDList.clear();
  for(G4int i =0; i < hitsCollection->entries(); ++i)
  {
  bool TrackListUnlocker = 1;
  LyonTrackHit* hit = (*hitsCollection)[i];
  for(unsigned j = 0; j < TrackIDList.size(); j++)
  {
  if(TrackIDList[j] == hit->GetTrackID())
  {
  TrackListUnlocker = 0;
  break;
  }
  }
  if(TrackListUnlocker) TrackIDList.push_back(hit->GetTrackID());
  }
  NbOfTracks = TrackIDList.size();
  
  hits.clear();
  hits.resize(NbOfTracks); 
  singleHits.clear();
  singleHits.resize(NbOfTracks);



  for(G4int k=0; k < NbOfTracks; k++)
  {
  hits[k].resize(Nb_Cal);
  singleHits[k].resize(Nb_Cal);
  for(unsigned i =0; i < hits[k].size(); ++i)
  {
  hits[k][i].resize(6);
  for(unsigned  j =0; j < hits[k][i].size(); ++j)
  hits[k][i][j].clear();
  }
  }
  
  
  for(G4int i =0; i < hitsCollection->entries(); ++i)
  {

  LyonTrackHit* hit = (*hitsCollection)[i];

  for(G4int k = 0; k < NbOfTracks; k++)
  {
  if(TrackIDList[k] == hit->GetTrackID())
  {
  hits[k][hit->GetCalo()][hit->GetChamber()].push_back(hit);
  break;
  }
  }
  }
*/
//  if(NbOfTracks != TrackIDList.size() && TrackIDList.size() != 0)
//    {
//      G4String* str = new G4String("Number of tracks doesn't match");
//      G4Exception("LyonAnalysisToolKit", "1", RunMustBeAborted, *str);
//    }

/*algorithm::Distance<CLHEP::Hep3Vector,G4double *> dist;
    for( unsigned int i=0 ; i<_positions.size() ; i++ ) {
    G4double d=dist.getDistance(_positions.at(i),_params);
    G4double mult = _clusterSize.at(i);
    G4double err = mult/sqrt(12/100.);
    _chi2 += (d/err)*(d/err);
    }
    _chi2=_chi2/(_positions.size()-1);*/
