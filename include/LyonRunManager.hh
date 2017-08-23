#ifndef LyonRunManager_h
#define LyonRunManager_h 1
#include "G4RunManager.hh"
#include "LyonPrimaryGeneratorAction.hh"

class LyonRunManager : public G4RunManager
{
  public:
	LyonRunManager(){lyonRunMan = this;}
	~LyonRunManager(){}
	LyonPrimaryGeneratorAction*  GetLyonPrimaryGeneratorAction(){return lyonPrimaryGeneratorAction;}
	void SetLyonPrimaryGeneratorAction(LyonPrimaryGeneratorAction* priAct){lyonPrimaryGeneratorAction = priAct;}
	static LyonRunManager* GetLyonRunMan(){return lyonRunMan;};
  private:
    LyonPrimaryGeneratorAction* lyonPrimaryGeneratorAction;
    static LyonRunManager* lyonRunMan;
};

#endif
