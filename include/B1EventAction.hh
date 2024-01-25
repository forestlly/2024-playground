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
/// \file B1EventAction.hh
/// \brief Definition of the B1EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include <G4ThreeVector.hh>
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

class B1RunAction;

/// Event action class
///

class B1EventAction : public G4UserEventAction
{
  public:
    B1EventAction(B1RunAction* runAction);
    virtual ~B1EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    // Used in scoring volume
    void AddEdep(G4double edep) { fEdep += edep; }
    void Record_Position(G4ThreeVector hit_position) { hit_positions.push_back(hit_position); }
    void Record_Time(G4double hit_time) { hit_times.push_back(hit_time); }
    void Record_Momentum(G4double fhit_velocity) {
      // from 3d-velocity to 3d-momentum
      G4double hit_velocity = fhit_velocity * pow(10,6) * CLHEP::m/CLHEP::s;
      G4double hit_momentum = 938.272 * CLHEP::MeV *
        hit_velocity / sqrt(sqr(CLHEP::c_light) - sqr(hit_velocity));

      // record momentum
      hit_momenta.push_back(hit_momentum);
    }

  private:
    B1RunAction* fRunAction;
    G4double     fEdep;
    std::vector<G4ThreeVector> hit_positions; // positions of particles hitting on dectectors
    std::vector<G4double> hit_times; // times of particles hitting on dectectors
    std::vector<G4double> hit_momenta; // from position and time to momentum
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
