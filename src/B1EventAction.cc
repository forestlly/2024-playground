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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"

#include <G4UnitsTable.hh>

#include "B1RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.; // energy deposition initialization
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event*)
{   
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);

  // calculate momentum
  for(G4int i=0; i<6; i++) {
    G4ThreeVector fspacing_three_vector = hit_positions.at(i+1) - hit_positions.at(i); // relative position vector
    G4double spacing_three_vector = fspacing_three_vector.mag(); // vector length
    G4double time_interval = hit_times.at(i+1) - hit_times.at(i); // time interval
    Record_Momentum(spacing_three_vector / time_interval);
  }

  // print the momentum at the end of each event
  G4cout
     // << " Cumulated dose per run, in scoring volume : "
     // // << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------Momentum before hitting on the------------------------------"
     << G4endl
     << G4endl;

  for(G4int j=0; j<6; j++) {
    G4double momentum = hit_momenta.at(j);
    G4cout
      << G4endl
      << j+2 << "-th dectector is: " << momentum << " MeV" // G4BestUnit(momentum,"Momentum")
      << G4endl;
  }

  G4cout
     // << " Cumulated dose per run, in scoring volume : "
     // // << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;

  // output momentum information into a root file
  // create a new logout file
  std::unique_ptr<TFile> Hit_Momenta_out( TFile::Open("Hit_Momenta.root","RECREATE") );

  G4double p;
  // create a TTree instance to load momenta information
  TTree* momenta_tree = new TTree("momenta_tree","Proton Hitting Momenta on the Last 6 Dectectors");
  momenta_tree->Branch("hit_momenta", &p);
  // fill entries into the tree
  for(G4int k=0; k<6; k++) {
    p = hit_momenta.at(k);
    momenta_tree->Fill();
  }
  // go into the logout file and write information
  Hit_Momenta_out->cd();
  momenta_tree->Write();
  momenta_tree->Print();
  Hit_Momenta_out->Close();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
