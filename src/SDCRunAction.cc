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
// $Id: SDCRunAction.cc 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file SDCRunAction.cc
/// \brief Implementation of the SDCRunAction class

#include "SDCRunAction.hh"
#include "SDCPrimaryGeneratorAction.hh"
#include "SDCRun.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "SDCAnalysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDCRunAction::SDCRunAction()
 : G4UserRunAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDCRunAction::~SDCRunAction()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* SDCRunAction::GenerateRun()
{ return new SDCRun; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDCRunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  

  // Create analysis manager
  // Notice: it must be done the same way in master and workers
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstNtupleId(1);
  
  /*
  //Create a one-column ntuple
  analysisManager->CreateNtuple("SDC", "Energy and doses");
  // 1) total energy released in the crystals (double), MeV
  analysisManager->CreateNtupleDColumn("Energy"); 
  //ok, done
  analysisManager->FinishNtuple();
  */

  analysisManager->CreateH1("h1","Energy",3000,0.,3000.);

  // Create a new output file
  analysisManager->OpenFile("prova.txt");

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDCRunAction::EndOfRunAction(const G4Run* run)
{
  //retrieve the number of events produced in the run
  G4int nofEvents = run->GetNumberOfEvent();

  //do nothing, if no events were processed
  if (nofEvents == 0) return;
  
  // Run conditions
  // This retrieves the UserPrimaryGeneratorAction object: it is retrieved through the 
  // G4RunManager. 
  //
  // Following the SDCActionInitialization, the UserPrimaryGeneratorAction 
  // exists for all workers (-> Build()) but not for the master (-> BuildForMaster()). The 
  // SDCRunAction instead exists for the master and for the worker. So, when the function is 
  // executed by the master, no pointer is found for the primary generator and 
  // generatorAction = NULL

  const SDCPrimaryGeneratorAction* generatorAction = static_cast<const SDCPrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
 

  G4String partName;
  if (generatorAction) 
  {
      // The GetParticleGun() is defined inside the SDCPrimaryGeneratorAction.hh file and
      // it returns a pointer to the current concrete implementation of the G4VPrimaryGeneratorAction
      // Note that GetParticleDefinition return a 'const' type in the case of the use of the G4ParticleGun but not
      // in the case of the G4ParticleGeneralSource. Hence the GetParticleGun in the SDCPrimaryGeneratorAction.hh
      // must not be defined 'const' in the case of the GeneralParticleSource
      //
    G4ParticleDefinition* particle = generatorAction -> GetParticleGun() -> GetParticleDefinition();
    partName = particle->GetParticleName();
  }  
  
  //results
  //
  const SDCRun* theRun = static_cast<const SDCRun*>(run);
  G4int nbGoodEvents = theRun->GetNGoodEvents();
        
  //print
  //
  if (IsMaster())
  {
    G4cout
     << "\n--------------------End of Global Run-----------------------"
     << " \n The run was " << nofEvents << " events ";
  }
  else
  {
    G4cout
     << "\n--------------------End of Local Run------------------------"
     << " \n The run was " << nofEvents << " "<< partName;
  }      
  G4cout
     << "; Nb of events with energy deposit: " << nbGoodEvents
     << "\n------------------------------------------------------------\n"
     << G4endl;


  //save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  man->Write();
  man->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
