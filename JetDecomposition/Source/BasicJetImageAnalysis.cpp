#include <fstream>
#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLorentzVector.h"

#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

#include "CommandLine.h"
#include "CustomAssert.h"
#include "Constants.h"

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName = CL.Get("input");
   string InputTreeName = CL.Get("tree", "DiscretizedJetTree");
   string OutputBase    = CL.Get("output", "JetImages");
   double MinJetPT      = CL.GetDouble("minjetpt", 300);
   int MaxEntry         = CL.GetInt("maxentry", -1);

   TFile InputFile(InputFileName.c_str());
   TTree *Tree = (TTree *)InputFile.Get(InputTreeName.c_str());
   Assert(Tree != NULL, "Error getting tree for jet images");
   
   ofstream OutputFile(OutputBase + ".dat");

   int HydjetEntry, PythiaEntry, AllJetCount, JetCount, JetIndex;
   float JetPT, JetEta, JetPhi;
   float Images[JetImageBinCount][JetImageBinCount][7];

   Tree->SetBranchAddress("HydjetEntry", &HydjetEntry);
   Tree->SetBranchAddress("PythiaEntry", &PythiaEntry);
   Tree->SetBranchAddress("AllJetCount", &AllJetCount);
   Tree->SetBranchAddress("JetCount", &JetCount);
   Tree->SetBranchAddress("JetIndex", &JetIndex);
   Tree->SetBranchAddress("JetPT", &JetPT);
   Tree->SetBranchAddress("JetEta", &JetEta);
   Tree->SetBranchAddress("JetPhi", &JetPhi);
   Tree->SetBranchAddress("Images", &Images);

   int Count = 0;
   TH2D HAverageCETSignal("HAverageCETSignal", "Signal charged ET;Eta Bin;Phi Bin",
      JetImageBinCount, 0, JetImageBinCount, JetImageBinCount, 0, JetImageBinCount);
   TH2D HAverageNETSignal("HAverageNETSignal", "Signal neutral ET;Eta Bin;Phi Bin",
      JetImageBinCount, 0, JetImageBinCount, JetImageBinCount, 0, JetImageBinCount);
   TH2D HAverageCETBackground("HAverageCETBackground", "Background charged ET;Eta Bin;Phi Bin",
      JetImageBinCount, 0, JetImageBinCount, JetImageBinCount, 0, JetImageBinCount);
   TH2D HAverageNETBackground("HAverageNETBackground", "Background neutral ET;Eta Bin;Phi Bin",
      JetImageBinCount, 0, JetImageBinCount, JetImageBinCount, 0, JetImageBinCount);

   HAverageCETSignal.SetStats(0);
   HAverageNETSignal.SetStats(0);
   HAverageCETBackground.SetStats(0);
   HAverageNETBackground.SetStats(0);

   double R0 = 0.4;

   int EntryCount = Tree->GetEntries();
   if(MaxEntry > 0)
      EntryCount = min(MaxEntry, EntryCount);
   for(int iE = 0; iE < EntryCount; iE++)
   {
      Tree->GetEntry(iE);

      if(JetPT < MinJetPT)
         continue;

      Count = Count + 1;

      vector<PseudoJet> Particles;

      for(int i = 0; i < JetImageBinCount; i++)
      {
         for(int j = 0; j < JetImageBinCount; j++)
         {
            HAverageCETSignal.Fill(i + 0.5, j + 0.5, Images[i][j][1]);
            HAverageNETSignal.Fill(i + 0.5, j + 0.5, Images[i][j][2]);
            HAverageCETBackground.Fill(i + 0.5, j + 0.5, Images[i][j][4]);
            HAverageNETBackground.Fill(i + 0.5, j + 0.5, Images[i][j][5]);

            double Eta = (i + 0.5 - 0.5 * JetImageBinCount) * JetImageBinSize;
            double Phi = (j + 0.5 - 0.5 * JetImageBinCount) * JetImageBinSize;
            
            double ETSignal = Images[i][j][1] + Images[i][j][2];
            double ETBackground = Images[i][j][4] + Images[i][j][5];

            if(ETSignal > 0.01)
            {
               PseudoJet Particle;
               Particle.reset_PtYPhiM(ETSignal, Eta, Phi, 0.0);
               Particle.set_user_index(10101);
               Particles.push_back(Particle);
            }
            if(ETBackground > 0.01)
            {
               PseudoJet Particle;
               Particle.reset_PtYPhiM(ETBackground, Eta, Phi, 0.0);
               Particle.set_user_index(10102);
               Particles.push_back(Particle);
            }
         }
      }

      JetDefinition Definition(antikt_algorithm, R0);
      ClusterSequence Sequence(Particles, Definition);
      vector<PseudoJet> Jets = sorted_by_pt(Sequence.inclusive_jets());

      PseudoJet SJet(0, 0, 0, 0);
      PseudoJet BJet(0, 0, 0, 0);
      PseudoJet Jet = Jets[0];
      
      double JetP = sqrt(Jet.modp2());
      PseudoJet JetAxis(Jet.px() / JetP, Jet.py() / JetP, Jet.pz() / JetP, 1);

      for(int iPart = 0; iPart < Jet.constituents().size(); iPart++)
      {
         if(Jet.constituents()[iPart].user_index() == 10101)
            SJet = SJet + Jet.constituents()[iPart];
         if(Jet.constituents()[iPart].user_index() == 10102)
            BJet = BJet + Jet.constituents()[iPart];
      }

      double BPlus   = BJet.e() - BJet.px() * JetAxis.px() - BJet.py() * JetAxis.py() - BJet.pz() * JetAxis.pz();
      double BMinus  = BJet.e() + BJet.px() * JetAxis.px() + BJet.py() * JetAxis.py() + BJet.pz() * JetAxis.pz();
      double SPlus   = SJet.e() - SJet.px() * JetAxis.px() - SJet.py() * JetAxis.py() - SJet.pz() * JetAxis.pz();
      double SMinus  = SJet.e() + SJet.px() * JetAxis.px() + SJet.py() * JetAxis.py() + SJet.pz() * JetAxis.pz();
      double ReconstructedMass = sqrt(SJet.m() * SJet.m() + 2 * SJet.e() * BPlus + BJet.m() * BJet.m());
      
      // double Product = BJet.e() * SJet.e() - BJet.px() * SJet.px() - BJet.py() * SJet.py() - BJet.pz() * SJet.pz();
      // double ReconstructedMass = sqrt(SJet.m() * SJet.m() + 2 * Product + BJet.m() * BJet.m());

      OutputFile << Jet.pt() << " " << Jet.m()
         << " " << SJet.pt() << " " << SJet.m()
         << " " << BJet.pt() << " " << BJet.m()
         << " " << ReconstructedMass << endl;
   }

   HAverageCETSignal.Scale(double(1) / Count);
   HAverageNETSignal.Scale(double(1) / Count);
   HAverageCETBackground.Scale(double(1) / Count);
   HAverageNETBackground.Scale(double(1) / Count);

   TCanvas Canvas;

   Canvas.Divide(2, 2);

   Canvas.cd(1);
   Canvas.cd(1)->SetGridx();
   Canvas.cd(1)->SetGridy();
   HAverageCETSignal.Draw("colz");
   Canvas.cd(2);
   Canvas.cd(2)->SetGridx();
   Canvas.cd(2)->SetGridy();
   HAverageNETSignal.Draw("colz");
   Canvas.cd(3);
   Canvas.cd(3)->SetGridx();
   Canvas.cd(3)->SetGridy();
   HAverageCETBackground.Draw("colz");
   Canvas.cd(4);
   Canvas.cd(4)->SetGridx();
   Canvas.cd(4)->SetGridy();
   HAverageNETBackground.Draw("colz");

   Canvas.SaveAs((OutputBase + ".pdf").c_str());

   OutputFile.close();
   InputFile.Close();

   return 0;
}





