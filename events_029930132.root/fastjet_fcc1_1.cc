#include <iostream>
#include <fstream>
#include <ROOT/RVec.hxx>
#include <vector>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/SDFlavPlugin.hh>
#include <fastjet/FlavInfo.hh>
#include <fastjet/tools/Filter.hh>
#include <fastjet/Selector.hh>
//#include <TRatioPlot.h>

ROOT::RVec<double> convertToDouble(const ROOT::RVec<float> &column) {

  ROOT::RVec<double> col_value_double(column.size());
    
  for (size_t i = 0; i < column.size(); ++i) {
      col_value_double[i] = static_cast<double>(column[i]);
  }
    
  return col_value_double;
}

std::vector<fastjet::PseudoJet> createPseudoJets(ROOT::RVec<double>& px, 
                                                 ROOT::RVec<double>& py, 
                                                 ROOT::RVec<double>& pz, 
                                                 ROOT::RVec<double>& e,
                                                 ROOT::RVec<int>& id) {
    std::vector<fastjet::PseudoJet> jets;
 
    for (size_t i = 0; i < px.size(); ++i) {

      fastjet::PseudoJet jet(px[i], py[i], pz[i], e[i]);
      jet.set_user_info(new fastjet::contrib::FlavHistory(id[i])); 
      
      jets.push_back(jet);
    }

    return jets;
}

std::vector<fastjet::PseudoJet> applyAntiKtAlgorithm(std::vector<fastjet::PseudoJet>& particles, double R) {
   
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
    // fastjet::ClusterSequence clust_seq(particles, jet_def);
    std::vector<fastjet::PseudoJet> jets = jet_def(particles);

    return sorted_by_pt(jets);
}

std::vector<fastjet::PseudoJet> applyFlavouredAntiKtAlgorithm(std::vector<fastjet::PseudoJet>& particles, double R) {

  fastjet::contrib::FlavRecombiner flav_recombiner; 
  

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R); 
  jet_def.set_recombiner(&flav_recombiner);

  //fastjet::ClusterSequence clust_seq(particles, jet_def);

  std::vector<fastjet::PseudoJet> flav_jets = jet_def(particles);

  //SDFlavourCalc sdFlavCalc;
  //sdFlavCalc(flav_jets[1]);

  //std::vector<fastjet::PseudoJet> sdflav_jets = jet_def(particles);
  //sdFlavCalc(sdflav_jets);

  //for (int i = 0; i < flav_jets.size(); i++) {
    //if (flav_jets[i].has_constituents()){
    //SDFlavourCalc sdFlavCalc;
    //sdFlavCalc(flav_jets[i]);
    //}
    //sdflav_jets.push_back(flav_jets[i]);
  //}

  return sorted_by_pt(flav_jets);
}

//std::vector<fastjet::PseudoJet> applyFilter(std::vector<fastjet::PseudoJet>& jets, double R, int N) {
//const Filter& f = Filter(JetDefinition(cambridge_algorithm, R), SelectorNHardest(N));
//std::vector<fastjet::PseudoJet> filtered_jets;
//for (int i; i < jets.size(); i++) {

//}
//}

void fastjet_fcc1_1()
{

  ROOT::RDataFrame df("events", "events_029930132.root");

  auto df0 = df.Define("Px", convertToDouble, {"Particle.momentum.x"});
  auto df1 = df0.Define("Py", convertToDouble, {"Particle.momentum.y"});
  auto df2 = df1.Define("Pz", convertToDouble, {"Particle.momentum.z"});
  auto df3 = df2.Define("M", "Particle.mass");
  auto df4 = df3.Define("E", "sqrt(Px*Px + Py*Py + Pz*Pz + M*M)");
  auto df5 = df4.Define("ID", "Particle.PDG");
  auto df6 = df5.Define("PJs", [](ROOT::RVec<double>& Px, 
                                  ROOT::RVec<double>& Py,
                                  ROOT::RVec<double>& Pz,
                                  ROOT::RVec<double>& E,
				  ROOT::RVec<int>& ID){
			  return createPseudoJets(Px, Py, Pz, E, ID);}, {"Px", "Py", "Pz", "E", "ID"});
  auto df7 = df6.Define("Normal_Jets", [](std::vector<fastjet::PseudoJet>& PJs){
				   return applyAntiKtAlgorithm(PJs, 0.8);}, {"PJs"});
  auto df8 = df7.Define("N_Normal_Jets", "Normal_Jets.size()");
  auto df9 = df8.Filter("Normal_Jets.size() > 1");
  auto df10 = df9.Define("Normal_Dijets_phi", [](std::vector<fastjet::PseudoJet>& Normal_Jets){
			               return (Normal_Jets[0] + Normal_Jets[1]).phi();}, {"Normal_Jets"});
  auto df11 = df10.Define("Normal_Dijets_Eta", [](std::vector<fastjet::PseudoJet>& Normal_Jets){
			               return (Normal_Jets[0] + Normal_Jets[1]).eta();}, {"Normal_Jets"});
  auto df12 = df11.Define("Normal_Dijets_Pt", [](std::vector<fastjet::PseudoJet>& Normal_Jets){
			               return (Normal_Jets[0] + Normal_Jets[1]).pt();}, {"Normal_Jets"});
 auto df13 = df12.Define("Flav_Jets", [](std::vector<fastjet::PseudoJet>& PJs) {
				   return applyFlavouredAntiKtAlgorithm(PJs, 0.8);}, {"PJs"});
  auto df14 = df13.Filter("Flav_Jets.size() > 1");
 auto df15 = df14.Define("SDFlav_Jets",[](std::vector<fastjet::PseudoJet>& Flav_Jets) {
			      std::vector<fastjet::PseudoJet> sdflav_jets; 
			      SDFlavourCalc sdFlavCalc;
                              for (int i = 0; i < 2; i++) {
				sdFlavCalc(Flav_Jets[i]);
				sdflav_jets.push_back(Flav_Jets[i]);};
			      return sdflav_jets;}, {"Flav_Jets"});


 //auto df15_1 = df15.Define("FSDFlav_Jets",[](std::vector<fastjet::PseudoJet>& SDFlav_Jets) {
 //		      std::vector<fastjet::PseudoJet> fil_sdflav_jets; 
 //		      fastjet::Filter f(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3))//;
 // for (int i = 0; i < 2; i++) {
 //			f(SDFlav_Jets[i]);
 //			fil_sdflav_jets.push_back(SDFlav_Jets[i]);};
 //		      return fil_sdflav_jets;}, {"SDFlav_Jets"});
 // auto df15_2 = df15_1.Define("FSDFlav_Dijets_pt", [](std::vector<fastjet::PseudoJet>& FSDFlav_Jets) {
 //		      return (FSDFlav_Jets[0] + FSDFlav_Jets[1]).pt();}, {"FSDFlav_Jets"});
 // auto RDF = df15_2.Define("pt_ratio", "FSDFlav_Dijets_pt/Normal_Dijets_Pt");

   //  auto df17 = df16.Define("SDFlav_Dijets_Eta", [](std::vector<fastjet::PseudoJet>& SDFlav_Jets) {
   //	    return (SDFlav_Jets[0] + SDFlav_Jets[1]).eta();}, {"SDFlav_Jets"});
   //auto df18 = df17.Define("SDFlav_Dijets_Pt", [](std::vector<fastjet::PseudoJet>& SDFlav_Jets) {		                
   //	     return (SDFlav_Jets[0] + SDFlav_Jets[1]).pt();}, {"SDFlav_Jets"});
   //auto df19 = df18.Define("Mass_ratio", "SDFlav_Dijets_M/Normal_Dijets_M");
   //auto df20 = df19.Define("pT_ratio", "SDFlav_Dijets_Pt/Normal_Dijets_Pt");
   //auto RDF1 = df20.Define("Eta_ratio", "SDFlav_Dijets_Eta/Normal_Dijets_Eta");

    auto df16 = df15.Define("SDFlav_Dijets_phi", [](std::vector<fastjet::PseudoJet>& SDFlav_Jets) {
			      return (SDFlav_Jets[0] + SDFlav_Jets[1]).phi();}, {"SDFlav_Jets"});
  auto df17 = df16.Define("SDFlav_Dijets_Eta", [](std::vector<fastjet::PseudoJet>& SDFlav_Jets) {
			    return (SDFlav_Jets[0] + SDFlav_Jets[1]).eta();}, {"SDFlav_Jets"});
  auto df18 = df17.Define("SDFlav_Dijets_Pt", [](std::vector<fastjet::PseudoJet>& SDFlav_Jets) {		                
			     return (SDFlav_Jets[0] + SDFlav_Jets[1]).pt();}, {"SDFlav_Jets"});
  auto RDF1 = df18.Define("Pt_diff", "Normal_Dijets_Pt - SDFlav_Dijets_Pt");
  

  //auto df19 = df18.Define("phi_ratio", "SDFlav_Dijets_phi/Normal_Dijets_phi");
  //auto df20 = df19.Define("pT_ratio", "SDFlav_Dijets_Pt/Normal_Dijets_Pt");
  //auto RDF1 = df20.Define("Eta_ratio", "SDFlav_Dijets_Eta/Normal_Dijets_Eta");


  //auto h0 = RDF1.Histo1D({"N_Jets", "Normal_Dijet", 30u, 0, 30}, "N_Normal_Jets");

  //TH1D *hist1 = new TH1D("hist1", "Histogram", 400, 0., 500.);
  //TH1D *hist2 = new TH1D("hist2", "Histogram", 400, 0., 500.);

  //auto hist3 = RDF1.Fill<double>(hist1, {"Normal_Dijets_Pt"});

  auto h0 = RDF1.Histo1D({"pT", "Normal Dijet", 400u, 0., 400.}, "Normal_Dijets_Pt");
  h0->GetXaxis()->SetTitle("pT");
  auto c0 = new TCanvas("c0");
  h0->DrawCopy();

  
  auto h1 = RDF1.Histo1D({"pT", "SDFlav Dijet", 400u, 0., 400.}, "SDFlav_Dijets_Pt");
  h1->GetXaxis()->SetTitle("pT");
  auto c1 = new TCanvas("c1");
  h1->DrawCopy();

  TH1D* hist0 = h0.GetPtr();
  TH1D* hist1 = h1.GetPtr();

  //print(type(hist0));

  TH1D* hist2 = (TH1D*)hist1->Clone("hist2");
  hist2->Divide(hist0);
  auto c2 = new TCanvas("c2");
  hist2->DrawCopy();


  //auto h2 = RDF1.Histo1D({"pT_diff", "Diff.", 50, 0., 1.}, "Pt_diff");
  //h2->GetXaxis()->SetTitle("pT_diff");
  //auto c2 = new TCanvas("c2");
  //h2->DrawCopy();

  //TCanvas *canvas = new TCanvas("canvas", "Canvas for TRatioPlot", 800, 600);

  //auto h0_ptr = h0.GetPtr();
  //auto h1_ptr = h1.GetPtr();

  























  //auto h1 = RDF1.Histo1D({"pT_ratio", "Dijet_comparisons", 20u, 0., 20.}, "pT_ratio");
  //auto h2 = RDF1.Histo1D({"phi_ratio", "Dijet_comparisons", 20u, 0., 20.}, "phi_ratio");
  //auto h3 = RDF1.Histo1D({"Eta_ratio", "Dijet_comparisons", 20u, -10., 10.}, "Eta_ratio");
  //auto h4 = RDF1.Histo1D({"pT", "Normal_Dijet", 200u, 0., 200.}, "Normal_Dijets_Pt");

  //auto h5 = RDF.Histo1D({"pT_ratio", "Dijet_sub_comparisons", 20u, 0., 20.}, "pt_ratio");
 
  //h2->GetXaxis()->SetTitle("pT");

  // h2->GetXaxis()->SetTitle("M_ratio");
  //h3->GetXaxis()->SetTitle("Eta_ratio");

  //h4->GetXaxis()->SetTitle("normal_pt");
  //h5->GetXaxis()->SetTitle("pt_ratio");
  
  // h2->Draw();

  //auto c2 = new TCanvas("c2");
  //h2->DrawCopy();

  //auto c3 = new TCanvas("c3");
  //h3->DrawCopy();

  //auto c4 = new TCanvas("c4");
  //h4->DrawCopy();

  //auto c5 = new TCanvas("c5");
  //h5->DrawCopy();

 
  //auto rdf8 = rdf7.Define("Flav_Dijets_M", [](std::vector<fastjet::PseudoJet>& Flav_Jets) {
			    //		    return (Flav_Jets[0] + Flav_Jets[1]).m();}, {"Flav_Jets"});
  //auto rdf9 = rdf8.Define("Flav_Dijets_Eta", [](std::vector<fastjet::PseudoJet>& Flav_Jets) {
			    //		     return (Flav_Jets[0] + Flav_Jets[1]).eta();}, {"Flav_Jets"});
  //auto rdf10 = rdf9.Define("Flav_Dijets_Pt", [](std::vector<fastjet::PseudoJet>& Flav_Jets) {		                
			     //return (Flav_Jets[0] + Flav_Jets[1]).pt();}, {"Flav_Jets"});

 
//auto rdf7 = rdf6.Define("N_SDFlav_Jets", "SDFlav_Jets.size()");

  //auto rdf12 = rdf11.Define("Flav_Description", [](std::vector<fastjet::PseudoJet>& jets) {
  //		     std::vector<std::string> descriptions;
  //                         for (int i = 0; i < 2; i++) {
  //		       descriptions.push_back(fastjet::contrib::FlavHistory::current_flavour_of(jets[i]).description());}
  //		     return descriptions;}, {"Flav_Jets"});
  //auto RDF2 = rdf12.Define("SDFlav_Description", [](std::vector<fastjet::PseudoJet>& jets) {
  //		     std::vector<std::string> descriptions;
  //                         for (int i = 0; i < 2; i++) {
  //		       descriptions.push_back(fastjet::contrib::FlavHistory::current_flavour_of(jets[i]).description());}
  //		     return descriptions;}, {"SDFlav_Jets"});

  //auto displayOutput_1 = RDF2.Display("Flav_Description");

  //std::ofstream outFile_1("Flav_Descriptions.txt");

  //if (outFile_1.is_open()) {

  //std::streambuf* coutBuffer = std::cout.rdbuf();
  //std::cout.rdbuf(outFile_1.rdbuf());

  //displayOutput_1->Print();

  //  std::cout.rdbuf(coutBuffer);

  //outFile_1.close();
  //} 

  //auto displayOutput_2 = RDF2.Display("SDFlav_Description");

  //std::ofstream outFile_2("SDFlav_Descriptions.txt");

  //if (outFile_2.is_open()) {

  //std::streambuf* coutBuffer = std::cout.rdbuf();
  //std::cout.rdbuf(outFile_2.rdbuf());

  //displayOutput_2->Print();

  //std::cout.rdbuf(coutBuffer);

  //outFile_2.close();
  //} 

  //auto h4 = RDF2.Histo1D({"N", "SDFlav_Dijet", 30u, 0, 30}, "N_SDFlav_Jets");

  //auto h5 = RDF2.Histo1D({"pT", "SDFlav_Dijet", 20u, 0., 20.}, "Flav_Dijets_Pt");
  //auto h6 = RDF2.Histo1D({"M", "SDFlav_Dijet", 20u, 0., 20.}, "Flav_Dijets_M");
  //auto h7 = RDF2.Histo1D({"Eta", "SDFlav_Dijet", 20u, -10., 10.}, "Flav_Dijets_Eta");

  // h4->GetXaxis()->SetTitle("N");

  //h5->GetXaxis()->SetTitle("pT[GeV]");
  //h6->GetXaxis()->SetTitle("M[GeV]");
  //h7->GetXaxis()->SetTitle("Eta");

  //auto c4 = new TCanvas("c4");
  //h4->DrawCopy();  

  //auto c5 = new TCanvas("c5");
  //h5->DrawCopy();

  //auto c6 = new TCanvas("c6");
  //h6->DrawCopy();

  //auto c7 = new TCanvas("c7");
  //h7->DrawCopy();


}










 

 
