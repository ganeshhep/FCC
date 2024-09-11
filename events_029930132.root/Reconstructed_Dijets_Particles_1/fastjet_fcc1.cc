#include <iostream>
#include <fstream>
#include <ROOT/RVec.hxx>
#include <vector>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/SDFlavPlugin.hh>
#include <fastjet/FlavInfo.hh>

ROOT::RVec<double> convertToDouble(const ROOT::RVec<float> &column) {

  ROOT::RVec<double> col_value_double(column.size());
    
  for (size_t i = 0; i < column.size(); ++i) {
      col_value_double[i] = static_cast<double>(column[i]);
  }
    
  return col_value_double;
}

std::vector<fastjet::PseudoJet> createPseudoJets(const ROOT::RVec<double>& px, 
                                                 const ROOT::RVec<double>& py, 
                                                 const ROOT::RVec<double>& pz, 
                                                 const ROOT::RVec<double>& e,
                                                 const ROOT::RVec<int>& id) {
    std::vector<fastjet::PseudoJet> jets;
 
    for (size_t i = 0; i < px.size(); ++i) {

      fastjet::PseudoJet jet(px[i], py[i], pz[i], e[i]);
      jet.set_user_info(new fastjet::contrib::FlavHistory(id[i])); 
      
      jets.push_back(jet);
    }

    return jets;
}

std::vector<fastjet::PseudoJet> applyAntiKtAlgorithm(const std::vector<fastjet::PseudoJet>& particles, double R) {
   
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
    fastjet::ClusterSequence clust_seq(particles, jet_def);

    return sorted_by_pt(clust_seq.inclusive_jets());
}

std::vector<fastjet::PseudoJet> applyFlavouredAntiKtAlgorithm(const std::vector<fastjet::PseudoJet>& particles, double R) {

  fastjet::contrib::FlavRecombiner flav_recombiner; 

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R); 
  jet_def.set_recombiner(&flav_recombiner);

  fastjet::ClusterSequence clust_seq(particles, jet_def);

  return sorted_by_pt(clust_seq.inclusive_jets());
}

void fastjet_fcc1()
{

  ROOT::RDataFrame df("events", "events_029930132.root");

  auto df0 = df.Define("Px", convertToDouble, {"Particle.momentum.x"});
  auto df1 = df0.Define("Py", convertToDouble, {"Particle.momentum.y"});
  auto df2 = df1.Define("Pz", convertToDouble, {"Particle.momentum.z"});
  auto df3 = df2.Define("M", "Particle.mass");
  auto df4 = df3.Define("E", "sqrt(Px*Px + Py*Py + Pz*Pz + M*M)");
  auto df5 = df4.Define("ID", "Particle.PDG");
  auto df6 = df5.Define("PJs", [](const ROOT::RVec<double>& Px, 
                                  const ROOT::RVec<double>& Py,
                                  const ROOT::RVec<double>& Pz,
                                  const ROOT::RVec<double>& E,
				  const ROOT::RVec<int>& ID){
			  return createPseudoJets(Px, Py, Pz, E, ID);}, {"Px", "Py", "Pz", "E", "ID"});

  auto rdf1 = df6.Define("Normal_Jets", [](const std::vector<fastjet::PseudoJet>& PJs){
				   return applyAntiKtAlgorithm(PJs, 0.8);}, {"PJs"});
  auto rdf2 = rdf1.Define("N_Normal_Jets", "Normal_Jets.size()");
  auto rdf3 = rdf2.Filter("Normal_Jets.size() > 1");
  auto rdf4 = rdf3.Define("Normal_Dijets_M", [](const std::vector<fastjet::PseudoJet>& Normal_Jets){
			               return (Normal_Jets[0] + Normal_Jets[1]).m();}, {"Normal_Jets"});
  auto rdf5 = rdf4.Define("Normal_Dijets_Eta", [](const std::vector<fastjet::PseudoJet>& Normal_Jets){
			               return (Normal_Jets[0] + Normal_Jets[1]).eta();}, {"Normal_Jets"});
  auto RDF1 = rdf5.Define("Normal_Dijets_Pt", [](const std::vector<fastjet::PseudoJet>& Normal_Jets){
			               return (Normal_Jets[0] + Normal_Jets[1]).pt();}, {"Normal_Jets"});

  auto h0 = RDF1.Histo1D({"N_Jets", "Normal_Dijet", 30u, 0, 30}, "N_Normal_Jets");
  auto h1 = RDF1.Histo1D({"pT", "Normal_Dijet", 400u, 0., 400.}, "Normal_Dijets_Pt");
  auto h2 = RDF1.Histo1D({"M", "Normal_Dijet", 800u, 0., 800.}, "Normal_Dijets_M");
  auto h3 = RDF1.Histo1D({"Eta", "Normal_Dijet", 20u, -10., 10.}, "Normal_Dijets_Eta");

  h0->GetXaxis()->SetTitle("N");
  h1->GetXaxis()->SetTitle("pT[GeV]");
  h2->GetXaxis()->SetTitle("M[GeV]");
  h3->GetXaxis()->SetTitle("Eta");

  auto c0 = new TCanvas("c0");
  h0->DrawCopy();

  auto c1 = new TCanvas("c1");
  h1->DrawCopy();

  auto c2 = new TCanvas("c2");
  h2->DrawCopy();

  auto c3 = new TCanvas("c3");
  h3->DrawCopy();

  auto rdf6 = df6.Define("Flav_Jets", [](const std::vector<fastjet::PseudoJet>& PJs){
				   return applyFlavouredAntiKtAlgorithm(PJs, 0.8);}, {"PJs"});
  auto rdf7 = rdf6.Define("N_Flav_Jets", "Flav_Jets.size()");
  auto rdf8 = rdf7.Filter("Flav_Jets.size() > 1");
  auto rdf9 = rdf8.Define("Flav_Dijets_M", [](const std::vector<fastjet::PseudoJet>& Flav_Jets){
			               return (Flav_Jets[0] + Flav_Jets[1]).m();}, {"Flav_Jets"});
  auto rdf10 = rdf9.Define("Flav_Dijets_Eta", [](const std::vector<fastjet::PseudoJet>& Flav_Jets){
			                 return (Flav_Jets[0] + Flav_Jets[1]).eta();}, {"Flav_Jets"});
  auto rdf11 = rdf10.Define("Flav_Dijets_Pt", [](const std::vector<fastjet::PseudoJet>& Flav_Jets){
			                return (Flav_Jets[0] + Flav_Jets[1]).pt();}, {"Flav_Jets"});
  auto RDF2 = rdf11.Define("Flav_Description", [](const std::vector<fastjet::PseudoJet>& jets) {
			     std::vector<std::string> descriptions;
                             for (const auto& jet : jets) {
			       descriptions.push_back(fastjet::contrib::FlavHistory::current_flavour_of(jet).description());
                             }
			     return descriptions;
			   }, {"Flav_Jets"});

  auto displayOutput_1 = RDF2.Display("Flav_Description");

  std::ofstream outFile_1("Flav_description.txt");

  if (outFile_1.is_open()) {

    std::streambuf* coutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(outFile_1.rdbuf());

    displayOutput_1->Print();

    std::cout.rdbuf(coutBuffer);

    outFile_1.close();
  } 

  auto h4 = RDF2.Histo1D({"N", "Flav_Dijet", 30u, 0, 30}, "N_Flav_Jets");
  auto h5 = RDF2.Histo1D({"pT", "Flav_Dijet", 400u, 0., 400.}, "Flav_Dijets_Pt");
  auto h6 = RDF2.Histo1D({"M", "Flav_Dijet", 800u, 0., 800.}, "Flav_Dijets_M");
  auto h7 = RDF2.Histo1D({"Eta", "FlavDijet", 20u, -10., 10.}, "Flav_Dijets_Eta");

  h4->GetXaxis()->SetTitle("N");
  h5->GetXaxis()->SetTitle("pT[GeV]");
  h6->GetXaxis()->SetTitle("M[GeV]");
  h7->GetXaxis()->SetTitle("Eta");

  auto c4 = new TCanvas("c4");
  h4->DrawCopy();  

  auto c5 = new TCanvas("c5");
  h5->DrawCopy();

  auto c6 = new TCanvas("c6");
  h6->DrawCopy();

  auto c7 = new TCanvas("c7");
  h7->DrawCopy();

}
