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

  std::vector<fastjet::PseudoJet> jets = jet_def(particles);

  return sorted_by_pt(jets);
}

std::vector<fastjet::PseudoJet> applyFlavouredAntiKtAlgorithm(std::vector<fastjet::PseudoJet>& particles, double R) {

  fastjet::contrib::FlavRecombiner flav_recombiner; 
  
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R); 
  jet_def.set_recombiner(&flav_recombiner);

  std::vector<fastjet::PseudoJet> flav_jets = jet_def(particles);

  return sorted_by_pt(flav_jets);
}

std::string GetFirstElement(const std::string &entry) {
    auto start = entry.find_first_of('[') + 1;
    auto end = entry.find_first_of(' ', start);
    if (end == std::string::npos) {
        end = entry.find_first_of(']', start);
    }
    return entry.substr(start, end - start);
}

void jet_flavour1()
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
  auto df8 = df7.Filter("Normal_Jets.size() > 1");
  auto df9 = df8.Define("Normal_Dijets_phi", [](std::vector<fastjet::PseudoJet>& Normal_Jets){
			               return (Normal_Jets[0] + Normal_Jets[1]).phi();}, {"Normal_Jets"});
  auto df10 = df9.Define("Normal_Dijets_Eta", [](std::vector<fastjet::PseudoJet>& Normal_Jets){
			               return (Normal_Jets[0] + Normal_Jets[1]).eta();}, {"Normal_Jets"});
  auto df11 = df10.Define("Normal_Dijets_Pt", [](std::vector<fastjet::PseudoJet>& Normal_Jets){
			               return (Normal_Jets[0] + Normal_Jets[1]).pt();}, {"Normal_Jets"});
  auto df12 = df11.Define("Flav_Jets", [](std::vector<fastjet::PseudoJet>& PJs) {
			   return applyFlavouredAntiKtAlgorithm(PJs, 0.8);}, {"PJs"});
  auto df13 = df12.Define("SDFlav_Jets",[](std::vector<fastjet::PseudoJet>& Flav_Jets) {
			      std::vector<fastjet::PseudoJet> sdflav_jets; 
			      SDFlavourCalc sdFlavCalc;
                              for (int i = 0; i < 2; i++) {
				sdFlavCalc(Flav_Jets[i]);
				sdflav_jets.push_back(Flav_Jets[i]);};
			      return sdflav_jets;}, {"Flav_Jets"});
  auto df14 = df13.Define("Flav_Description", [](std::vector<fastjet::PseudoJet>& jets) {
			      std::vector<std::string> descriptions;
                              for (const auto& jet : jets) {
  		                descriptions.push_back(fastjet::contrib::FlavHistory::current_flavour_of(jet).description());}
			      return descriptions;}, {"Flav_Jets"});
  auto df15 = df14.Define("SDFlav_Description", [](std::vector<fastjet::PseudoJet>& jets) {
			     std::vector<std::string> descriptions;
                             for (const auto& jet : jets) {
  		               descriptions.push_back(fastjet::contrib::FlavHistory::current_flavour_of(jet).description());}
			     return descriptions;}, {"SDFlav_Jets"});
  auto df16 = df15.Define("SDFJ", [](std::vector<std::string>& Flav_Description, std::vector<std::string>& SDFlav_Description, std::vector<fastjet::PseudoJet>& Flav_Jets) {  
      std::vector<fastjet::PseudoJet> jet;
      std::unordered_set<int> matchedIndices;
      for (size_t i = 0; i < SDFlav_Description.size(); i++) {
        std::string targetElement = GetFirstElement(SDFlav_Description[i]);
	int matchedEntry = -1;
        for (size_t j = 0; j < Flav_Description.size(); j++) {
            if (matchedIndices.count(j) == 0 && GetFirstElement(Flav_Description[j]) == targetElement) {
                matchedEntry = j; // Store the matched entry number
		matchedIndices.insert(j);
                break; // Stop searching after finding the first match
            }
        }
	jet.push_back(Flav_Jets[matchedEntry]);
      }; return jet;}, {"Flav_Description", "SDFlav_Description", "Flav_Jets"});
auto df17 = df16.Define("SDFlav_Dijets_Pt", [](const std::vector<fastjet::PseudoJet>& SDFlav_Jets){
			  return (SDFlav_Jets[0] + SDFlav_Jets[1]).pt();}, {"SDFlav_Jets"});
auto df18 = df17.Define("SDF_Dijets_Pt", [](const std::vector<fastjet::PseudoJet>& SDFJ){
			                return (SDFJ[0] + SDFJ[1]).pt();}, {"SDFJ"});
auto df19 = df18.Define("SDFlav_Dijets_Phi", [](const std::vector<fastjet::PseudoJet>& SDFlav_Jets){
			  return (SDFlav_Jets[0] + SDFlav_Jets[1]).phi();}, {"SDFlav_Jets"});
auto df20 = df19.Define("SDF_Dijets_Phi", [](const std::vector<fastjet::PseudoJet>& SDFJ){
			                return (SDFJ[0] + SDFJ[1]).phi();}, {"SDFJ"});
auto df21 = df20.Define("SDFlav_Dijets_Eta", [](const std::vector<fastjet::PseudoJet>& SDFlav_Jets){
			  return (SDFlav_Jets[0] + SDFlav_Jets[1]).eta();}, {"SDFlav_Jets"});
auto df22 = df21.Define("SDF_Dijets_Eta", [](const std::vector<fastjet::PseudoJet>& SDFJ){
			                return (SDFJ[0] + SDFJ[1]).eta();}, {"SDFJ"});
auto df23 = df22.Define("SDFlav_Dijets_M", [](const std::vector<fastjet::PseudoJet>& SDFlav_Jets){
			  return (SDFlav_Jets[0] + SDFlav_Jets[1]).m();}, {"SDFlav_Jets"});
auto rdf = df23.Define("SDF_Dijets_M", [](const std::vector<fastjet::PseudoJet>& SDFJ){
			                return (SDFJ[0] + SDFJ[1]).m();}, {"SDFJ"});



auto h1 = rdf.Histo1D({"pT", "SDFlav_Dijet", 400u, 0., 400.}, "SDFlav_Dijets_Pt");
auto h2 = rdf.Histo1D({"pT", "SDF_Dijet", 400u, 0., 400.}, "SDF_Dijets_Pt");
auto h3 = rdf.Histo1D({"phi", "SDFlav_Dijet", 10u, 0., 10.}, "SDFlav_Dijets_Phi");
auto h4 = rdf.Histo1D({"phi", "SDF_Dijet", 10u, 0., 10.}, "SDF_Dijets_Phi");
auto h5 = rdf.Histo1D({"eta", "SDFlav_Dijet", 10u, 0., 10.}, "SDFlav_Dijets_Eta");
auto h6 = rdf.Histo1D({"eta", "SDF_Dijet", 10u, 0., 10.}, "SDF_Dijets_Eta");
auto h7 = rdf.Histo1D({"M", "SDFlav_Dijet", 300u, 0., 300.}, "SDFlav_Dijets_M");
auto h8 = rdf.Histo1D({"M", "SDF_Dijet", 300u, 0., 300.}, "SDF_Dijets_M");

h1->GetXaxis()->SetTitle("pT [GeV]");
h2->GetXaxis()->SetTitle("pT [GeV]");
h3->GetXaxis()->SetTitle("#phi [radians]");
h4->GetXaxis()->SetTitle("#phi [radians]");
h5->GetXaxis()->SetTitle("#eta");
h6->GetXaxis()->SetTitle("#eta");
h7->GetXaxis()->SetTitle("M [GeV]");
h8->GetXaxis()->SetTitle("M [GeV]");

auto c1 = new TCanvas("c1");
h1->DrawCopy();  

auto c2 = new TCanvas("c2");
h2->DrawCopy();

auto c3 = new TCanvas("c3");
h3->DrawCopy();

auto c4 = new TCanvas("c4");
h4->DrawCopy();

auto c5 = new TCanvas("c5");
h5->DrawCopy();

auto c6 = new TCanvas("c6");
h6->DrawCopy();

auto c7 = new TCanvas("c7");
h7->DrawCopy();

auto c8 = new TCanvas("c8");
h8->DrawCopy();

TH1D* hist1 = h1.GetPtr();
TH1D* hist2 = h2.GetPtr();

  TH1D* H1 = (TH1D*)hist2->Clone("H1");
  H1->Divide(hist1);
  H1->SetTitle("pT ratio [R = 0.8]");
  H1->GetYaxis()->SetTitle("Ratio");
  auto C1 = new TCanvas("C1");
  H1->DrawCopy();

TH1D* hist3 = h3.GetPtr();
TH1D* hist4 = h4.GetPtr();

  TH1D* H2 = (TH1D*)hist4->Clone("H2");
  H2->Divide(hist3);
  H2->SetTitle("#phi ratio [R = 0.8]");
  H2->GetYaxis()->SetTitle("Ratio");
  auto C2 = new TCanvas("C2");
  H2->DrawCopy();

TH1D* hist5 = h5.GetPtr();
TH1D* hist6 = h6.GetPtr();

  TH1D* H3 = (TH1D*)hist6->Clone("H3");
  H3->Divide(hist5);
  H3->SetTitle("#eta ratio [R = 0.8]");
  H3->GetYaxis()->SetTitle("Ratio");
  auto C3 = new TCanvas("C3");
  H3->DrawCopy();

TH1D* hist7 = h7.GetPtr();
TH1D* hist8 = h8.GetPtr();

  TH1D* H4 = (TH1D*)hist8->Clone("H4");
  H4->Divide(hist7);
  H4->SetTitle("M ratio [R = 0.8]");
  H4->GetYaxis()->SetTitle("Ratio");
  auto C4 = new TCanvas("C4");
  H4->DrawCopy();


  auto displayOutput_1 = rdf.Display("Flav_Description");

  std::ofstream outFile_1("Flavour_0.8.txt");

  if (outFile_1.is_open()) {

  std::streambuf* coutBuffer = std::cout.rdbuf();
  std::cout.rdbuf(outFile_1.rdbuf());

  displayOutput_1->Print();

  std::cout.rdbuf(coutBuffer);

  outFile_1.close();
  } 

  auto displayOutput_2 = rdf.Display("SDFlav_Description");

  std::ofstream outFile_2("SDFlavour_0.8.txt");

  if (outFile_2.is_open()) {

  std::streambuf* coutBuffer = std::cout.rdbuf();
  std::cout.rdbuf(outFile_2.rdbuf());

  displayOutput_2->Print();

  std::cout.rdbuf(coutBuffer);

  outFile_2.close();
  } 

}










 

 
