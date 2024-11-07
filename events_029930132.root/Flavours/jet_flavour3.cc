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

std::vector<fastjet::PseudoJet> applyFlavouredAntiKtAlgorithm(std::vector<fastjet::PseudoJet>& particles, double R) {

  fastjet::contrib::FlavRecombiner flav_recombiner; 
  
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R); 
  jet_def.set_recombiner(&flav_recombiner);

  std::vector<fastjet::PseudoJet> flav_jets = jet_def(particles);

  return sorted_by_pt(flav_jets);
}

bool ContainsSandSbar(const std::string &entry) {
    return (entry.find(" s ") != std::string::npos) || (entry.find(" sbar ") != std::string::npos) ||
           (entry.find("sbar]") != std::string::npos) || (entry.find("[s ") != std::string::npos) || (entry.find("s]") != std::string::npos);
}

bool CompareColumns(const std::vector<int>& column1, const std::vector<int>& column2) {
    // Check if both columns have the same size
    if (column1.size() != column2.size()) {
        return false; // Sizes don't match, so columns cannot match
    }

    // Compare corresponding elements
    for (size_t i = 0; i < column1.size(); ++i) {
        if (column1[i] != column2[i]) {
            return false; // Found a mismatch
        }
    }

    return true; // All elements match
}

void jet_flavour3()
{

  ROOT::RDataFrame df("events", "events_029930132.root");

  auto df1 = df.Define("Px", convertToDouble, {"Particle.momentum.x"});
  auto df2 = df1.Define("Py", convertToDouble, {"Particle.momentum.y"});
  auto df3 = df2.Define("Pz", convertToDouble, {"Particle.momentum.z"});
  auto df4 = df3.Define("M", "Particle.mass").Define("E", "sqrt(Px*Px + Py*Py + Pz*Pz + M*M)").Define("ID", "Particle.PDG");
  auto df5 = df4.Define("Status", "Particle.generatorStatus").Define("Stable_Particles", [](ROOT::RVec<int>& Status, ROOT::RVec<double>& Px, ROOT::RVec<double>& Py, ROOT::RVec<double>& Pz, ROOT::RVec<double>& E, ROOT::RVec<int>& ID){
								       ROOT::RVec<double> Stable_Px, Stable_Py, Stable_Pz, Stable_E; ROOT::RVec<int> Stable_ID;
			  for (size_t i = 0; i < Status.size(); i++){
			    if (Status[i] == 1) 
			      {Stable_Px.push_back(Px[i]);
				Stable_Py.push_back(Py[i]);
				Stable_Pz.push_back(Pz[i]);
				Stable_E.push_back(E[i]);
				Stable_ID.push_back(ID[i]);} } return std::make_tuple(Stable_Px, Stable_Py, Stable_Pz, Stable_E, Stable_ID);}, {"Status", "Px", "Py", "Pz", "E", "ID"}).Define("Stable_Px", "std::get<0>(Stable_Particles)").Define("Stable_Py", "std::get<1>(Stable_Particles)").Define("Stable_Pz", "std::get<2>(Stable_Particles)").Define("Stable_E", "std::get<3>(Stable_Particles)").Define("Stable_ID", "std::get<4>(Stable_Particles)");
  auto df6 = df5.Define("PJs", [](ROOT::RVec<double>& Px, 
                                  ROOT::RVec<double>& Py,
                                  ROOT::RVec<double>& Pz,
                                  ROOT::RVec<double>& E,
				  ROOT::RVec<int>& ID){
			  return createPseudoJets(Px, Py, Pz, E, ID);}, {"Stable_Px", "Stable_Py", "Stable_Pz", "Stable_E", "Stable_ID"});
  auto df7 = df6.Define("Flav_Jets", [](std::vector<fastjet::PseudoJet>& PJs) {
			    return applyFlavouredAntiKtAlgorithm(PJs, 0.8);}, {"PJs"}).Filter("Flav_Jets.size() > 1");
  auto df8 = df7.Define("SDFlav_Jets",[](std::vector<fastjet::PseudoJet>& Flav_Jets) {
			      std::vector<fastjet::PseudoJet> sdflav_jets; 
			      SDFlavourCalc sdFlavCalc;
                              for (int i = 0; i < Flav_Jets.size(); i++) {
				sdFlavCalc(Flav_Jets[i]);
				sdflav_jets.push_back(Flav_Jets[i]);};
			      return sdflav_jets;}, {"Flav_Jets"});
  auto df9 = df8.Define("Flav_Description", [](std::vector<fastjet::PseudoJet>& jets) {
			      std::vector<std::string> descriptions;
                              for (const auto& jet : jets) {
  		                descriptions.push_back(fastjet::contrib::FlavHistory::current_flavour_of(jet).description());}
			      return descriptions;}, {"Flav_Jets"});
  auto df10 = df9.Define("SDFlav_Description", [](std::vector<fastjet::PseudoJet>& jets) {
			     std::vector<std::string> descriptions;
                             for (const auto& jet : jets) {
  		               descriptions.push_back(fastjet::contrib::FlavHistory::current_flavour_of(jet).description());}
			     return descriptions;}, {"SDFlav_Jets"});
  auto df11 = df10.Define("N_s_Flav_Jets", [](std::vector<std::string>& Flav_Description) {
			     int totalCount = 0;
                             for (const auto& entry : Flav_Description) {
  		               // Check if "s" or "sbar" is present
        if (ContainsSandSbar(entry)) {
            totalCount += 1; // Add 1 if "s" or "sbar" is found
        } }
			     return totalCount;}, {"Flav_Description"}).Define("N_s_SDFlav_Jets", [](std::vector<std::string>& SDFlav_Description) {
			     int totalCount = 0;
                             for (const auto& entry : SDFlav_Description) {
  		               // Check if "s" or "sbar" is present
        if (ContainsSandSbar(entry)) {
            totalCount += 1; // Add 1 if "s" or "sbar" is found
        } }
			     return totalCount;}, {"SDFlav_Description"});
  auto df12 = df11.Filter("N_s_Flav_Jets > 1 && N_s_SDFlav_Jets > 1");
  auto df13 = df12.Define("Label_S_Flav_Jets", [](std::vector<std::string>& Flav_Description) {
			     std::vector<int> labels;
                             for (size_t i = 0; i < Flav_Description.size(); i++) {
  		               // Check if "s" or "sbar" is present
        if (ContainsSandSbar(Flav_Description[i])) {
	    labels.push_back(i);
        } }
			     return labels;}, {"Flav_Description"}).Define("Label_S_SDFlav_Jets", [](std::vector<std::string>& SDFlav_Description) {
			     std::vector<int> labels;
                             for (size_t i = 0; i < SDFlav_Description.size(); i++) {
  		               // Check if "s" or "sbar" is present
        if (ContainsSandSbar(SDFlav_Description[i])) {
	    labels.push_back(i);
        } }
			     return labels;}, {"SDFlav_Description"});
  auto df14 = df13.Filter("Label_S_Flav_Jets.size() == Label_S_SDFlav_Jets.size()");
  auto rdf = df14.Define("Label_Match", [](std::vector<int> Label_S_Flav_Jets, std::vector<int> Label_S_SDFlav_Jets) { return CompareColumns(Label_S_Flav_Jets, Label_S_SDFlav_Jets);}, {"Label_S_Flav_Jets", "Label_S_SDFlav_Jets"});


auto h = rdf.Histo1D({"Label", "Strange Jets Map", 2u, 0, 2}, "Label_Match");

h->GetXaxis()->SetTitle("Label association");

auto c = new TCanvas("c");
h->DrawCopy();

auto displayOutput_1 = rdf.Display("Label_S_Flav_Jets");

std::ofstream outFile_1("Label_S_Flav_Jets_0.8.txt");

if (outFile_1.is_open()) {

 std::streambuf* coutBuffer = std::cout.rdbuf();
 std::cout.rdbuf(outFile_1.rdbuf());

 displayOutput_1->Print();

 std::cout.rdbuf(coutBuffer);

 outFile_1.close();
  } 

  auto displayOutput_2 = rdf.Display("Label_S_SDFlav_Jets");

  std::ofstream outFile_2("Label_S_SDFlav_Jets_0.8.txt");

  if (outFile_2.is_open()) {

  std::streambuf* coutBuffer = std::cout.rdbuf();
  std::cout.rdbuf(outFile_2.rdbuf());

  displayOutput_2->Print();

  std::cout.rdbuf(coutBuffer);

  outFile_2.close();
  } 

}










 

 
