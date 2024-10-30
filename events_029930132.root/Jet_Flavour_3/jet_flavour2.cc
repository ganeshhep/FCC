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

std::string GetLastElement(const std::string &entry) {
    auto end = entry.find_last_of(']');
    auto lastSpace = entry.find_last_of(' ', end - 1);

    // Check if there's no space, meaning it's a single-element list
    if (lastSpace == std::string::npos) {
        lastSpace = entry.find_first_of('[') + 1;
    } else {
        lastSpace += 1; // Move to the beginning of the last element
    }

    return entry.substr(lastSpace, end - lastSpace);
}

void jet_flavour2()
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
  auto rdf = df15.Define("Matched_Flav", [](std::vector<std::string>& Flav_Description, std::vector<std::string>& SDFlav_Description) {  
      std::vector<std::string> flav_descriptions;
      std::unordered_set<int> matchedIndices;
      for (size_t i = 0; i < SDFlav_Description.size(); i++) {
	std::string targetEntry = SDFlav_Description[i];
        std::string targetElement = GetLastElement(SDFlav_Description[i]);
        int matchedEntry = -1; // Initialize with -1 to indicate no match by default
	for (size_t j = 0; j < Flav_Description.size(); j++){
	  if (matchedIndices.count(j) == 0 && Flav_Description[j] == targetEntry){
	    matchedEntry = j;
	    matchedIndices.insert(j);
	    break;
	  }
	}
	if (matchedEntry == -1) {
        for (size_t j = 0; j < Flav_Description.size(); j++) {
            if (matchedIndices.count(j) == 0 && GetLastElement(Flav_Description[j]) == targetElement) {
                matchedEntry = j; // Store the matched entry number
		matchedIndices.insert(j);
                break; // Stop searching after finding the first match
            }
        }
	}
	flav_descriptions.push_back(Flav_Description[matchedEntry]);
      }; return flav_descriptions;}, {"Flav_Description", "SDFlav_Description"});

auto displayOutput_3 = rdf.Display("Matched_Flav");

  std::ofstream outFile_3("Matched_Flav_0.8.txt");

  if (outFile_3.is_open()) {

  std::streambuf* coutBuffer = std::cout.rdbuf();
  std::cout.rdbuf(outFile_3.rdbuf());

  displayOutput_3->Print();

  std::cout.rdbuf(coutBuffer);

  outFile_3.close();
  }

}










 

 
