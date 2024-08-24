#include <ROOT/RVec.hxx>
#include <vector>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>

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
                                                 const ROOT::RVec<double>& e) {
    std::vector<fastjet::PseudoJet> jets;
 
    for (size_t i = 0; i < px.size(); ++i) {
        jets.emplace_back(px[i], py[i], pz[i], e[i]);
    }

    return jets;
}

std::vector<fastjet::PseudoJet> applyAntiKtAlgorithm(const std::vector<fastjet::PseudoJet>& particles, double R) {
   
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
    fastjet::ClusterSequence clust_seq(particles, jet_def);

    return sorted_by_pt(clust_seq.inclusive_jets());
}

void fastjet_fcc1()
{

  ROOT::RDataFrame df("events", "events_029930132.root");

  auto df0 = df.Define("Px", convertToDouble, {"ReconstructedParticles.momentum.x"});
  auto df1 = df0.Define("Py", convertToDouble, {"ReconstructedParticles.momentum.y"});
  auto df2 = df1.Define("Pz", convertToDouble, {"ReconstructedParticles.momentum.z"});
  auto df3 = df2.Define("E", convertToDouble, {"ReconstructedParticles.energy"});
  auto df4 = df3.Define("PJs", [](const ROOT::RVec<double>& Px, 
                                  const ROOT::RVec<double>& Py,
                                  const ROOT::RVec<double>& Pz,
                                  const ROOT::RVec<double>& E){
                                  return createPseudoJets(Px, Py, Pz, E);}, {"Px", "Py", "Pz", "E"});
  auto df5 = df4.Define("Jets", [](const std::vector<fastjet::PseudoJet>& PJs){
				   return applyAntiKtAlgorithm(PJs, 0.8);}, {"PJs"});
  auto df6 = df5.Filter("Jets.size() > 1");
  auto df7 = df6.Define("Dijets_m", [](const std::vector<fastjet::PseudoJet>& Jets){
			               return (Jets[0] + Jets[1]).m();}, {"Jets"});
  auto df8 = df7.Define("Dijets_eta", [](const std::vector<fastjet::PseudoJet>& Jets){
			                 return (Jets[0] + Jets[1]).eta();}, {"Jets"});
  auto rdf = df8.Define("Dijets_pt", [](const std::vector<fastjet::PseudoJet>& Jets){
			                return (Jets[0] + Jets[1]).pt();}, {"Jets"});
  
  auto h1 = rdf.Histo1D({"pT", "Dijet", 200u, 0., 200.}, "Dijets_pt");
  auto h2 = rdf.Histo1D({"M", "Dijet", 250u, 0., 250.}, "Dijets_m");
  auto h3 = rdf.Histo1D({"Eta", "Dijet", 20u, -10., 10.}, "Dijets_eta");

  h1->GetXaxis()->SetTitle("pT");
  h2->GetXaxis()->SetTitle("M");
  h3->GetXaxis()->SetTitle("Eta");

  auto c1 = new TCanvas("c1");
  h1->DrawCopy();

  auto c2 = new TCanvas("c2");
  h2->DrawCopy();

  auto c3 = new TCanvas("c3");
  h3->DrawCopy();

}
