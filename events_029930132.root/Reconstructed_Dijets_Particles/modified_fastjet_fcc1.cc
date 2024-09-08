
#include <ROOT/RVec.hxx>
#include <vector>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/FlavKtPlugin.hh>  // Added for flavored anti-kt algorithm
#include <iostream>

// Function to convert ROOT::RVec<float> to ROOT::RVec<double>
ROOT::RVec<double> convertToDouble(const ROOT::RVec<float> &column) {
    ROOT::RVec<double> col_value_double(column.size());
    for (size_t i = 0; i < column.size(); ++i) {
        col_value_double[i] = static_cast<double>(column[i]);
    }
    return col_value_double;
}

// Function to create PseudoJets from px, py, pz, and energy
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

// Function to apply the standard anti-kt algorithm
std::vector<fastjet::PseudoJet> applyAntiKtAlgorithm(const std::vector<fastjet::PseudoJet>& particles, double R) {
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
    fastjet::ClusterSequence cs(particles, jet_def);
    return fastjet::sorted_by_pt(cs.inclusive_jets());
}

// Function to apply flavored anti-kt algorithm
std::vector<fastjet::PseudoJet> applyFlavoredAntiKtAlgorithm(const std::vector<fastjet::PseudoJet>& particles, double R) {
    fastjet::contrib::FlavKtPlugin flav_plugin(R);  // Create the flavored anti-kt plugin
    fastjet::JetDefinition jet_def(&flav_plugin);   // Define the jet definition with the plugin
    fastjet::ClusterSequence cs(particles, jet_def);
    return fastjet::sorted_by_pt(cs.inclusive_jets());
}

int main() {
    // Example usage: reconstruct particles from some input data (px, py, pz, e)

    // Placeholder for actual particle data
    ROOT::RVec<double> px = {10.0, 20.0, 30.0};  // Example data
    ROOT::RVec<double> py = {5.0, 15.0, 25.0};
    ROOT::RVec<double> pz = {2.0, 12.0, 22.0};
    ROOT::RVec<double> e = {50.0, 60.0, 70.0};

    // Create PseudoJets
    std::vector<fastjet::PseudoJet> particles = createPseudoJets(px, py, pz, e);

    // Apply standard anti-kt algorithm
    double R = 0.4;
    std::vector<fastjet::PseudoJet> antiKtJets = applyAntiKtAlgorithm(particles, R);
    std::cout << "Anti-Kt jets: " << antiKtJets.size() << std::endl;

    // Apply flavored anti-kt algorithm
    std::vector<fastjet::PseudoJet> flavAntiKtJets = applyFlavoredAntiKtAlgorithm(particles, R);
    std::cout << "Flavored Anti-Kt jets: " << flavAntiKtJets.size() << std::endl;

    return 0;
}
