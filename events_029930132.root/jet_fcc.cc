void jet_fcc()
{
  //TFile *f = TFile::Open("events_029930132.root");
  //TTree *tree = f->Get<TTree>("events/Jet"); 
  using namespace ROOT;
  ROOT::RDataFrame df("events", "events_029930132.root");
  //ROOT::RDataFrame df(*tree);


  // Defining a new column "Jet_pt" in the root dataframe 'df'.
  auto df0 = df.Filter("Jet.momentum.x.size() > 1");
  auto df1 = df0.Define("LeadingJets_px", "auto myvec(Jet.momentum.x); myvec.resize(2); return myvec;");
  auto df2 = df1.Define("LeadingJets_py", "auto myvec(Jet.momentum.y); myvec.resize(2); return myvec;");
  auto df3 = df2.Define("LeadingJets_pz", "auto myvec(Jet.momentum.z); myvec.resize(2); return myvec;");
  auto df4 = df3.Define("LeadingJets_energy", "auto myvec(Jet.energy); myvec.resize(2); return myvec;");
  auto df5 = df4.Define("LeadingJets_P4", "auto P4 = Construct<ROOT::Math::PxPyPzEVector>(LeadingJets_px, LeadingJets_py, LeadingJets_pz, LeadingJets_energy); return P4;");
  auto df6 = df5.Define("Dijet_P4", "auto Q4 = LeadingJets_P4[0] + LeadingJets_P4[1]; return Q4;");
  auto df7 = df6.Define("Dijet_pt", "auto pt = Dijet_P4.Pt(); return pt;");
  auto df8 = df7.Define("Dijet_eta", "auto eta = Dijet_P4.Eta(); return eta;");
  auto rdf = df8.Define("Dijet_mass", "auto mass = Dijet_P4.M(); return mass;");
  



  auto h1 = rdf.Histo1D({"pT", "Dijet", 100u, 0., 100.}, "Dijet_pt");
  auto h2 = rdf.Histo1D({"M", "Dijet", 100u, 50., 250.}, "Dijet_mass");
  auto h3 = rdf.Histo1D({"Eta", "Dijet", 20u, -10., 10.}, "Dijet_eta");

  h1->GetXaxis()->SetTitle("Dijet transverse momentum");
  h2->GetXaxis()->SetTitle("Dijet invariant mass");
  h3->GetXaxis()->SetTitle("Dijet eta");

  auto c1 = new TCanvas("c1");
  h1->DrawCopy();

  //rdf.Display({"Dijet_P4", "Dijet_energy"})->Print();
  auto c2 = new TCanvas("c2");
  h2->DrawCopy();

  auto c3 = new TCanvas("c3");
  h3->DrawCopy();

}
