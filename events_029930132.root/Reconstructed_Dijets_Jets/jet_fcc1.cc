void jet_fcc1()
{ 
  using namespace ROOT;

  ROOT::RDataFrame df("events", "events_029930132.root");


  //Defining columns of root dataframe df
  auto df0 = df.Filter("Jet.momentum.x.size() > 1");
  auto df1 = df0.Define("Jetpt", "sqrt(Jet.momentum.x*Jet.momentum.x + Jet.momentum.y*Jet.momentum.y)");
  auto df2 = df1.Define("SortedJets_px", "auto sortIndices = ROOT::VecOps::Argsort(Jetpt, [](double x, double y) {return x > y;}); auto values = ROOT::VecOps::Take(Jet.momentum.x, sortIndices); return values;");
  auto df3 = df2.Define("SortedJets_py", "auto sortIndices = ROOT::VecOps::Argsort(Jetpt, [](double x, double y) {return x > y;}); auto values = ROOT::VecOps::Take(Jet.momentum.y, sortIndices); return values;");
  auto df4 = df3.Define("SortedJets_pz", "auto sortIndices = ROOT::VecOps::Argsort(Jetpt, [](double x, double y) {return x > y;}); auto values = ROOT::VecOps::Take(Jet.momentum.z, sortIndices); return values;");
  auto df5 = df4.Define("SortedJets_energy", "auto sortIndices = ROOT::VecOps::Argsort(Jetpt, [](double x, double y) {return x > y;}); auto values = ROOT::VecOps::Take(Jet.energy, sortIndices); return values;");
  auto df6 = df5.Define("Jets_P4", "auto P4 = Construct<ROOT::Math::PxPyPzEVector>(SortedJets_px, SortedJets_py, SortedJets_pz, SortedJets_energy); return P4;");
  auto df7 = df6.Define("Dijet_P4", "auto Q4 = Jets_P4[0] + Jets_P4[1]; return Q4;");
  auto df8 = df7.Define("Dijet_pt", "auto pt = Dijet_P4.Pt(); return pt;");
  auto df9 = df8.Define("Dijet_eta", "auto eta = Dijet_P4.Eta(); return eta;");
  auto rdf = df9.Define("Dijet_mass", "auto mass = Dijet_P4.M(); return mass;");
  
  auto h1 = rdf.Histo1D({"pT", "Dijet", 20u, 0., 20.}, "Dijet_pt");
  auto h2 = rdf.Histo1D({"M", "Dijet", 250u, 0., 250.}, "Dijet_mass");
  auto h3 = rdf.Histo1D({"Eta", "Dijet", 20u, -10., 10.}, "Dijet_eta");

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
