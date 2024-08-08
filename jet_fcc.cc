void jet_fcc()
{
  //TFile *f = TFile::Open("events_029930132.root");
  //TTree *tree = f->Get<TTree>("events/Jet");
  using namespace ROOT;
  ROOT::RDataFrame df("events", "events_029930132.root");
  //ROOT::RDataFrame df(*tree);

  auto rdf = df.Define("Jet_pt", "sqrt((Jet.momentum.x)*(Jet.momentum.x) + (Jet.momentum.y)*(Jet.momentum.y))");

  auto h1 = rdf.Histo1D("Jet.energy");
  auto h2 = rdf.Histo1D("Jet.mass");
  auto h3 = rdf.Histo1D("Jet_pt");

  h1->GetXaxis()->SetTitle("Jet energy");
  h2->GetXaxis()->SetTitle("Jet mass");
  h3->GetXaxis()->SetTitle("Jet transverse momentum");

  auto c1 = new TCanvas("c1");
  h1->DrawCopy();

  auto c2 = new TCanvas("c2");
  h2->DrawCopy();

  auto c3 = new TCanvas("c3");
  h3->DrawCopy();

}
