void file()
{
  using namespace ROOT;
  ROOT::RDataFrame df("events", "events_029930132.root");

  auto df0 = df.Define("Jet_pt", "sqrt((Jet.momentum.x)*(Jet.momentum.x) + (Jet.momentum.y)*(Jet.momentum.y))");
  auto rdf = df0.Define("Jet_n", "Jet.momentum.x.size()");

  auto h1 = rdf.Histo1D("Jet_pt");
  auto h2 = rdf.Histo1D({"N", "Jet", 7u, 0, 7},"Jet_n");

  h1->GetXaxis()->SetTitle("Jet pT");
  h2->GetXaxis()->SetTitle("Number of jets");
 

  auto c1 = new TCanvas("c1");
  h1->DrawCopy();  

  auto c2 = new TCanvas("c2");  
  h2->DrawCopy();

}
