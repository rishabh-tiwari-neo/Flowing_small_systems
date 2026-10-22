#include "TFile.h"
#include "TParameter.h"
void add_g(const char* fname, double g) {
    TFile f(fname, "UPDATE");
    TParameter<double> p("g", g);
    p.Write();
    f.Close();
}
