#include <iostream>
#include "species.h"
#include "pitzer.h"
#include "solver.h"
#include "speciation.h"

int main () {
  std::vector<AqueousSpeies> species = loadspecies("database.txt");

for (auto& sp : species) {
  if (sp.name == "UO2++" || sp.name == "F-") {
  sp.isPrimary = true;
  }
}

Solution solution;
TotalInput input;
input.totals["U"] = 1e-3;
input.totals["F"] = 2e-3;

solution.concentrations["UO2++"] = 5e-4;
solutions.concentrations"F-" = 1e-3;

calculateActivityCoefficients(species, solution);
solveSpeciation(species, solution, input);

return 0;
