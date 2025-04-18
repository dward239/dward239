#include "pitzer.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
// Load Pitzer binary interaction parameters from a file
std::vector<PitzerPair> loadPitzerParameters() {
   std::vector<PitzerPair> params;
   std::ifstream file("pitzer_params.txt");
   std::string line;
   if (!file.is_open()) {
       std::cerr << "Error: Could not open pitzer_params.txt\n";
       return params;
   }
   while (getline(file, line)) {
       if (line.empty() || line[0] == '#') continue;
       std::istringstream iss(line);
       PitzerPair p;
       iss >> p.ion1 >> p.ion2 >> p.beta0 >> p.beta1 >> p.Cphi;
       params.push_back(p);
   }
   return params;
}
// Calculate ionic strength I = 0.5 * sum(ci * zi^2)
double calculateIonicStrength(const Solution& solution, const std::map<std::string, double>& charges) {
   double I = 0.0;
   for (const auto& pair : solution.concentrations) {
       const std::string& ion = pair.first;
       double conc = pair.second;
       auto it = charges.find(ion);
       if (it != charges.end()) {
           double z = it->second;
           I += conc * z * z;
       }
   }
   return 0.5 * I;
}
// Apply a simplified version of Pitzer activity coefficient calculation (placeholder)
void calculateActivityCoefficients(std::vector<AqueousSpecies>& species, Solution& solution) {
   auto pitzerParams = loadPitzerParameters();
   // Build a quick charge map
   std::map<std::string, double> charges;
   for (const auto& sp : species) {
       charges[sp.name] = sp.charge;
   }
   // Calculate ionic strength first
   solution.ionicStrength = calculateIonicStrength(solution, charges);
   std::cout << "Ionic strength (Pitzer): " << solution.ionicStrength << "\n";
   // For now, use gamma = 1 as a placeholder
   for (auto& sp : species) {
       sp.activityCoefficient = 1.0;
   }
   // (Optional) Apply full Pitzer equation later here
}
