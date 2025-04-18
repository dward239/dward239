#include "speciation.h"
#include <iostream>
#include <cmath>
#include <limits>
#include <set>
#include <map>
// Speciation solver stub: will implement equilibrium & mass balance later
void solveSpeciation(std::vector<AqueousSpecies>& species,
                     Solution& solution,
                     const TotalInput& input) {
    std::cout << "\n=== Speciation Solver ===\n";
    // STEP 1: Placeholder printout of total inputs
    std::cout << "Total element inputs (mol/L):\n";
    for (const auto& pair : input.totals) {
        std::cout << "  " << pair.first << ": " << pair.second << "\n";
    }
    // STEP 2: Detect elements and build stoichiometric matrix
    auto elements = detectUniqueElements(species);
    auto stoichMatrix = buildStoichiometricMatrix(species);
    // STEP 3: Compute residuals from current concentrations
    auto residuals = computeMassBalanceResiduals(species, solution, input, stoichMatrix);
    std::cout << "\nMass balance residuals:\n";
    for (const auto& r : residuals) {
        std::cout << "  " << r.first << ": " << r.second << "\n";
    }
    // STEP 4: TODO - Solve for equilibrium species concentrations (Newton-Raphson or iterative update)
}
// Extract elemental components from a species name or stoichiometry key
// For now, assume top-level species names or stoichiometry keys are elements or known labels
std::set<std::string> detectUniqueElements(const std::vector<AqueousSpecies>& species) {
    std::set<std::string> elements;
    for (const auto& sp : species) {
        for (const auto& comp : sp.stoichiometry) {
            const std::string& label = comp.first;
            // Crude heuristic: if the label contains no digits or symbols, assume it's an element
            bool isElement = true;
            for (char c : label) {
                if (isdigit(c) || c == '+' || c == '-' || c == '(' || c == ')') {
                    isElement = false;
                    break;
                }
            }
            if (isElement) {
                elements.insert(label);
            }
        }
    }
    return elements;
}
// Builds a stoichiometric coefficient matrix: element -> species -> coefficient
std::map<std::string, std::map<std::string, int>>
buildStoichiometricMatrix(const std::vector<AqueousSpecies>& species) {
    std::map<std::string, std::map<std::string, int>> matrix;
    for (const auto& sp : species) {
        for (const auto& comp : sp.stoichiometry) {
            const std::string& element = comp.first;
            int coeff = comp.second;
            matrix[element][sp.name] += coeff;
        }
    }
    return matrix;
}
// Compute mass balance residuals: element-wise difference between total input and sum of species contributions
std::map<std::string, double>
computeMassBalanceResiduals(const std::vector<AqueousSpecies>& species,
                            const Solution& solution,
                            const TotalInput& input,
                            const std::map<std::string, std::map<std::string, int>>& stoichMatrix) {
    std::map<std::string, double> residuals;
    for (const auto& elementPair : stoichMatrix) {
        const std::string& element = elementPair.first;
        double total = input.totals.count(element) ? input.totals.at(element) : 0.0;
        double sum = 0.0;
        for (const auto& speciesPair : elementPair.second) {
            const std::string& speciesName = speciesPair.first;
            int coeff = speciesPair.second;
            auto it = solution.concentrations.find(speciesName);
            if (it != solution.concentrations.end()) {
                sum += coeff * it->second;
            }
        }
        residuals[element] = sum - total;
    }
    return residuals;
}
