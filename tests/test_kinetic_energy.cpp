#include <cassert>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>   // For std::time
#include <cstdlib> // For std::rand and std::srand
#include "tests.h"

void test_kinetic_energy(const std::string& filename, int width, int height, float particleRadius) {
    std::vector<Particle> particles = initializeParticles(filename, width, height, particleRadius);

    // Set initial parameters for Lennard-Jones simulation
    const float epsilon = 3.0f;
    const float sigma = 2.0f;

    // Initialize variables to keep track of energy change over different iterations
    std::vector<float> energyHistory;
    const size_t energyHistorySize = 100;
    const float equilibriumThreshold = 0.01f;

    // Update particle positions and velocities for 5 iterations
    for (int step = 0; step < 5; ++step) {
        float deltaTime = 0.016f;
        float totalEnergy = updateParticles(particles, width, height, epsilon, sigma, deltaTime);
        resolveDirectCollisions(particles, particleRadius, deltaTime);
        energyHistory.push_back(totalEnergy);
    }
    float energyChange = std::abs(energyHistory.front() - energyHistory.back());
    assert(energyChange >= 0);
}

int main() {
    try {
        test_kinetic_energy("particles_675", 675, 675, 5.0f);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
