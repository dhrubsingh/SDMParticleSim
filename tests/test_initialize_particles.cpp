#include <cassert>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime> 
#include <cstdlib> 
#include <chrono>
#include <thread>
#include "tests.h"


std::vector<Particle> initializeParticles(const std::string& filename, int width, int height, float radius);

void test_initialize_particles(int n, const std::string& filename, int width, int height, float max_velocity, float radius) {
    
    std::vector<Particle> particles = initializeParticles(filename, width, height, radius);
    assert(particles.size() == n);
}

int main() {
    try {
        test_initialize_particles(675, "particles_675", 10, 10, 10, 2);
    } catch (const std::exception& e) {
        std::cerr << "Standard Error: " << e.what() << '\n';
    } catch (...) {
    std::cerr << "An unknown error occurred.\n";
}

    return 0;
}
