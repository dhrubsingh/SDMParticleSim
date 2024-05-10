#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
using namespace std;

struct Particle {
    float x, y;       // Position
    float vx, vy;     // Velocity
};

// Check if particles are a certain distance away from each other
bool isPlaceable(const std::vector<Particle>& particles, float x, float y, float min_distance) {
    for (const auto& particle : particles) {
        float dx = particle.x - x;
        float dy = particle.y - y;
        if (dx * dx + dy * dy < min_distance * min_distance) {
            return false;
        }
    }
    return true;
}

// Function to generate particles for a given 2D space
std::vector<Particle> initializeParticles(int n, int width, int height, float radius) {
    // Initialization for writing particles to file
    ofstream ParticleFile("particles_500");
    std::vector<Particle> particles;
    particles.reserve(n);
    float centerX = width / 2.0f;
    float centerY = height / 2.0f;
    
    // Random seed initialization
    srand(static_cast<unsigned int>(time(nullptr)));

    // Generate n particles
    for (int i = 0; i < n; ++i) {
        bool placed = false;
        int attempts = 0;
        while (!placed && attempts < 10000) {
            // Generate random position for given particle
            float x = (static_cast<float>(rand()) / RAND_MAX) * (width - 2 * radius) + radius;
            float y = (static_cast<float>(rand()) / RAND_MAX) * (height - 2 * radius) + radius;

            // Check if generated position is sufficient distance from current particles
            if (isPlaceable(particles, x, y, radius * 2)) {
                // Calculate distance from center of 2D space to set initial velocities for particles
                float dx = x - centerX;
                float dy = y - centerY;
                float distance = sqrt(dx * dx + dy * dy);
                float distanceFactor = distance / (sqrt(centerX * centerX + centerY * centerY));
                float baseSpeed = 50.0f * distanceFactor;

                // Generate velocity for particle and create particle
                Particle particle = {
                    x, y,
                    (static_cast<float>(rand()) / RAND_MAX - 0.5f) * 2.0f * baseSpeed,
                    (static_cast<float>(rand()) / RAND_MAX - 0.5f) * 2.0f * baseSpeed,
                };

                // Write the particle to the file
                ParticleFile << particle.x << ' ' << particle.y << ' ' << particle.vx << ' ' << particle.vy << '\n';

                // Add particle to vector
                particles.push_back(particle);
                placed = true;
            } else {
                ++attempts;
            }
        }

        // Break if unable to place particle
        if (!placed) {
            std::cerr << "Failed to place particle: " << i << std::endl;
            break;
        }
    }

    // Close file with particle positions and velocities
    ParticleFile.close();

    return particles;
}

int main() {
    srand(static_cast<unsigned int>(time(nullptr)));

    // Set initial parameters for generating particles
    const int width = 500, height = 500;
    const float particleRadius = 5.0f;
    const int particleCount = 500;

    // Initialize particles
    std::vector<Particle> particles = initializeParticles(particleCount, width, height, particleRadius);
}
