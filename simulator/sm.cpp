#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>

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

// Read in particle positions and velocities from input file into particles vector
std::vector<Particle> initializeParticles(const std::string& filename, int width, int height, float radius) {
    std::vector<Particle> particles;
    std::ifstream file(filename);

    if (file.is_open()) {
        float x, y, vx, vy;
        while (file >> x >> y >> vx >> vy) {
            Particle particle = {x, y, vx, vy};
            if (isPlaceable(particles, x, y, radius * 2)) {
                particles.push_back(particle);
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    return particles;
}

// Function to update particle positions and velocities based on Lennard-Jones Potential
float updateParticles(std::vector<Particle>& particles, int width, int height, float epsilon, float sigma, float dt) {
    // Initialize acceleration vectors
    std::vector<float> ax(particles.size(), 0), ay(particles.size(), 0);
    std::vector<float> newAx(particles.size(), 0), newAy(particles.size(), 0);
    
    // Initialize potential and kinetic energy for given iteration
    float totalPotentialEnergy = 0.0f;
    float totalKineticEnergy = 0.0f;

    // Calculate initial accelerations based on Lennard-Jones formula 
    // Parallelized using dynamic scheduling
    #pragma omp parallel for schedule(dynamic) default(none) shared(particles, ax, ay, epsilon, sigma) reduction(+:totalPotentialEnergy)
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Particle& p1 = particles[i];
            Particle& p2 = particles[j];

            float dx = p2.x - p1.x;
            float dy = p2.y - p1.y;
            float rSquared = dx * dx + dy * dy;

            // Go to next iteration if particles are too far from each other or unrealistically close
            if (rSquared == 0 || sqrt(rSquared) > 2.5 * sigma) continue;

            float r2_inv = 1.0 / rSquared;
            float r6_inv = r2_inv * r2_inv * r2_inv;
            float forceScalar = 24 * epsilon * r6_inv * (2 * r6_inv - 1) * r2_inv;

            // Update accelerations for both particles involved
            ax[i] += dx * forceScalar;
            ay[i] += dy * forceScalar;
            ax[j] -= dx * forceScalar;
            ay[j] -= dy * forceScalar;

            // Compute potential energy of given iteration
            totalPotentialEnergy += 4 * epsilon * (r6_inv * r6_inv / r2_inv - r6_inv);
        }
    }

    // Update positions using dynamic scheduling
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].x += (particles[i].vx * dt) + (0.5f * ax[i] * dt * dt);
        particles[i].y += (particles[i].vy * dt) + (0.5f * ay[i] * dt * dt);
    }

    // Take care of collisions with boundary, essentially bouncing off the boundary
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < particles.size(); ++i) {
        if (particles[i].x < 0 || particles[i].x > width) {
            particles[i].vx *= -1;
            particles[i].x = std::max(0.0f, std::min(particles[i].x, static_cast<float>(width)));
        }
        if (particles[i].y < 0 || particles[i].y > height) {
            particles[i].vy *= -1;
            particles[i].y = std::max(0.0f, std::min(particles[i].y, static_cast<float>(height)));
        }
    }

    // Compute accelerations based on new positions using dynamic scheduling
    #pragma omp parallel for schedule(dynamic) default(none) shared(particles, newAx, newAy, epsilon, sigma)
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Particle& p1 = particles[i];
            Particle& p2 = particles[j];

            float dx = p2.x - p1.x;
            float dy = p2.y - p1.y;
            float rSquared = dx * dx + dy * dy;
            if (rSquared == 0 || sqrt(rSquared) > 2.5 * sigma) continue;

            float r2_inv = 1.0 / rSquared;
            float r6_inv = r2_inv * r2_inv * r2_inv;
            float forceScalar = 24 * epsilon * r6_inv * (2 * r6_inv - 1) * r2_inv;

            newAx[i] += dx * forceScalar;
            newAy[i] += dy * forceScalar;
            newAx[j] -= dx * forceScalar;
            newAy[j] -= dy * forceScalar;
        }
    }

    // Update velocities based on updated accelerations, parallelized with dynamic scheduling
    #pragma omp parallel for schedule(dynamic) reduction(+:totalKineticEnergy)
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].vx += 0.5f * (ax[i] + newAx[i]) * dt;
        particles[i].vy += 0.5f * (ay[i] + newAy[i]) * dt;

        // Compute kinetic energy of given iteration
        totalKineticEnergy += 0.5f * (particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy);
    }

    // Return total energy of system
    return totalKineticEnergy + totalPotentialEnergy;
}

// Function to ensure that particles do not overlap after their positions and velocities have been updated
void resolveDirectCollisions(std::vector<Particle>& particles, float radius, float dt) {
    for (size_t i = 0; i < particles.size(); i++) {
        for (size_t j = i + 1; j < particles.size(); j++) {
            Particle& p1 = particles[i];
            Particle& p2 = particles[j];

            float dx = p2.x - p1.x;
            float dy = p2.y - p1.y;
            float distanceSquared = dx * dx + dy * dy;
            float distance = sqrt(distanceSquared);
            float minDistance = 2 * radius;

            // Ensure that particles do not overlap
            if (distance < minDistance && distance > 0) {
                // Normalize vector between particles
                float nx = dx / distance;
                float ny = dy / distance;

                // Difference in particle velocities
                float vx = p1.vx - p2.vx;
                float vy = p1.vy - p2.vy;

                // Update particle positions if there is overlap
                float overlap = (minDistance - distance) / 2;
                p1.x -= overlap * nx;
                p1.y -= overlap * ny;
                p2.x += overlap * nx;
                p2.y += overlap * ny;

                // Particles should bounce off each other and move in opposite directions
                if (vx * nx + vy * ny < 0) {
                    float dot = (vx * nx + vy * ny) / (distance * distance);
                    p1.vx -= 2 * dot * dx;
                    p1.vy -= 2 * dot * dy;
                    p2.vx += 2 * dot * dx;
                    p2.vy += 2 * dot * dy;
                }
            }
        }
    }
}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <particle_size> <particle_file>" << std::endl;
        return 1;
    }

    // Read in number of particles and particle file
    int particleSize = std::atoi(argv[1]);
    std::string particleFilePath = argv[2];

    // Set initial parameters for Lennard-Jones simulation
    const int width = particleSize, height = particleSize;
    const float epsilon = 3.0f;
    const float sigma = 2.0f;
    const float particleRadius = 5.0f;

    // Initialize particles vector based on file input
    const std::string particleFile = particleFilePath;
    std::vector<Particle> particles = initializeParticles(particleFile, width, height, particleRadius);

    // Initialize variables to keep track of energy change over different iterations
    std::vector<float> energyHistory;
    const size_t energyHistorySize = 100;
    const float equilibriumThreshold = 0.01f;

    auto startTime = std::chrono::steady_clock::now();
    bool equilibriumReached = false;
    std::ofstream energyFile("total_energies.txt");

    // Update particle positions and velocities for 100 iterations
    for (int step = 0; step < 100; ++step) {
        float deltaTime = 0.016f;

        float totalEnergy = updateParticles(particles, width, height, epsilon, sigma, deltaTime);
        resolveDirectCollisions(particles, particleRadius, deltaTime);

        energyHistory.push_back(totalEnergy);
        if (energyHistory.size() > energyHistorySize) {
            energyHistory.erase(energyHistory.begin());
            float energyChange = std::abs(energyHistory.front() - energyHistory.back());

            if (energyChange < equilibriumThreshold) {
                equilibriumReached = true;
                std::cout << "Equilibrium reached after " << step << " steps." << std::endl;
            }
        }
        energyFile << step << ", " << totalEnergy << std::endl;
    }

    energyFile.close();
    auto endTime = std::chrono::steady_clock::now();
    std::chrono::duration<double> runtime = endTime - startTime;
    std::cout << "Simulation completed in " << runtime.count() << " seconds." << std::endl;

    return 0;
}
