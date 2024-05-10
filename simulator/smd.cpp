#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <algorithm>

struct Particle {
    float x, y;       // Position
    float vx, vy;     // Velocity
    float r, g, b;    // Color
    float radius; // particle radius
};

float updateSectionParticles(std::vector<Particle>& particles, int sectionId, int width, int height, float epsilon, float sigma, float dt);
void exchangeParticleData(std::vector<Particle>& particles, int sectionId, int width, int height, float particleRadius);

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

// Function to initialize particles, reading the input file and filtering for the specific quadrant
std::vector<Particle> initializeParticles(int sectionId, int numSections, const std::string& filename, int width, int height, float radius) {
    std::vector<Particle> particles;

    // Initialize the boundaries of the quadrant
    int midWidth = width / 2;
    int midHeight = height / 2;
    int startX, endX, startY, endY;

    // Set horizontal boundaries for quadrant
    if (sectionId % 2 == 0) {
        startX = 0;
        endX = midWidth;
    } else {
        startX = midWidth;
        endX = width;
    }

    // Set vertical boundaries for quadrant
    if (sectionId > 1) {
        startY = 0;
        endY = midHeight;
    } else {
        startY = midHeight;
        endY = height;
    }

    // Read in particle positions and velocities from input file into particles vector
    std::ifstream file(filename);
    if (file.is_open()) {
        float x, y, vx, vy;
        while (file >> x >> y >> vx >> vy) {
            // Check if particle is within quadrant boundaries
            if (x >= startX && x < endX && y >= startY && y < endY) {
                Particle particle = {x, y, vx, vy};
                particles.push_back(particle);
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    return particles;
}

// Function to update particle positions and velocities
float updateParticles(std::vector<Particle>& particles, int sectionId, int width, int height, float epsilon, float sigma, float dt, float particleRadius) {
    float energy;

    // Update particles within quadrant
    energy = updateSectionParticles(particles, sectionId, width, height, epsilon, sigma, dt);

    // Exchange particle data with neighboring sections as necessary
    exchangeParticleData(particles, sectionId, width, height, particleRadius);

    // Return total energy of system
    return energy;
}

// Function to update particle positions and velocities based on Lennard-Jones Potential
float updateSectionParticles(std::vector<Particle>& particles, int sectionId, int width, int height, float epsilon, float sigma, float dt) {
    // Initialize acceleration vectors and potential energy
    std::vector<float> ax(particles.size(), 0), ay(particles.size(), 0);
    float totalPotentialEnergy = 0.0f;

    // Calculate initial accelerations based on Lennard-Jones formula, parallelized using dynamic scheduling
    #pragma omp parallel for schedule(dynamic) default(none) shared(particles, ax, ay, epsilon, sigma) reduction(+:totalPotentialEnergy)
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Particle& p1 = particles[i];
            Particle& p2 = particles[j];

            float dx = p2.x - p1.x;
            float dy = p2.y - p1.y;
            float rSquared = dx * dx + dy * dy;
            if (rSquared == 0 || sqrt(rSquared) > 2.5 * sigma) continue; // Avoid division by zero and distant particles

            float r2_inv = 1.0 / rSquared;
            float r6_inv = r2_inv * r2_inv * r2_inv;
            float forceScalar = 24 * epsilon * r6_inv * (2 * r6_inv - 1) * r2_inv;

            ax[i] += dx * forceScalar;
            ay[i] += dy * forceScalar;
            ax[j] -= dx * forceScalar;
            ay[j] -= dy * forceScalar;

            totalPotentialEnergy += 4 * epsilon * (r6_inv * r6_inv / r2_inv - r6_inv);
        }
    }

    // Update positions using dynamic scheduling
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].x += particles[i].vx * dt + 0.5f * ax[i] * dt * dt;
        particles[i].y += particles[i].vy * dt + 0.5f * ay[i] * dt * dt;
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
    std::vector<float> newAx(particles.size(), 0), newAy(particles.size(), 0);
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
    float totalKineticEnergy = 0.0f;
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

// Function to send particle data if particle moves to a neighboring quadrant
void exchangeParticleData(std::vector<Particle>& particles, int sectionId, int width, int height, float particleRadius) {
    // Initialize neighbor variables
    int neighborIds[2];
    int numNeighbors = 0;

    // Determine neighbor section ID based on quadrant-based layout
    if (sectionId == 0) {
        neighborIds[numNeighbors++] = 1; // Right neighbor
        neighborIds[numNeighbors++] = 2; // Bottom neighbor
    } else if (sectionId == 1) {
        neighborIds[numNeighbors++] = 0; // Left neighbor
        neighborIds[numNeighbors++] = 3; // Bottom-right neighbor
    } else if (sectionId == 2) {
        neighborIds[numNeighbors++] = 0; // Top neighbor
        neighborIds[numNeighbors++] = 3; // Right neighbor
    } else if (sectionId == 3) {
        neighborIds[numNeighbors++] = 1; // Top-left neighbor
        neighborIds[numNeighbors++] = 2; // Left neighbor
    }

    std::vector<MPI_Request> requests;

    // Initialize vectors for particles based on which neighbor it belongs in
    std::vector<Particle> sendParticles_neighbor1;
    std::vector<Particle> sendParticles_neighbor2;

    // Determine which particles belong to which neighboring quadrant
    for (const auto& particle : particles) {
        if (sectionId == 0) {
            if (particle.x >= width / 2) {
                sendParticles_neighbor1.push_back(particle); // Right neighbor
            }
            if (particle.y >= height / 2) {
                sendParticles_neighbor2.push_back(particle); // Bottom neighbor
            }
        } else if (sectionId == 1) {
            if (particle.x < width / 2) {
                sendParticles_neighbor1.push_back(particle); // Left neighbor
            }
            if (particle.y >= height / 2) {
                sendParticles_neighbor2.push_back(particle); // Bottom-right neighbor
            }
        } else if (sectionId == 2) {
            if (particle.y < height / 2) {
                sendParticles_neighbor1.push_back(particle); // Top neighbor
            }
            if (particle.x >= width / 2) {
                sendParticles_neighbor2.push_back(particle); // Right neighbor
            }
        } else if (sectionId == 3) {
            if (particle.y < height / 2) {
                sendParticles_neighbor1.push_back(particle); // Top-left neighbor
            }
            if (particle.x < width / 2) {
                sendParticles_neighbor2.push_back(particle); // Left neighbor
            }
        }
    }

    // Create MPI datatype for particles
    MPI_Datatype particleType;
    MPI_Type_contiguous(7, MPI_FLOAT, &particleType);
    MPI_Type_commit(&particleType);

    // Initialization to send number of particles to be sent to each neighboring quadrant
    int sendCounts[2] = {static_cast<int>(sendParticles_neighbor1.size()), static_cast<int>(sendParticles_neighbor2.size())};
    int recvCounts[2];
    
    // Communicating number of particles to neighboring quadrants
    for (int i = 0; i < numNeighbors; ++i) {
        MPI_Request request;
        MPI_Isend(&sendCounts[i], 1, MPI_INT, neighborIds[i], 0, MPI_COMM_WORLD, &request);
        requests.push_back(request);
        MPI_Irecv(&recvCounts[i], 1, MPI_INT, neighborIds[i], 0, MPI_COMM_WORLD, &request);
        requests.push_back(request);
    }

    // Wait for particle counts to be communicated before continuing
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    requests.clear();

    // Initialize receive buffers
    std::vector<Particle> recvBuffers[2];
    for (int i = 0; i < numNeighbors; ++i) {
        recvBuffers[i].resize(recvCounts[i]);
    }

    // Communicate particle positions and velocities to neighboring quadrants
    for (int i = 0; i < numNeighbors; ++i) {
        MPI_Request request;
        if (i == 0) {
            MPI_Isend(sendParticles_neighbor1.data(), sendCounts[i], particleType, neighborIds[i], 1, MPI_COMM_WORLD, &request);
        } else {
            MPI_Isend(sendParticles_neighbor2.data(), sendCounts[i], particleType, neighborIds[i], 1, MPI_COMM_WORLD, &request);
        }
        requests.push_back(request);
        MPI_Irecv(recvBuffers[i].data(), recvCounts[i], particleType, neighborIds[i], 1, MPI_COMM_WORLD, &request);
        requests.push_back(request);
    }

    // Wait for particle positions and velocities to be communicated before continuing
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    // Remove particles that are no longer in the quadrant
    particles.erase(std::remove_if(particles.begin(), particles.end(), [&](const Particle& particle) {
        for (const auto& sendParticle : sendParticles_neighbor1) {
            if (particle.x == sendParticle.x && particle.y == sendParticle.y) {
                return true;
            }
        }
        for (const auto& sendParticle : sendParticles_neighbor2) {
            if (particle.x == sendParticle.x && particle.y == sendParticle.y) {
                return true;
            }
        }
        return false;
    }), particles.end());

    // Add particles sent from neighboring quadrants
    for (const auto& recvBuffer : recvBuffers) {
        particles.insert(particles.end(), recvBuffer.begin(), recvBuffer.end());
    }

    MPI_Type_free(&particleType);
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

                // Update particle positions if there is overlap
                float overlap = (minDistance - distance) / 2;
                p1.x -= overlap * nx;
                p1.y -= overlap * ny;
                p2.x += overlap * nx;
                p2.y += overlap * ny;

                // Difference in particle velocities
                float vx = p1.vx - p2.vx;
                float vy = p1.vy - p2.vy;

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

// Function to take care of collisions that occur between two particles in neighboring quadrants
void resolveQuadrantBoundaryCollisions(std::vector<Particle>& particles, int sectionId, int width, int height, float dt) {
    // Initialize neighbor variables
    int neighborIds[2];
    int numNeighbors = 0;

    // Determine neighbor section ID based on quadrant-based layout
    if (sectionId == 0) {
        neighborIds[numNeighbors++] = 1; // Right neighbor
        neighborIds[numNeighbors++] = 2; // Bottom neighbor
    } else if (sectionId == 1) {
        neighborIds[numNeighbors++] = 0; // Left neighbor
        neighborIds[numNeighbors++] = 3; // Bottom-right neighbor
    } else if (sectionId == 2) {
        neighborIds[numNeighbors++] = 0; // Top neighbor
        neighborIds[numNeighbors++] = 3; // Right neighbor
    } else if (sectionId == 3) {
        neighborIds[numNeighbors++] = 1; // Top-left neighbor
        neighborIds[numNeighbors++] = 2; // Left neighbor
    }

    // Resolve collisions with particles in boundary regions of neighboring quadrants
    for (int i = 0; i < numNeighbors; ++i) {
        // Initialize vector for particles in boundary regions
        std::vector<Particle> boundaryParticles;

        // Determine particles in the boundary region with the current neighbor
        for (const auto& particle : particles) {
            if (sectionId == 0) {
                if ((i == 0 && particle.x >= width / 2 - particle.radius) ||
                    (i == 1 && particle.y >= height / 2 - particle.radius)) {
                    boundaryParticles.push_back(particle);
                }
            } else if (sectionId == 1) {
                if ((i == 0 && particle.x < particle.radius) ||
                    (i == 1 && particle.y >= height / 2 - particle.radius)) {
                    boundaryParticles.push_back(particle);
                }
            } else if (sectionId == 2) {
                if ((i == 0 && particle.y < particle.radius) ||
                    (i == 1 && particle.x >= width / 2 - particle.radius)) {
                    boundaryParticles.push_back(particle);
                }
            } else if (sectionId == 3) {
                if ((i == 0 && particle.y < particle.radius) ||
                    (i == 1 && particle.x < particle.radius)) {
                    boundaryParticles.push_back(particle);
                }
            }
        }

        // Initialization to send number of particles to be sent to each neighboring quadrant
        int neighborId = neighborIds[i];
        int sendCount = boundaryParticles.size();
        int recvCount;

        // Sending number of particles to neighboring quadrant
        MPI_Request requests[2];
        MPI_Isend(&sendCount, 1, MPI_INT, neighborId, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(&recvCount, 1, MPI_INT, neighborId, 0, MPI_COMM_WORLD, &requests[1]);
        MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

        // Initialize receive buffer
        std::vector<Particle> recvParticles(recvCount);

        // Create MPI datatype for particles
        MPI_Datatype particleType;
        MPI_Type_contiguous(7, MPI_FLOAT, &particleType);
        MPI_Type_commit(&particleType);

        // Communicate particle positions and velocities to neighboring quadrants
        MPI_Isend(boundaryParticles.data(), sendCount, particleType, neighborId, 1, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(recvParticles.data(), recvCount, particleType, neighborId, 1, MPI_COMM_WORLD, &requests[1]);
        
        // Wait for particle positions and velocities to be communicated before continuing
        MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

        MPI_Type_free(&particleType);

        // Ensure boundary particles and received particles don't overlap
        for (auto& p1 : boundaryParticles) {
            for (const auto& p2 : recvParticles) {
                float dx = p2.x - p1.x;
                float dy = p2.y - p1.y;
                float distanceSquared = dx * dx + dy * dy;
                float distance = std::sqrt(distanceSquared);
                float minDistance = p1.radius + p2.radius;

                // Checking for overlap between particles
                if (distance < minDistance && distance > 0) {
                    // Normalize vector between particles
                    float overlap = (minDistance - distance) / 2;
                    float nx = dx / distance;
                    float ny = dy / distance;

                    // Update particle positions if there is overlap
                    p1.x -= overlap * nx;
                    p1.y -= overlap * ny;

                    // Difference in particle velocities
                    float vx = p1.vx - p2.vx;
                    float vy = p1.vy - p2.vy;

                    // Particles should bounce off each other and move in opposite directions
                    if (vx * nx + vy * ny < 0) {
                        float dot = (vx * nx + vy * ny) / (distance * distance);
                        p1.vx -= 2 * dot * dx;
                        p1.vy -= 2 * dot * dy;
                    }
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {
    // MPI Initialization
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <particle_size> <particle_file>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    // Read in number of particles and particle file
    int particleSize = std::atoi(argv[1]);
    std::string particleFilePath = argv[2];

    // Determine section ID based on rank and number of nodes
    int sectionId = rank % 4;

    // Set initial parameters for Lennard-Jones simulation
    const float epsilon = 3.0f;
    const float sigma = 2.0f;
    const float particleRadius = 5.0f;

    // Initialize particles vector based on file input
    std::vector<Particle> particles = initializeParticles(sectionId, 4, particleFilePath, particleSize, particleSize, particleRadius);

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

        float totalEnergy = updateParticles(particles, sectionId, particleSize, particleSize, epsilon, sigma, deltaTime, particleRadius);
        resolveDirectCollisions(particles, particleRadius, deltaTime);
        resolveQuadrantBoundaryCollisions(particles, sectionId, particleSize, particleSize, deltaTime);

        std::vector<float> allEnergies(size);
        MPI_Allgather(&totalEnergy, 1, MPI_FLOAT, allEnergies.data(), 1, MPI_FLOAT, MPI_COMM_WORLD);

        if (rank == 0) {
            energyHistory.push_back(totalEnergy);
            if (energyHistory.size() > energyHistorySize) {
                energyHistory.erase(energyHistory.begin());
                float energyChange = std::abs(energyHistory.front() - energyHistory.back());

                if (energyChange < equilibriumThreshold) {
                    equilibriumReached = true;
                    std::cout << "Equilibrium reached after " << step << " steps." << std::endl;
                }
            }

            if (step % 1 == 0) {
                energyFile << step << ", " << totalEnergy << std::endl;
            }
        }
    }

    energyFile.close();
    auto endTime = std::chrono::steady_clock::now();
    std::chrono::duration<double> runtime = endTime - startTime;
   
    if (rank == 0) {
        std::cout << "Simulation completed in " << runtime.count() << " seconds." << std::endl;
    }

    MPI_Finalize();
    return 0;
}
