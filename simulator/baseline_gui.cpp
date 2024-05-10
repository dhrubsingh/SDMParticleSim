#define GL_SILENCE_DEPRECATION
#include "../glfw-3.4/include/GLFW/glfw3.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <chrono>

struct Particle {
    float x, y;       // Position
    float vx, vy;     // Velocity
    float r, g, b;    // Color
};

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

std::vector<Particle> initializeParticles(int n, int width, int height, float radius) {
    std::vector<Particle> particles;
    particles.reserve(n);
    float centerX = width / 2.0f;
    float centerY = height / 2.0f;

    for (int i = 0; i < n; ++i) {
        int attempts = 0;
        bool placed = false;
        while (!placed && attempts < 10000) {
            float x = (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) * (width - 2 * radius) + radius;
            float y = (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) * (height - 2 * radius) + radius;

            if (isPlaceable(particles, x, y, radius * 3)) {
                float dx = x - centerX;
                float dy = y - centerY;
                float distance = sqrt(dx * dx + dy * dy);
                float distanceFactor = distance / (sqrt(centerX * centerX + centerY * centerY));
                float baseSpeed = 50.0f;  // Adjust this value as needed

                Particle p = {
                x, y,  // Position
                (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) - 0.5f) * 2.0f * baseSpeed,  // x velocity
                (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) - 0.5f) * 2.0f * baseSpeed,  // y velocity
                static_cast<float>(rand()) / static_cast<float>(RAND_MAX),  // Red color component
                static_cast<float>(rand()) / static_cast<float>(RAND_MAX),  // Green color component
                static_cast<float>(rand()) / static_cast<float>(RAND_MAX)   // Blue color component
            };

                particles.push_back(p);
                placed = true;
            } else {
                ++attempts;
            }
        }

        if (!placed) {
            std::cerr << "Unable to place all particles without overlap." << std::endl;
            break;
        }
    }

    return particles;
}


void drawCircle(float cx, float cy, float r, int num_segments) {
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(cx, cy); // center of circle
    for (int ii = 0; ii <= num_segments; ++ii) {
        float theta = 2.0f * 3.1415926f * static_cast<float>(ii) / static_cast<float>(num_segments); // Current angle
        float x = r * cosf(theta); // Calculate the x component
        float y = r * sinf(theta); // Calculate the y component
        glVertex2f(x + cx, y + cy); // Output vertex
    }
    glEnd();
}

void handleCollision(Particle& p1, Particle& p2, float radius) {
    // Calculate displacement from p1 to p2
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;

    // Calculate distance between particles
    float distance = sqrt(dx * dx + dy * dy);

    // Normalizing vector components (dx, dy) to get direction of impact
    float nx = dx / distance;
    float ny = dy / distance;

    // Calculate velocity component along the normal direction (dot product)
    float p1n = nx * p1.vx + ny * p1.vy;
    float p2n = nx * p2.vx + ny * p2.vy;

    // Conservation of momentum in 1D for each component
    float p1nFinal = p2n;
    float p2nFinal = p1n;

    // Update velocities
    p1.vx += (p1nFinal - p1n) * nx;
    p1.vy += (p1nFinal - p1n) * ny;
    p2.vx += (p2nFinal - p2n) * nx;
    p2.vy += (p2nFinal - p2n) * ny;

    float overlap = 2 * radius - distance; // distance is the distance between particles
    p1.x -= overlap / 2 * nx;
    p1.y -= overlap / 2 * ny;
    p2.x += overlap / 2 * nx;
    p2.y += overlap / 2 * ny;
}


float updateParticles(std::vector<Particle>& particles, int width, int height, float epsilon, float sigma, float dt) {
    std::vector<float> ax(particles.size(), 0), ay(particles.size(), 0);
    float totalPotentialEnergy = 0.0f; // Initialize total potential energy for this step

    // Step 1: Calculate initial accelerations based on Lennard-Jones potential
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
            float r6 = rSquared * rSquared * rSquared;
            float r12 = r6 * r6;
            float sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
            float sigma12 = sigma6 * sigma6;

            ax[i] += dx * forceScalar;
            ay[i] += dy * forceScalar;
            ax[j] -= dx * forceScalar;
            ay[j] -= dy * forceScalar;

            // Calculate and accumulate potential energy

            
             
            totalPotentialEnergy += 4 * epsilon * (sigma12 / r12 - sigma6 / r6);
        }
    }

    // Step 2: Update positions based on initial velocities and accelerations
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].x += particles[i].vx * dt + 0.5f * ax[i] * dt * dt;
        particles[i].y += particles[i].vy * dt + 0.5f * ay[i] * dt * dt;
    }

    // Boundary conditions (Reflection)
    for (auto& p : particles) {
        if (p.x < 0 || p.x > width) {
            p.vx *= -1;
            p.x = std::max(0.0f, std::min(p.x, static_cast<float>(width)));
        }
        if (p.y < 0 || p.y > height) {
            p.vy *= -1;
            p.y = std::max(0.0f, std::min(p.y, static_cast<float>(height)));
        }
    }

    // Recalculate accelerations for new positions (Velocity Verlet Step 3)
    std::vector<float> newAx(particles.size(), 0), newAy(particles.size(), 0);
    // Similar loop to calculate newAx, newAy based on updated positions

    // Step 4: Update velocities with new accelerations
    float totalKineticEnergy = 0.0f;
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].vx += 0.5f * (ax[i] + newAx[i]) * dt;
        particles[i].vy += 0.5f * (ay[i] + newAy[i]) * dt;

        // Calculate kinetic energy for each particle
        totalKineticEnergy += 0.5 * (particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy);
    }

    std::cout << "Kinetic: " << totalKineticEnergy << std::endl;
    std::cout << "Potential: " << totalPotentialEnergy << std::endl;
    // Return the total energy (kinetic + potential)
    return totalKineticEnergy + totalPotentialEnergy;
}



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

            // Check for overlap
            if (distance < minDistance && distance > 0) {
                // Normal vector between particles
                float nx = dx / distance;
                float ny = dy / distance;

                // Resolve penetration
                float overlap = (minDistance - distance) / 2;
                p1.x -= overlap * nx;
                p1.y -= overlap * ny;
                p2.x += overlap * nx;
                p2.y += overlap * ny;

                // Calculate relative velocity
                float vx = p1.vx - p2.vx;
                float vy = p1.vy - p2.vy;

                // Perform reflection if particles are moving towards each other
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



int main() {
    srand(static_cast<unsigned int>(time(nullptr)));

    if (!glfwInit()) {
        return -1;
    }

    GLFWwindow* window = glfwCreateWindow(640, 480, "Particle Visualization", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    const float epsilon = 3.0f;
    const float sigma = 2.0f;
    const float particleRadius = 10.0f;

    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    std::vector<Particle> particles = initializeParticles(5000, width, height, particleRadius);

    std::vector<float> energyHistory;
    const size_t energyHistorySize = 100; // Number of steps to consider for equilibrium
    const float equilibriumThreshold = 0.01f; // Energy change threshold
    bool equilibriumReached = false;

    double lastTime = glfwGetTime();

    using namespace std::chrono;
    auto t0 = steady_clock::now();

    while (!glfwWindowShouldClose(window)) {
        double currentTime = glfwGetTime();
        float deltaTime = static_cast<float>(currentTime - lastTime);
        lastTime = currentTime;

        float totalEnergy = updateParticles(particles, width, height, epsilon, sigma, deltaTime);
        resolveDirectCollisions(particles, particleRadius, deltaTime);

        // Debugging: Print the total energy of the system
        std::cout << "Total Energy: " << totalEnergy << std::endl;

        energyHistory.push_back(totalEnergy);
        if (energyHistory.size() > energyHistorySize) {
            energyHistory.erase(energyHistory.begin()); // Keep history size fixed

            // Check for equilibrium
            float energyChange = std::abs(energyHistory.front() - energyHistory.back());
            
            // Debugging: Print the change in energy
            std::cout << "Energy Change (over last " << energyHistorySize << " steps): " << energyChange << std::endl;

            if (energyChange < equilibriumThreshold) {
                std::cout << "Equilibrium reached." << std::endl;
                break;
            }
        }

        glClear(GL_COLOR_BUFFER_BIT);
        glViewport(0, 0, width, height);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0.0f, static_cast<float>(width), static_cast<float>(height), 0.0f, -1.0f, 1.0f);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        for (const auto& particle : particles) {
            glColor3f(particle.r, particle.g, particle.b);
            drawCircle(particle.x, particle.y, particleRadius, 20);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();

    auto t1 = steady_clock::now();
    const double t = duration_cast<microseconds>(t1 - t0).count() / 1000000.0;

    printf("\tRuntime: %e s\n", t);

    return 0;
}




//g++ -std=c++11 -I../glfw-3.4/include/GLFW/ baseline.cpp -o baseline -framework OpenGL -L../glfw-3.4/lib-arm64/ -lglfw3 -framework Cocoa -framework IOKit -framework CoreVideo

