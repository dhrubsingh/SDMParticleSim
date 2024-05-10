#include <vector>    // Correct header for std::vector
#include <fstream>   // For std::ifstream
#include <iostream>  // For std::cerr
#include <string>    // For std::string

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
