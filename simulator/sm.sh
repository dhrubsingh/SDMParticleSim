#!/usr/bin/env bash
#SBATCH --job-name=sm
#SBATCH --output=sm_%j.out
#SBATCH --error=sm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36

# Particle file names
particle_files=("../particle_files/particles_675" "../particle_files/particles_1250" "../particle_files/particles_2500" "../particle_files/particles_5000" "../particle_files/particles_10000" "../particle_files/particles_25000" "../particle_files/particles_50000" "../particle_files/particles_100000")

# Particle sizes
particle_sizes=(675 1250 2500 5000 10000 25000 50000 100000)

# Executable name
executable="particle_simulation"

# Compile the code
mpic++ -fopenmp -std=c++11 -o $executable sm.cpp

# Run simulations for different particle files and sizes
for i in "${!particle_files[@]}"; do
    particle_file="${particle_files[$i]}"
    particle_size="${particle_sizes[$i]}"

    echo "Running simulation for ${particle_file} (${particle_size} particles)"
    srun -n 1 -c 18 ./$executable $particle_size $particle_file
done
