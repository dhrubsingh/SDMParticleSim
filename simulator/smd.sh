#!/usr/bin/env bash

#SBATCH --job-name=smd
#SBATCH --output=smd_%j.out
#SBATCH --error=smd_%j.err
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=9

# Load the necessary modules
spack load openmpi@4.1.6
module load openmpi

# Particle file names
particle_files=("../particle_files/particles_675" "../particle_files/particles_1250" "../particle_files/particles_2500" "../particle_files/particles_5000" "../particle_files/particles_10000" "../particle_files/particles_25000" "../particle_files/particles_50000" "../particle_files/particles_100000")

# Particle sizes
particle_sizes=(675 1250 2500 5000 10000 25000 50000 100000)

# Executable name
executable="particle_simulation"

# Compile the code
mpic++ -fopenmp -std=c++11 -o $executable smd.cpp

# Run simulations for different particle files and sizes
for i in "${!particle_files[@]}"; do
    particle_file="${particle_files[$i]}"
    particle_size="${particle_sizes[$i]}"
    
    echo "Running simulation for ${particle_file} (${particle_size} particles)"
    srun ./$executable $particle_size $particle_file
done
