## Simulating Thermal Equilibrium

<a name="readme-top"></a>


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#folder-structure">Folder Structure</a></li>
      </ul>
    </li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

Understanding molecular dynamics, and more broadly particle interaction, is important in a variety of fields including drug discovery, material science, biological mechanisms, and more. This project simulates and visualizes the dynamics of particles in a two-dimensional space using the Lennard-Jones potential to model inter-particle forces. The simulation calculates forces between particles based on the Lennard-Jones potential, which comprises attractive and repulsive components. Particle collisions are detected and resolved to ensure realistic motion, and the system's kinetic energy is tracked over time to monitor equilibrium. Importantly, this research is aimed at understanding the suitability of this problem scope to high performance computing (HPC). We evaluate 3 implementations, namely, a baseline sequential implementation, shared memory implementation, and a combined shared and distributed memory implementation, for various particle problem sizes. Overall, we find a performance speedup of 2.1 for our shared memory model and 6.2 for our combined shared and distributed memory model relative to the baseline, showing the benefit of applying HPC to this setting, and supporting the ability for future research at larger problem scales to be conducted.



### Folder Structure

* `glfw-3.4`
  * GUI Interface Library
* `milestone`
  * Presentations/Final Report
* `particle_files`
  * Generated particles data (Position and Velocity)
* `simulator`
  * Baseline, OMP, MPI + OMP Simulation
* `tests`
  * Unit Tests

