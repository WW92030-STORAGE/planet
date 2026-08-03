# PLANET

Lightweight planetary orbit simulator. Supports both static and gravity-influenced point masses.

# definitions

- `Planet` - a mass that is under the influence of external forces, and exerts forces on other masses.
- `Star` - a mass that only exerts forces on other masses but does not get influenced by external forces.

# information

Three integrators are provided: Forward Euler, Slightly Better Forward Euler, and Leapfrog. By default we use Leapfrog.
For a planet on a circular orbit, Forward Euler leaks energy, SBFE gains energy, and Leapfrog has varying but otherwise consistent energy.
