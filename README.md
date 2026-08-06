# PLANET

Lightweight planetary orbit simulator. Supports both static and gravity-influenced point masses.

# definitions

- `Planet` - a mass that is under the influence of external forces, and exerts forces on other masses.
- `Star` - a mass that only exerts forces on other masses but does not get influenced by external forces.

By these definitions a binary star system will consist of two `Planet` objects.

# information

Three integrators are provided: Forward Euler, Slightly Better Forward Euler, and Leapfrog. By default we use Leapfrog (since it is the most stable for this kind of simulation under sufficiently small timesteps). Be wary of energy instability due to precision errors.  

# miscellaneous

- Source for the Moore preset data - https://astronomy.stackexchange.com/questions/50297/initial-state-for-a-3-body-problem-to-create-figure-8-restricted-to-2d
- If you want real world data/estimates - https://ssd.jpl.nasa.gov/horizons/app.html#/