# PLANET

Lightweight planetary orbit simulator. Supports both static and gravity-influenced point masses.

# definitions

- `Planet` - a mass that is under the influence of external forces, and exerts forces on other masses.
- `Star` - a mass that only exerts forces on other masses but does not get influenced by external forces.

# information

Three integrators are provided: Forward Euler, Slightly Better Forward Euler, and Leapfrog. By default we use Leapfrog (since it is the most precise). Be wary of energy instability due to precision errors.

# miscellaneous

https://www.math.uni-bielefeld.de/~rehmann/ECM/cdrom/3ecm/pdfs/pant3/simo.pdf