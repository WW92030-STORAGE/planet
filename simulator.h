#ifndef PLANETSIM_SIMULATOR_NORMEXT
#define PLANETSIM_SIMULATOR_NORMEXT

#include "planet.h"
#include "star.h"
#include "util.h"

#include <vector>
#include <string>

struct PlanetSimulator {
    std::vector<Planet> planets;
    std::vector<PSUtil::point> forces;
    std::vector<Star> stars;

    PSUtil::numeric BOUNDARY = 1024;

    void removePlanet(int index) {
        if (index >= planets.size() || index < 0) return;
        std::swap(planets[index], planets[planets.size() - 1]);
        planets.pop_back();

        // experimental: also pop the force
        if (index >= forces.size() || forces.size() != planets.size()) return;
        std::swap(forces[index], forces[forces.size() - 1]);
        forces.pop_back();
    }

    void removeStar(int index) {
        if (index >= stars.size() || index < 0) return;
        std::swap(stars[index], stars[stars.size() - 1]);
        stars.pop_back();
    }

    void computeForces() {
        int P = planets.size();
        int S = stars.size();

        // zero out
        while (forces.size() > P) forces.pop_back();
        for (int i = 0; i < forces.size(); i++) forces[i] = {0, 0};
        while (forces.size() < P) forces.push_back({0, 0});

        // compute. right now it's O(P * P + P * S)

        // compute planet-star forces

        for (int p = 0; p < P; p++) {
            for (int s = 0; s < S; s++) {
                Planet planet = planets[p];
                Star star = stars[s];

                PSUtil::point disp = PSUtil::sub(star.x, planet.x);
                PSUtil::numeric len = PSUtil::length(disp);
                if (PSUtil::zero(len)) continue;
                PSUtil::numeric force_magnitude = PSUtil::gravity(planet.mass, star.mass, len);

                // unit vector + add
                forces[p] = PSUtil::add(forces[p], PSUtil::scale(disp, force_magnitude / len));
            }
        }

        // compute planet-planet forces. Here, we use "star" to denote what is influencing the planet.

        for (int p = 0; p < P; p++) {
            for (int s = 0; s < P; s++) {
                if (s == p) continue;
                Planet planet = planets[p];
                Planet star = planets[s];

                PSUtil::point disp = PSUtil::sub(star.x, planet.x);
                PSUtil::numeric len = PSUtil::length(disp);
                if (PSUtil::zero(len)) continue;
                PSUtil::numeric force_magnitude = PSUtil::gravity(planet.mass, star.mass, len);

                // unit vector + add
                forces[p] = PSUtil::add(forces[p], PSUtil::scale(disp, force_magnitude / len));
            }
        }
    }

    // Store the correct accelerations. Assumes you've computed forces beforehand.
    void updateAccelerations() {
        int P = planets.size();
        for (int p = 0; p < P && p < forces.size(); p++) {
            if (!PSUtil::zero(planets[p].mass)) planets[p].a = PSUtil::scale(forces[p], 1.0 / planets[p].mass);
        }
    }

    // assumes you've computed forces beforehand. otherwise use the iterate method below. There is no guarantee that the forces array is preserved after this.
    void update(PSUtil::numeric delta) {
        updateAccelerations();

        updateLeapfrog(delta);
    }

    void updateEuler(PSUtil::numeric delta) {
        int P = planets.size();
        for (int p = 0; p < P; p++) {
            planets[p].v = PSUtil::add(planets[p].v, PSUtil::scale(planets[p].a, delta));
            planets[p].x = PSUtil::add(planets[p].x, PSUtil::scale(planets[p].v, delta));
        }
    }

    void updateBetter(PSUtil::numeric delta) {
        auto hd2 = delta * delta * 0.5;
        int P = planets.size();

        for (int p = 0; p < P; p++) {
            auto displacement = PSUtil::add(PSUtil::scale(planets[p].v, delta), PSUtil::scale(planets[p].a, hd2));
            planets[p].x = PSUtil::add(planets[p].x, displacement);
            planets[p].v = PSUtil::add(planets[p].v, PSUtil::scale(planets[p].a, delta));
        }
    }

    // WARNING: This integrator will compute forces again during the step.
    void updateLeapfrog(PSUtil::numeric delta) {
        auto hd = delta * 0.5;
        int P = planets.size();
        for (int p = 0; p < P; p++) {
            planets[p].v = PSUtil::add(planets[p].v, PSUtil::scale(planets[p].a, hd));
            planets[p].x = PSUtil::add(planets[p].x, PSUtil::scale(planets[p].v, delta));
        }
        computeForces();
        updateAccelerations();
        for (int p = 0; p < P; p++) {
            planets[p].v = PSUtil::add(planets[p].v, PSUtil::scale(planets[p].a, hd));
        }
    }

    void removeOutOfBounds() {
        for (int i = planets.size(); i >= 0; i--) {
            auto pos = planets[i].x;
            if (abs(pos.first) > BOUNDARY || abs(pos.second) > BOUNDARY) removePlanet(i);
        }
    }

    // Perform one full iteration: compute forces and update. There is no guarantee that the forces array after running this represents the computed forces before running.
    inline void iterate(PSUtil::numeric delta, bool removeOOB = true) {
        computeForces();
        update(delta);

        if (removeOOB) removeOutOfBounds();
    }

    inline std::string sprintf() {
        std::string res = "PlanetSimulator:";
        res += "\nPlanets:\n";
        for (auto p : planets) res += p.to_string() + "\n";
        res += "Stars:\n";
        for (auto s : stars) res += s.to_string() + "\n";
        return res + "/PlanetSimulator";

    }
};


#endif
