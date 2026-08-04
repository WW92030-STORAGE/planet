#ifndef PLANETSIM_SIMULATOR_NORMEXT
#define PLANETSIM_SIMULATOR_NORMEXT

#include "planet.h"
#include "star.h"
#include "util.h"

#include <vector>
#include <string>

struct PlanetSimulator {
    PSUtil::numeric grav = PSUtil::G;
    PSUtil::numeric BOUNDARY = 1024;
    
    std::vector<Planet> planets;
    std::vector<Star> stars;

    void removePlanet(int index) {
        if (index >= planets.size() || index < 0) return;
        std::swap(planets[index], planets[planets.size() - 1]);
        planets.pop_back();
    }

    void removeStar(int index) {
        if (index >= stars.size() || index < 0) return;
        std::swap(stars[index], stars[stars.size() - 1]);
        stars.pop_back();
    }

    // actually computes accelerations
    void computeForces() {
        int P = planets.size();
        int S = stars.size();

        for (int p = 0; p < P; p++) planets[p].a = {0, 0};

        // compute. right now it's O(P * P + P * S)

        // compute planet-star forces

        for (int p = 0; p < P; p++) {
            for (int s = 0; s < S; s++) {
                Planet planet = planets[p];
                Star star = stars[s];

                PSUtil::point displacement = PSUtil::sub(star.x, planet.x);
                PSUtil::numeric len = PSUtil::length(displacement);
                if (PSUtil::zero(len)) continue;
                PSUtil::numeric force_magnitude = PSUtil::gravity(star.mass, len, grav);

                // unit vector + add
                planets[p].a = PSUtil::axpy(force_magnitude / len, displacement, planets[p].a);
            }
        }

        // compute planet-planet forces. Here, we use "star" to denote what is influencing the planet.

        for (int p = 0; p < P; p++) {
            for (int s = 0; s < P; s++) {
                if (s == p) continue;
                Planet planet = planets[p];
                Planet star = planets[s];

                PSUtil::point displacement = PSUtil::sub(star.x, planet.x);
                PSUtil::numeric len = PSUtil::length(displacement);
                if (PSUtil::zero(len)) continue;
                PSUtil::numeric force_magnitude = PSUtil::gravity(star.mass, len, grav);

                // unit vector + add
                planets[p].a = PSUtil::axpy(force_magnitude / len, displacement, planets[p].a);

                // std::cout << p << " " << s << ":" << force_magnitude << " | " << PSUtil::to_string(displacement) << " | " << PSUtil::to_string(forces[p]) << "\n";
            }
        }
    }

    // assumes you've computed forces beforehand. otherwise use the iterate method below. There is no guarantee that the forces array is preserved after this.
    void update(PSUtil::numeric delta) {
        updateLeapfrog(delta);
    }

    void updateEuler(PSUtil::numeric delta) {
        int P = planets.size();
        for (int p = 0; p < P; p++) {
            planets[p].v = PSUtil::axpy(delta, planets[p].a, planets[p].v);
            planets[p].x = PSUtil::axpy(delta, planets[p].v, planets[p].x);
        }
    }

    void updateBetter(PSUtil::numeric delta) {
        auto hd2 = delta * delta * 0.5;
        int P = planets.size();

        for (int p = 0; p < P; p++) {
            auto displacement = PSUtil::add(PSUtil::scale(planets[p].v, delta), PSUtil::scale(planets[p].a, hd2));
            planets[p].x = PSUtil::add(planets[p].x, displacement);
            planets[p].v = PSUtil::axpy(delta, planets[p].a, planets[p].v);
        }
    }

    // WARNING: This integrator will compute forces again during the step.
    void updateLeapfrog(PSUtil::numeric delta) {
        auto hd = delta * 0.5;
        int P = planets.size();
        for (int p = 0; p < P; p++) {
            planets[p].v = PSUtil::axpy(hd, planets[p].a, planets[p].v);
            planets[p].x = PSUtil::axpy(delta, planets[p].v, planets[p].x);
        }
        computeForces();
        for (int p = 0; p < P; p++) {
            planets[p].v = PSUtil::axpy(hd, planets[p].a, planets[p].v);
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
        res += "G=" + std::to_string(grav) + "\n";
        return res + "/PlanetSimulator";
    }

    std::vector<PSUtil::point> accelerations() {
        std::vector<PSUtil::point> res(planets.size());
        for (int p = 0; p < planets.size(); p++) res[p] = planets[p].a;
        return res;
    }
};


#endif
