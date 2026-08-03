#include <bits/stdc++.h>
using namespace std;

#include "util.h"
#include "planet.h"
#include "star.h"
#include "simulator.h"

#include <assert.h>

void runTests(bool verbose = false) {
    cout << "G = " << PSUtil::G << endl;

    // 1. Test basics (planet and simulator information storage and display)

    Planet p;
    cout << p.to_string() << endl;

    Planet p1{100000, {0, 0}};
    Planet p2{100000, {1, 0}};

    auto grav = PSUtil::gravity(p1.mass, p2.mass, 1);
    assert(PSUtil::equals(grav, PSUtil::G * 1e10));

    // 2. Test force computations

    PlanetSimulator s;
    s.planets.push_back(p1);
    s.planets.push_back(p2);

    s.computeForces();
    if (verbose) cout << PSUtil::disp(s.forces) << endl;

    assert(PSUtil::equals(s.forces[0], {PSUtil::G * 1e10, 0}));
    assert(PSUtil::equals(s.forces[1], {-1 * PSUtil::G * 1e10, 0}));

    // 3. Test a simple star system (one star and one planet, computed for circular orbit.)

    PlanetSimulator s2;
    s2.stars.push_back(Star{1e+10, {0, 0}});
    s2.planets.push_back(p2);

    s2.computeForces();
    if (verbose) cout << PSUtil::disp(s2.forces) << endl;

    auto velocity = PSUtil::circular(s2.stars[0].mass, PSUtil::distance(s2.planets[0].x, s2.stars[0].x));
    if (verbose) cout << velocity << endl;
    s2.planets[0].v = {0, velocity};

    assert(PSUtil::equals(s2.planets[0].v, {0, sqrt(PSUtil::G * 1e10)}));

    for (int i = 0; i < 1000; i++) {
        s2.iterate(0.01);
        if (i % 100 == 0) cout << "ITER " << (i + 1) << ":" << s2.planets[0].to_string() << ":" << PSUtil::length(s2.planets[0].x) << endl;
        assert(abs(PSUtil::length(s2.planets[0].x) - 1) < 0.001);
    }

    // 4. Out of bounds

    PlanetSimulator ps3;
    ps3.planets.push_back(Planet{1, {1000, 0}, {1, 0}});
    ps3.planets.push_back(Planet{1, {-999, 0}, {-1, 0}});
    for (int i = 0; i < 400; i++) {
        ps3.iterate(0.1);
        for (auto p : ps3.planets) assert(abs(p.x.first) <= 1024 && abs(p.x.second) <= 1024);
    }
}