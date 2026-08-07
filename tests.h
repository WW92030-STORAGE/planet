#include <bits/stdc++.h>
using namespace std;

#include "src/util.h"
#include "src/planet.h"
#include "src/star.h"
#include "src/simulator.h"
#include "src/presets.h"

#include <assert.h>

/*

Most of these tests do not account for compounding of integration errors. 
For small enough timestamps the tests should withstand thousands of iterations.

*/

// Test basics (planet and simulator information storage and display)
void test_basic(bool verbose = false) {
    cout << "test_basic" << endl;
    Planet p;
    cout << p.to_string() << endl;

    Planet p1{100000, {0, 0}};
    Planet p2{100000, {1, 0}};

    auto grav = PSUtil::gravity(p2.mass, 1);
    assert(PSUtil::equals(grav, PSUtil::G * 1e5));
}

// Test force computations
void test_forces(bool verbose = false) {
    cout << "test_forces" << endl;
    Planet p1{100000, {0, 0}};
    Planet p2{100000, {1, 0}};
    PlanetSimulator s;
    s.planets.push_back(p1);
    s.planets.push_back(p2);

    s.computeForces();
    if (verbose) cout << PSUtil::disp(s.accelerations()) << endl;

    assert(PSUtil::equals(s.accelerations()[0], {PSUtil::G * 1e5, 0}));
    assert(PSUtil::equals(s.accelerations()[1], {-1 * PSUtil::G * 1e5, 0}));
}

// Test circular orbit
void test_circular_orbit(bool verbose = false) {
    cout << "test_circular_orbit" << endl;
    Planet p2{100000, {1, 0}};
    PlanetSimulator s2;
    s2.stars.push_back(Star{1e+10, {0, 0}});
    s2.planets.push_back(p2);

    s2.computeForces();
    if (verbose) cout << PSUtil::disp(s2.accelerations()) << endl;

    auto velocity = PSUtil::circular(s2.stars[0].mass, PSUtil::distance(s2.planets[0].x, s2.stars[0].x));
    if (verbose) cout << velocity << endl;
    s2.planets[0].v = {0, velocity};

    assert(PSUtil::equals(s2.planets[0].v, {0, sqrt(PSUtil::G * 1e10)}));

    for (int i = 0; i < 8192; i++) {
        s2.iterate(0.01);
        if (verbose && (i + 1) % 100 == 0) cout << "ITER " << (i + 1) << ":" << s2.planets[0].to_string() << ":" << PSUtil::length(s2.planets[0].x) << endl;
        assert(abs(PSUtil::length(s2.planets[0].x) - 1) < 0.001);
    }
}

// Test out of bounds
void test_oob(bool verbose = false) {
    cout << "test_oob" << endl;
    PlanetSimulator ps3;
    ps3.planets.push_back(Planet{1, {1000, 0}, {1, 0}});
    ps3.planets.push_back(Planet{1, {-999, 0}, {-1, 0}});
    for (int i = 0; i < 8192; i++) {
        ps3.iterate(0.1);
        for (auto p : ps3.planets) assert(abs(p.x.first) <= 1024 && abs(p.x.second) <= 1024);
    }
}

// Moore's 3 body solution
// https://astronomy.stackexchange.com/questions/50297/initial-state-for-a-3-body-problem-to-create-figure-8-restricted-to-2d
void test_moore(bool verbose = false) {
    cout << "test_moore" << endl;
    Planet p0{1, {0.9700436, -0.24308753}, {0.466203685, 0.43236573}};
    Planet p1{1, PSUtil::scale(p0.x, -1), p0.v};
    Planet p2{1, {0, 0}, PSUtil::scale(p0.v, -2)};

    PlanetSimulator ps4;
    ps4.grav = 1;
    ps4.planets.push_back(p0);
    ps4.planets.push_back(p1);
    ps4.planets.push_back(p2);

    PlanetSimulator ps5 = moore();

    for (int i = 0; i < 8192; i++) {
        ps4.iterate(0.01);
        ps5.iterate(0.01);
        if (verbose && (i + 1) % 100 == 0) {
            cout << "G = 1:         " << PSUtil::to_string(ps4.planets[0].x) << endl;
            cout << "G = PSUtil::G: " << PSUtil::to_string(ps4.planets[0].x) << endl;
        }
        assert(abs(ps4.planets[0].x.first) < 1.1 && abs(ps4.planets[0].x.second) < 0.36);
        assert(PSUtil::equals(ps4.planets[0].x, ps5.planets[0].x));
    }    
}

// Arbitrary ring of planets
void test_ring(bool verbose = false) {
    cout << "test_ring" << endl;
    PlanetSimulator ps = ring(1, 5, 7);
    PlanetSimulator ps2 = ring(0.5, 5, 7);
    if (verbose) cout << ps.planets[0].v.second << " " << ps2.planets[0].v.second << endl;
    for (int i = 0; i < 8192; i++) {
        ps.iterate(0.01);
        ps2.iterate(0.01);
        if (verbose && (i + 1) % 100 == 0) {
            cout << PSUtil::to_string(ps.planets[0].x) << " | " << PSUtil::to_string(PSUtil::rect2polar(ps.planets[0].x)) << endl;
            cout << PSUtil::to_string(ps2.planets[0].x) << " | " << PSUtil::to_string(PSUtil::rect2polar(ps2.planets[0].x)) << endl;
        }

        assert(abs(PSUtil::length(ps.planets[0].x) - 7) < 0.001);
        assert(abs(PSUtil::length(ps2.planets[0].x) - 7) < 0.001);
        assert(PSUtil::equals(ps.planets[0].x, ps2.planets[0].x));
    }    
}

void test_moore_scaled(bool verbose = false) {
    cout << "test_moore_scaled" << endl;
    auto scale = 2.0;
    PlanetSimulator ps5 = moore(0.1, scale);

    for (int i = 0; i < 8192; i++) {
        ps5.iterate(0.01);
        if (verbose && (i + 1) % 100 == 0) {
            cout << PSUtil::to_string(ps5.planets[0].x) << "\n";
        }
        assert(abs(ps5.planets[0].x.first) < 1.1 * scale && abs(ps5.planets[0].x.second) < 0.36 * scale);
    }    
}

// Test circular orbit on a massless planet
void test_circular_orbit_massless(bool verbose = false) {
    cout << "test_circular_orbit_massless" << endl;
    Planet p2{0, {1, 0}};
    PlanetSimulator s2;
    s2.stars.push_back(Star{1e+10, {0, 0}});
    s2.planets.push_back(p2);

    s2.computeForces();
    if (verbose) cout << PSUtil::disp(s2.accelerations()) << endl;

    auto velocity = PSUtil::circular(s2.stars[0].mass, PSUtil::distance(s2.planets[0].x, s2.stars[0].x));
    if (verbose) cout << velocity << endl;
    s2.planets[0].v = {0, velocity};

    assert(PSUtil::equals(s2.planets[0].v, {0, sqrt(PSUtil::G * 1e10)}));

    for (int i = 0; i < 8192; i++) {
        s2.iterate(0.01);
        if (verbose && (i + 1) % 100 == 0) cout << "ITER " << (i + 1) << ":" << s2.planets[0].to_string() << ":" << PSUtil::length(s2.planets[0].x) << "\n";
        assert(abs(PSUtil::length(s2.planets[0].x) - 1) < 0.001);
    }
}

// Test vis-viva
void test_visviva(bool verbose = false) {
    cout << "test_visviva" << endl;
    PlanetSimulator sim = singleplanet(1, 2, 1, 1, 16);
    cout << PSUtil::to_string(sim.planets[0].x) << " | " << PSUtil::to_string(sim.planets[0].v) << endl;

    for (int i = 0; i < 1<<16; i++) {
        sim.iterate(0.0001);
        if (verbose && i % 100 == 0) cout << PSUtil::to_string(sim.planets[0].x) << endl;
        auto x = sim.planets[0].x.first + sqrt(3);
        auto y = sim.planets[0].x.second;
        // cout << x << " " << y << " " << x * x + 4 * y * y << endl;
        assert(abs(x * x + 4 * y * y - 4) < 0.001);
    }
}

// Test test_binary_equal
void test_binary_equal(bool verbose = false) {
    cout << "test_binary_equal" << endl;
    PlanetSimulator sim = binary_equal(1, 2, 1, 1);
    cout << PSUtil::to_string(sim.planets[0].x) << " | " << PSUtil::to_string(sim.planets[0].v) << endl;

    for (int i = 0; i < 1<<16; i++) {
        sim.iterate(0.001);
        if (verbose && i % 100 == 0) cout << PSUtil::to_string(sim.planets[0].x) << endl;
        
        auto x = sim.planets[0].x.first + sqrt(3);
        auto y = sim.planets[0].x.second;
        // cout << x << " " << y << " " << x * x + 4 * y * y << endl;
        assert(abs(x * x + 4 * y * y - 4) < 0.001);

        assert(PSUtil::equals(sim.planets[0].x.first, -1 * sim.planets[1].x.first));
        assert(PSUtil::equals(sim.planets[0].x.second, -1 * sim.planets[1].x.second));
    }
}

void runTests(bool verbose = false) {
    cout << "G = " << PSUtil::G << endl;
    
    test_basic(verbose);
    test_forces(verbose);
    test_circular_orbit(verbose);
    test_oob(verbose);
    test_moore(verbose);
    test_ring(verbose);
    test_moore_scaled(verbose);
    test_circular_orbit_massless(verbose);
    test_binary_equal(false);
}