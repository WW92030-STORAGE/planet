#include <bits/stdc++.h>
using namespace std;

#include "util.h"
#include "planet.h"
#include "star.h"
#include "simulator.h"

int main() {
    cout << "G = " << PSUtil::G << endl;

    Planet p;
    cout << p.to_string() << endl;

    Planet p1{100000, {0, 0}};
    Planet p2{100000, {1, 0}};

    cout << PSUtil::gravity(100000, 100000, 1) << endl;

    PlanetSimulator s;
    s.planets.push_back(p1);
    s.planets.push_back(p2);

    s.computeForces();
    cout << PSUtil::disp(s.forces) << endl;

    PlanetSimulator s2;
    s2.stars.push_back(Star{1e+10, {0, 0}});
    s2.planets.push_back(p2);

    s2.computeForces();
    cout << PSUtil::disp(s2.forces) << endl;

    auto velocity = PSUtil::circular(s2.stars[0].mass, PSUtil::distance(s2.planets[0].x, s2.stars[0].x));
    cout << velocity << endl;
    s2.planets[0].v = {0, velocity};

    for (int i = 0; i < 1000; i++) {
        cout << "ITER " << (i + 1) << ":";
        s2.iterate(0.01);
        cout << s2.planets[0].to_string() << ":" << PSUtil::length(s2.planets[0].x) << endl;
    }

    return 0;
}