#include <bits/stdc++.h>
using namespace std;

#include "src/util.h"
#include "src/planet.h"
#include "src/star.h"
#include "src/simulator.h"

#include "tests.h"

int main() {
    runTests(false);

    PlanetSimulator s1;
    s1.planets.push_back(Planet{1, {0, 0}});
    s1.stars.push_back(Star{100, {100, 0}});
    PlanetSimulator s2;
    s2.planets.push_back(Planet{2, {0, 2}});
    s2.stars.push_back(Star{1000, {1000, 0}});
    PlanetSimulator s3 = s2.join(s1);
    cout << s3.sprintf() << endl;

    s3.rotate(M_PI / 3);
    cout << s3.sprintf() << endl;

    return 0;
}