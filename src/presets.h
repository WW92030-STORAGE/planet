#ifndef PLANETSIM_PRESETS_NORMEXT
#define PLANETSIM_PRESETS_NORMEXT

#include "planet.h"
#include "star.h"
#include "simulator.h"
#include "util.h"

#include <cmath>

/*

Each of these presets should be invariant under the value of G.

*/

// Star and planet
PlanetSimulator singleplanet(PSUtil::numeric grav = PSUtil::G, PSUtil::numeric R = 1, PSUtil::numeric mp = 1, PSUtil::numeric ms = 16) {
    if (R <= 0) R = 1;
    PlanetSimulator sim{grav};
    sim.stars.push_back(Star{ms, {0, 0}});
    PSUtil::numeric velocity = PSUtil::circular(ms, R, grav);
    sim.planets.push_back(Planet{mp, {R, 0}, {0, velocity}});
    return sim;
}

// Circle of planets
PlanetSimulator ring(PSUtil::numeric grav = PSUtil::G, int N = 3, PSUtil::numeric R = 1, PSUtil::numeric mass = -1) {
    if (N <= 0) N = 3;
    if (mass <= 0) mass = 1.0 / grav;

    /*

    In a coordinate frame where the center is (0, 0) and the victim (i.e. any planet) is (1, 0)
    
    sum of accelerations = v^2 / R
    v = sqrt(sum of accelerations * R)

    what is the sum of the accelerations?
    let theta = 2pi/N
    then we are summing from i : 1 to N - 1:
        grav * (mass / distance^2) * cos(eccentric angle)
    = (grav * mass) * (sum of 1 / dist^2 * cos(phi))

    phi is the angle at the point (1, 0) subtended by the origin and the planet exerting the force (some other planet).

    what is dist^2 for each value of i?
    Form an isosceles triangle of legs R and apex angle theta' = theta * i, the distance is the third side, and phi = 90 - (theta' / 2).
    The angle bisector from the apex splits the distance in two. 
    thus 0.5 *distance = R * sin(theta' / 2)
    distance^2 = 4R^2 * sin^2(theta' / 2)

    so the sum is (grav * mass) * (1.0 / 4R^2) * sum[csc^2(theta' / 2) * cos(phi) = csc(theta' / 2)]

    v = sqrt(a * R) = sqrt(grav * mass * 1/4R * ...)
    */

    PSUtil::numeric half_theta = M_PI / N;

    PSUtil::numeric velocity = 0;
    for (int i = 1; i < N; i++) {
        PSUtil::numeric value = sin(half_theta * i);
        velocity += 1.0 / (value);
    }

    velocity = sqrt(velocity * grav * mass * (0.25 / R));

    PlanetSimulator ps;
    ps.grav = grav;
    PSUtil::numeric theta = (2 * M_PI) / N;
    for (int i = 0; i < N; i++) {
        PSUtil::numeric arg = theta * i;
        Planet p{mass, {R * cos(arg), R * sin(arg)}, {-1 * sin(arg) * velocity, cos(arg) * velocity}};
        ps.planets.push_back(p);
    }
    return ps;
}

// Moore figure 8 solution
PlanetSimulator moore(PSUtil::numeric grav = PSUtil::G, PSUtil::numeric scale = 1) {
    auto cubed = scale * scale * scale; // The force (acceleration) and distances scale up by scale. This contributes scale^3 to the mass.
    Planet p0{cubed / grav, {0.9700436 * scale, -0.24308753 * scale}, {0.466203685 * scale, 0.43236573 * scale}};
    Planet p1{cubed / grav, PSUtil::scale(p0.x, -1), p0.v};
    Planet p2{cubed / grav, {0, 0}, PSUtil::scale(p0.v, -2)};

    PlanetSimulator ps4;
    ps4.grav = grav;
    ps4.planets.push_back(p0);
    ps4.planets.push_back(p1);
    ps4.planets.push_back(p2);
    return ps4;
}

#endif