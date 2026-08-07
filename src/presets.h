#ifndef PLANETSIM_PRESETS_NORMEXT
#define PLANETSIM_PRESETS_NORMEXT

#include "planet.h"
#include "star.h"
#include "simulator.h"
#include "util.h"

#include <cmath>

/*

Each of these presets should be invariant under the value of G.
Most of these have relevant computations in comments.

*/

#define GRAV_PARAM PSUtil::numeric grav = PSUtil::G

// Star and planet in a circular orbit. 
PlanetSimulator singleplanet_circular(GRAV_PARAM, PSUtil::numeric R = 1, PSUtil::numeric mp = 1, PSUtil::numeric ms = 16) {
    // In this simplest case, we use the centripetal velocity solver in PSUtil to compute the velocity.
    if (R <= 0) R = 1;
    PlanetSimulator sim{grav};
    sim.stars.push_back(Star{ms, {0, 0}});
    PSUtil::numeric velocity = PSUtil::circular(ms, R, grav);
    sim.planets.push_back(Planet{mp, {R, 0}, {0, velocity}});
    return sim;
}

// Star and planet. The planet orbits around in an elliptical manner.
PlanetSimulator singleplanet(GRAV_PARAM, PSUtil::numeric a = 1, PSUtil::numeric b = 0.5, PSUtil::numeric mp = 1, PSUtil::numeric ms = 16) {
    if (b > a) std::swap(a, b);
    PSUtil::numeric c = sqrt(a * a - b * b);

    PlanetSimulator sim{grav};
    sim.planets.push_back(Planet{mp, {a - c, 0}});
    sim.stars.push_back(Star{ms, {0, 0}});

    /*

    E = U + K = const (since gravity, like other inverse square fields, is conservative)
    U + K = -GM(mp)/r + 0.5(mp)v^2 = const (where G is the gravitational constant.)

    where r is the distance between the bodies and a is the semimajor axis and M = ms. Dividing by (mp) on both sides:

    -GM/r + 0.5v^2 = const. Analyzing at two points: the perihelion and the aphelion:
    -GM/r + 0.5v^2 = -GM/r' + 0.5v'^2 (where r, v are taken at perihelion and r', v' at the aphelion. Naturally, r' > r.)
    0.5v^2 - 0.5v'^2 = GM/r - GM/r'
    v^2 - v'^2 = 2GM[1/r - 1/r']

    Conservation of angular momentum: rv = r'v'. Thus v' = rv/r'.
    v^2 - [rv/r']^2 = v^2[1 - [r/r']^2] = v^2[r'^2 - r^2] / [r']^2
    and 2GM[1/r - 1/r'] = 2GM[r' - r] / [r * r']

    Thus v^2[r'^2 - r^2] / r'^2 = 2GM[r' - r] / [r * r']
    v^2[r' + r][r' - r] / r'^2 = 2GM[r' - r] / [r * r']
    v^2[r' + r] / r' = 2GM / r. Since 2a = r' + r we obtain:
    v^2[2a/r'] = 2GM/r
    v^2 = 2GMr' / 2ar. Since r' = 2a - r:
    v^2 = 2GM[2a - r] / 2ar = 2GM[1/r - 1/2a] = GM[2/r - 1/a]
    
    The end result is the VIS-VIVA EQUATION which gives us the perihelion velocity:

    v^2 = GM[2/r - 1/a]
    
    *We derived the result for ellipses, where 1/a > 0. 1/a = 0 in parabolas and 1/a < 0 in hyperbolas.
    
    */

    PSUtil::numeric visviva = (grav * ms) * (2.0 / (a - c) - 1.0 / a);
    sim.planets[0].v = {0, sqrt(visviva)};

    return sim;
}


// Binary star system (actually a binary planet system) with two identical masses and orbiting on ellipses of 1/2major axis (axis) and 1/2minor axis (scale)
PlanetSimulator binary_equal(GRAV_PARAM, PSUtil::numeric a = 1, PSUtil::numeric b = 1, PSUtil::numeric mass = 1) {
    if (b > a) std::swap(a, b);
    PSUtil::numeric c = sqrt(a * a - b * b);

    PlanetSimulator sim{grav};
    sim.planets.push_back(Planet{mass, {a - c, 0}});
    sim.planets.push_back(Planet{mass, {c - a, 0}});

    /*
    
    Alright, let's get this party started. To solve this we will use a thing we call REDUCED MASS. 
    REDUCED MASS performs a coordinate shift from absolute coordinates to relative (irrotational) coordinates around one of the masses in the binary system:

    X --> X'
    m1' = m1 + m2
    m2' = m1 * m2 / m1'
    x1' = (0, 0)
    x2' = x2 - x1 <--- This is the relative position of mass 2 with respect to mass 1

    Now we're in the general one-star elliptical system. The VIS-VIVA equation is used here:

    v2' = sqrt[G * m1' * [2/r - 1/a']]
    what is r? r = (a' - c')

    But remember that v2' is relative velocity. At the perihelions, where the planets are closest to the common barycenter, their tangential velocities are opposing each other.
    Thus v2 = -v1 = 0.5v2'.

    BUT WAIT!!! It's not over yet!

    Remember that everything we just did was in terms of the relative orbit. To ascertain values from the absolute orbits we must do some scaling.

    At the perihelion, the two planets are 2 * (a - c) apart. At the aphelion, thte planets are 2 * (a + c) apart.
    Thus a' = 2a, the sum of these values.

    Similarly, when the planets are at the minor axis points, when planet 2 is on a minor axis point relative to planet 1, the vertical separation is 2b. Thus b' = 2b.
    Thus c' = sqrt(a'^2 - b'^2) = 2c.
    
    */
    
    PSUtil::numeric aprime = 2 * a;
    PSUtil::numeric bprime = 2 * b;
    PSUtil::numeric cprime = sqrt(aprime * aprime - bprime * bprime);
    PSUtil::numeric vprime = 0.5 * sqrt(2 * grav * mass * (2.0 / (aprime - cprime) - 1.0 / aprime));
    sim.planets[0].v = {0, vprime};
    sim.planets[1].v = {0, -vprime};

    return sim;
}

// Circle of planets
PlanetSimulator ring(GRAV_PARAM, int N = 3, PSUtil::numeric R = 1, PSUtil::numeric mass = -1) {
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

// Moore figure 8 solution. See the README for the initial state.
PlanetSimulator moore(GRAV_PARAM, PSUtil::numeric scale = 1) {
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