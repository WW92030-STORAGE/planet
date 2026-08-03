#ifndef PLANETSIM_UTIL_NORMEXT
#define PLANETSIM_UTIL_NORMEXT

#include <utility>
#include <algorithm>

namespace PSUtil {

typedef double numeric;
typedef std::pair<numeric, numeric> point;

constexpr numeric G = 6.67387E-11;
constexpr numeric EPSILON = 1e-5;

/*

Utility functions and constants.

Newton's Law: F = G(m1 * m2) / R^2

*/

// COMPARE FLOATS AND POINTS

inline bool equals(numeric a, numeric b) {
    return abs(a - b) < EPSILON;
}

inline bool equals(point a, point b) {
    return equals(a.first, b.first) && equals(a.second, b.second);
}

inline bool zero(numeric n) {
    return abs(n) < EPSILON;
}

// COMPUTE GRAVITATIONAL FORCE

inline numeric gravity(numeric m1, numeric m2, numeric r) {
    if (abs(m1) < abs(m2)) std::swap(m1, m2);
    return ((G * m1) / (r * r)) * m2;
}

// OPERATIONS ON POINTS

inline numeric lengthSquared(point p) {
    return p.first * p.first + p.second * p.second;
}

inline numeric length(point p) {
    return sqrt(lengthSquared(p));
}

inline point normalize(point p) {
    numeric len = length(p);
    if (equals(len, 0)) return p;
    numeric recip = 1.0 / len;
    return {p.first * recip, p.second * recip};
}

inline point add(point a, point b) {
    return {a.first + b.first, a.second + b.second};
}

inline point sub(point a, point b) {
    return {a.first - b.first, a.second - b.second};
}

inline point scale(point p, numeric x) {
    return {p.first * x, p.second * x};
}

inline numeric distanceSquared(point a, point b) {
    return lengthSquared(sub(a, b));
}

inline numeric distance(point a, point b) {
    return length(sub(a, b));
}

// PRINT STUFF AS STRING

std::string to_string(point p) {
    return "(" + std::to_string(p.first) + ", " + std::to_string(p.second) + ")";
}

inline std::string disp(std::vector<point> v) {
    std::string res = "[";
    for (int i = 0; i < v.size(); i++) {
        if (i) res = res + ", ";
        res = res + to_string(v[i]);
    }
    return res + "]";
}

// UTILITY FUNCTIONS

// Compute tangential velocity magnitude for a circular orbit
constexpr PSUtil::numeric circular(PSUtil::numeric mstar, PSUtil::numeric r) {
    /*
        
    a = v^2 / r
    v^2 = ar
    v = sqrt(ar)
    a = G(mstar) / (r * r)
    thus v = sqrt(G * mstar / r)

    */

    return sqrt(G * mstar / r);
}

}
#endif