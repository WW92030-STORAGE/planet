#ifndef PLANETSIM_PLANET_NORMEXT
#define PLANETSIM_PLANET_NORMEXT

#include "util.h"
#include <string>
#include <algorithm>

struct Planet {
    PSUtil::numeric mass = 1;
    PSUtil::point x = {0, 0};
    PSUtil::point v = {0, 0};
    PSUtil::point a = {0, 0};
    std::string name = "Planet";

    bool operator<(const Planet& other) {
        if (!PSUtil::equals(mass, other.mass)) return mass < other.mass;
        if (!PSUtil::equals(x, other.x)) return x < other.x;
        if (!PSUtil::equals(v, other.v)) return v < other.v;
        if (!PSUtil::equals(a, other.a)) return a < other.a;
        if (name != other.name) return name < other.name;
        return false;
    }

    bool operator==(const Planet& other) {
        if (!PSUtil::equals(mass, other.mass)) return false;
        if (!PSUtil::equals(x, other.x)) return false;
        if (!PSUtil::equals(v, other.v)) return false;
        if (!PSUtil::equals(a, other.a)) return false;
        if (name != other.name) return false;
        return true;
    }

    std::string to_string() {
        return "Planet[mass=" + std::to_string(mass) + ", x=" + PSUtil::to_string(x) + ", v=" + PSUtil::to_string(v) + ", a=" + PSUtil::to_string(a) + ", name=" + name + "]";
    }
};

#endif