#ifndef PLANETSIM_STAR_NORMEXT
#define PLANETSIM_STAR_NORMEXT

// A star is defined as a body that is not influenced by external forces.

#include "util.h"
#include <string>

struct Star {
    PSUtil::numeric mass = 1;
    PSUtil::point x = {0, 0};
    std::string name = "Star";

    bool operator<(const Star& other) {
        if (!PSUtil::equals(mass, other.mass)) return mass < other.mass;
        if (!PSUtil::equals(x, other.x)) return x < other.x;
        if (name != other.name) return name < other.name;
        return false;
    }

    bool operator==(const Star& other) {
        if (!PSUtil::equals(mass, other.mass)) return false;
        if (!PSUtil::equals(x, other.x)) return false; 
        if (name != other.name) return false;
        return true;
    }

    std::string to_string() {
        return "Star[mass=" + std::to_string(mass) + ", x=" + PSUtil::to_string(x) + ", name=" + name + "]";
    }
};

#endif