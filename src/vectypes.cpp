#include "vectypes.h"

point3f convert_to_point3f(const vec3f& p) {
    point3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}

point3f convert_to_point3f(const normal3f& p) {
    point3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}

vec3f convert_to_vec3f(const point3f& p) {
    vec3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}

vec3f convert_to_vec3f(const normal3f& p) {
    vec3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}

normal3f convert_to_normal3f(const vec3f& p) {
    normal3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}

normal3f convert_to_normal3f(const point3f& p) {
    normal3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}
