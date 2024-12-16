#ifndef __raytracer__
#define __raytracer__

#include "ppm.h"
#include "parser.h"
#include <iostream>
#include <limits>
#include <cmath>
#include <algorithm>
#include <pthread.h>

using namespace parser;

//Useful structures for rendering the image.
typedef struct Ray{
    Vec3f ray_origin;
    Vec3f ray_direction;
    bool ray_is_shadow;
    int ray_depth;
} ray;

enum Object_Type{
	MESH, TRIANGLE, SPHERE
};


typedef struct Hit{
	Object_Type object_type;
	int material_id;
	int object_index;
	float time;
	Vec3f intersectionPoint;
	Vec3f surfaceNormal;
    int face_index;
    bool hitHappened;
    bool in_shadow;
}hit;

typedef struct ThreadStruct{
    int start_row;
    int end_row;
    int width;
    int height;
    unsigned char *image;
    Scene *scene;
    int camera_index;
}thread_struct;



//Operator overloaders for 3D vectors.
Vec3f operator+(const Vec3f &vector1, const Vec3f &vector2){
    Vec3f solution;
    solution.x = vector1.x + vector2.x;
    solution.y = vector1.y + vector2.y;
    solution.z = vector1.z + vector2.z;

    return solution;
}

Vec3f operator-(const Vec3f &vector1){
    Vec3f solution;
    solution.x = - vector1.x;
    solution.y = - vector1.y;
    solution.z = - vector1.z;

    return solution;
}

Vec3f operator-(const Vec3f &vector1, const Vec3f &vector2){
    Vec3f solution;
    solution.x = vector1.x - vector2.x;
    solution.y = vector1.y - vector2.y;
    solution.z = vector1.z - vector2.z;

    return solution;
}

Vec3f operator*(const Vec3f &vector1, const Vec3f &vector2){
    Vec3f solution;
    solution.x = vector1.x * vector2.x;
    solution.y = vector1.y * vector2.y;
    solution.z = vector1.z * vector2.z;

    return solution;
}

Vec3f operator*(const Vec3f &vector1, const float multiplier){
    Vec3f solution;
    solution.x = vector1.x * multiplier;
    solution.y = vector1.y * multiplier;
    solution.z = vector1.z * multiplier;

    return solution;
}


Vec3f operator/(const Vec3f &vector1, const float divisor){
    Vec3f solution;
    solution.x = vector1.x / divisor;
    solution.y = vector1.y / divisor;
    solution.z = vector1.z / divisor;

    return solution;
}

inline void assign_hit(Hit &hit1, const Hit &hit2){
    hit1.object_type = hit2.object_type;
	hit1.material_id = hit2.material_id;
	hit1.object_index = hit2.object_index;
	hit1.time = hit2.time;
	hit1.intersectionPoint = hit2.intersectionPoint;
	hit1.surfaceNormal = hit2.surfaceNormal;
    hit1.face_index = hit2.face_index;
    hit1.hitHappened = hit2.hitHappened;
}

unsigned char clamp_function(float value){
    if (value > 255.0f) return 255;
    if(value < 0) return 0;
    return static_cast<unsigned char>(round(value));
}



//Other functions for calculation with 3D vectors.

inline float find_length(const Vec3f &vector){
    return sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

inline float find_distance(const Vec3f &vector1, const Vec3f &vector2){
    Vec3f temp_result = vector1 - vector2;
    return sqrt(pow(temp_result.x, 2) + pow(temp_result.y, 2) + pow(temp_result.z, 2));
}

inline float dot_product_calculator(const Vec3f &vector1, const Vec3f &vector2){
	return ((vector1.x * vector2.x) + (vector1.y * vector2.y) + (vector1.z * vector2.z));
}

inline Vec3f cross_product_calculator(const Vec3f &vector1, const Vec3f &vector2){
    Vec3f result;
    result.x = (vector1.y * vector2.z) - (vector1.z * vector2.y);
    result.y = (vector1.z * vector2.x) - (vector1.x * vector2.z);
    result.z = (vector1.x * vector2.y) - (vector1.y * vector2.x);

    return result;
}

inline Vec3f normalise_vector(Vec3f vector){
    return vector / find_length(vector);
}


//Function for determinant calculation.
inline float determinant_calculator(Vec3f vector1, Vec3f vector2, Vec3f vector3){
    
    return (vector1.x * ((vector2.y * vector3.z) - (vector2.z * vector3.y))) -
           (vector2.x * ((vector1.y * vector3.z) - (vector1.z * vector3.y))) +
           (vector3.x * ((vector1.y * vector2.z) - (vector1.z * vector2.y)));

}



#endif // __raytracer__
