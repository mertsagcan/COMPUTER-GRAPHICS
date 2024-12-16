#include "raytracer.h"


using namespace std;

#define EPSILON (numeric_limits<float>::epsilon())

typedef unsigned char RGB[3];

Vec3f pixel_color_calculator(const Ray &ray,const Scene &scene);


Ray ray_generator(float i, float j, const Camera &camera){

    float plane_left = camera.near_plane.x;
    float plane_right = camera.near_plane.y;
    float plane_bottom = camera.near_plane.z;
    float plane_top = camera.near_plane.w;

    //Find v. It is the up vector.
    Vec3f v = normalise_vector(camera.up);

    //Normal gaze.
    Vec3f normal_gaze = normalise_vector(camera.gaze);

    // w is the opposite of gaze vector.
    Vec3f w = normal_gaze * -1;

    //Find u. It is the cross product of v and w.
    Vec3f u = normalise_vector(cross_product_calculator(v, w));

    //Find m.
    Vec3f m = camera.position + (normal_gaze * camera.near_distance);

    //Find q.
    Vec3f q = m + (u * plane_left) + (v * plane_top);

    //Calculate su and sv.
    float su = (i + 0.5) * (plane_right - plane_left) / camera.image_width;
    float sv = (j + 0.5) * (plane_top - plane_bottom) / camera.image_height;

    //Calculate s. This is the target point.
    Vec3f s = q + u * su - v * sv;

    //Subtract origin from target point and normalise.This will give the direction.
    Vec3f direction = normalise_vector(s - camera.position);

    Ray ray;
    ray.ray_direction = direction;
    ray.ray_origin = camera.position;
    ray.ray_is_shadow = false;  //Since this is a pirmary ray.
    ray.ray_depth = 0;          //Since this is a primary ray.

    return ray;
}

Hit find_sphere_hit(const Ray &ray, const Scene &scene){
	Hit sphere_hit;
	sphere_hit.hitHappened = false;
	float t;

	int sphere_num = scene.spheres.size();
	int min_time_sphere_index = -1;
	float min_time = INFINITY;

	for(int current_sphere_index = 0; current_sphere_index < sphere_num; current_sphere_index++){
		Sphere sphere = scene.spheres[current_sphere_index];

		Vec3f o_c = ray.ray_origin - scene.vertex_data[sphere.center_vertex_id - 1];
		float A = dot_product_calculator(ray.ray_direction, ray.ray_direction);
		float B = 2 * dot_product_calculator(ray.ray_direction, o_c);
		float C = dot_product_calculator(o_c, o_c) - (sphere.radius * sphere.radius);

		float discriminant = B * B - 4 * A * C;

		if(discriminant < -EPSILON){
			continue;
		}

		else if(discriminant > -EPSILON && discriminant < EPSILON){ 
			t = -B / dot_product_calculator(ray.ray_direction,ray.ray_direction);

			if(t > EPSILON && t < min_time){
				min_time = t;
				min_time_sphere_index = current_sphere_index;
				sphere_hit.hitHappened = 1;
			}                   
		}

		else{
			float t1, t2;
			t1 = (-1 * B + sqrtf(discriminant))
							/ (A * 2);
			t2 = (-1 * B - sqrtf(discriminant))
							/ (A * 2);

            if(t1 > EPSILON && t2 <= EPSILON){
                t = t1;
            }
            else if(t2 > EPSILON && t1 <= EPSILON){
                t = t2;
            }
            else if(t1 > EPSILON && t2 > EPSILON){
                t = min(t1, t2);
            }
            else{
                continue;
            }

			if(t > EPSILON && t < min_time){
				min_time = t;
				min_time_sphere_index = current_sphere_index;
				sphere_hit.hitHappened = 1;
			}
		}
	}
	
	if(min_time_sphere_index != -1){
		sphere_hit.object_type = SPHERE;
		sphere_hit.material_id = scene.spheres[min_time_sphere_index].material_id;
		sphere_hit.object_index = min_time_sphere_index;
		sphere_hit.time = min_time;
        sphere_hit.face_index = -1;
		sphere_hit.intersectionPoint = ray.ray_origin + ray.ray_direction * min_time;
		sphere_hit.surfaceNormal = normalise_vector((sphere_hit.intersectionPoint - scene.vertex_data[scene.spheres[min_time_sphere_index].center_vertex_id - 1]) / scene.spheres[min_time_sphere_index].radius);
	}
	
	return sphere_hit;

}

Hit find_triangle_hit(const Ray &ray, const Scene &scene){
	Hit triangle_hit;
	triangle_hit.hitHappened = false;
	float t;

	int triangle_num = scene.triangles.size();
	int min_time_triangle_index = -1;
	float min_time = INFINITY;

	for(int current_triangle_index = 0; current_triangle_index < triangle_num; current_triangle_index++){
		Triangle triangle = scene.triangles[current_triangle_index];

		Vec3f a = scene.vertex_data[triangle.indices.v0_id - 1];
        Vec3f b = scene.vertex_data[triangle.indices.v1_id - 1];
        Vec3f c = scene.vertex_data[triangle.indices.v2_id - 1];


        
		float A = determinant_calculator(a - b, a - c, ray.ray_direction);
		float beta = determinant_calculator(a - ray.ray_origin, a - c, ray.ray_direction) / A;
		float gama = determinant_calculator(a - b, a - ray.ray_origin, ray.ray_direction) / A;
		t = determinant_calculator(a - b, a - c, a - ray.ray_origin) / A;
        


		if(beta + gama > 1 + EPSILON || beta < -EPSILON || gama < -EPSILON){
			continue;
		}

		if(t > EPSILON && t < min_time){
			min_time = t;
			min_time_triangle_index = current_triangle_index;
			triangle_hit.hitHappened = 1;
		}
	}

	if(min_time_triangle_index != -1){
        Vec3f a1 = scene.vertex_data[scene.triangles[min_time_triangle_index].indices.v0_id - 1];
        Vec3f b1 = scene.vertex_data[scene.triangles[min_time_triangle_index].indices.v1_id - 1];
        Vec3f c1 = scene.vertex_data[scene.triangles[min_time_triangle_index].indices.v2_id - 1];
		triangle_hit.object_type = TRIANGLE;
		triangle_hit.material_id = scene.triangles[min_time_triangle_index].material_id;
		triangle_hit.object_index = min_time_triangle_index;
		triangle_hit.time = min_time;
        triangle_hit.face_index = -1;
		triangle_hit.intersectionPoint = ray.ray_origin + ray.ray_direction * min_time;
		triangle_hit.surfaceNormal = normalise_vector(cross_product_calculator(b1 - a1, c1 - a1));
	}

	return triangle_hit;
}

Hit find_mesh_hit(const Ray &ray, const Scene &scene){
	Hit mesh_hit;
	mesh_hit.hitHappened = false;
	float t;
	int mesh_num = scene.meshes.size();
	int min_time_triangle_index_i = -1;
    int min_time_triangle_index_j = -1;
	float min_time = INFINITY;

	for(int i = 0; i < mesh_num; i++){
        int cur_mesh_size = scene.meshes[i].faces.size();
        for(int j = 0; j < cur_mesh_size; j++){

            Vec3f a = scene.vertex_data[scene.meshes[i].faces[j].v0_id - 1];
            Vec3f b = scene.vertex_data[scene.meshes[i].faces[j].v1_id - 1];
            Vec3f c = scene.vertex_data[scene.meshes[i].faces[j].v2_id - 1];

            
            float A = determinant_calculator(a - b, a - c, ray.ray_direction);
            float beta = determinant_calculator(a - ray.ray_origin, a - c, ray.ray_direction) / A;
            float gama = determinant_calculator(a - b, a - ray.ray_origin, ray.ray_direction) / A;
            t = determinant_calculator(a - b, a - c, a - ray.ray_origin) / A;

            if(beta + gama > 1 + EPSILON || beta < -EPSILON || gama < -EPSILON){
                continue;
            }

            if(t > EPSILON && t < min_time){
                min_time = t;
                min_time_triangle_index_i = i;
                min_time_triangle_index_j = j;
                mesh_hit.hitHappened = 1;
            }
        }
		
	}

	if(min_time_triangle_index_i != -1){
        Vec3f a1 = scene.vertex_data[scene.meshes[min_time_triangle_index_i].faces[min_time_triangle_index_j].v0_id - 1];
        Vec3f b1 = scene.vertex_data[scene.meshes[min_time_triangle_index_i].faces[min_time_triangle_index_j].v1_id - 1];
        Vec3f c1 = scene.vertex_data[scene.meshes[min_time_triangle_index_i].faces[min_time_triangle_index_j].v2_id - 1];
		mesh_hit.object_type = MESH;
		mesh_hit.material_id = scene.meshes[min_time_triangle_index_i].material_id;
		mesh_hit.object_index = min_time_triangle_index_i;
        mesh_hit.face_index = min_time_triangle_index_j;
		mesh_hit.time = min_time;
		mesh_hit.intersectionPoint = ray.ray_origin + ray.ray_direction * min_time;
		mesh_hit.surfaceNormal = normalise_vector(cross_product_calculator(b1 - a1, c1 - a1));
	}

	return mesh_hit;
}


Hit find_the_closest_hit(const Ray &ray, const Scene &scene){
    //Trace all shapes and find the closest intersection if there is one.
    Hit closest_hit;
    closest_hit.hitHappened = false;
    float min_time = INFINITY;
    Object_Type type;

    Hit sphere_hit = find_sphere_hit(ray, scene);
    Hit triangle_hit = find_triangle_hit(ray, scene);
    Hit mesh_hit = find_mesh_hit(ray, scene);

    if(sphere_hit.hitHappened){
        if(sphere_hit.time < min_time){
            min_time = sphere_hit.time;
            assign_hit(closest_hit, sphere_hit);
        }
    }
    if(triangle_hit.hitHappened){
        if(triangle_hit.time < min_time){
            min_time = triangle_hit.time;
            assign_hit(closest_hit, triangle_hit);
        }
    }
    if(mesh_hit.hitHappened){
        if(mesh_hit.time < min_time){
            min_time = mesh_hit.time;
            assign_hit(closest_hit, mesh_hit);
        }
    }

    return closest_hit;
}


Vec3f difuse_color_calculator(const Ray &ray, Hit hit, const Scene &scene, int current_light){
    Vec3f difuse_color = {0,0,0};
    Vec3f light_direction = normalise_vector(scene.point_lights[current_light].position - hit.intersectionPoint);
    float dot_product = max(0.0f, dot_product_calculator(light_direction, hit.surfaceNormal));
    
    difuse_color = difuse_color + scene.materials[hit.material_id - 1].diffuse * (scene.point_lights[current_light].intensity / pow(find_distance(scene.point_lights[current_light].position, hit.intersectionPoint) ,2) * dot_product);

    return difuse_color;
}

Vec3f specular_color_calculator(const Ray &ray, Hit hit, const Scene &scene, int current_light){
    Vec3f specular_color = {0,0,0};
    Vec3f light_direction = normalise_vector(scene.point_lights[current_light].position - hit.intersectionPoint) - ray.ray_direction;
    float cos_p = dot_product_calculator(hit.surfaceNormal, normalise_vector(light_direction));
    float pow_cos_p = pow(cos_p, scene.materials[hit.material_id - 1].phong_exponent);
    if(pow_cos_p > EPSILON){
        specular_color = specular_color + scene.materials[hit.material_id - 1].specular * (scene.point_lights[current_light].intensity / pow(find_distance(scene.point_lights[current_light].position, hit.intersectionPoint) ,2) * pow_cos_p);
    }
    return specular_color;
}

Ray reflection_ray_calculator(Ray ray, Hit hit, const Scene &scene){
    //Reflect the ray.
    Ray reflection_ray;

    reflection_ray.ray_origin = hit.intersectionPoint + hit.surfaceNormal * scene.shadow_ray_epsilon;
    reflection_ray.ray_direction = normalise_vector(hit.surfaceNormal * 2 * dot_product_calculator(-ray.ray_direction, hit.surfaceNormal) + ray.ray_direction);
    reflection_ray.ray_is_shadow = false;
    return reflection_ray;
}

Vec3f compute_color(const Ray &ray,const Scene &scene, Hit hit){
    Vec3f color = {0,0,0};
    Material material = scene.materials[hit.material_id - 1];
    int number_of_point_lights = scene.point_lights.size();
    color = color + material.ambient * scene.ambient_light;


    if(material.is_mirror){
        //add the color from the mirror direction.
        Ray reflection_ray = reflection_ray_calculator(ray, hit, scene);
        reflection_ray.ray_depth = ray.ray_depth + 1;
        Vec3f reflection_color = pixel_color_calculator(reflection_ray, scene);
        color = color + reflection_color * material.mirror;
     }
    for(int current_light = 0; current_light < number_of_point_lights; current_light++){
        //Hit shadow_hit = find_the_closest_hit({hit.intersectionPoint + hit.surfaceNormal * scene.shadow_ray_epsilon, normalise_vector(scene.point_lights[current_light].position - hit.intersectionPoint), true, 0}, scene);
        Vec3f light_ray = scene.point_lights[current_light].position - hit.intersectionPoint;
        Vec3f shadow_ray_dir = normalise_vector(scene.point_lights[current_light].position - (hit.intersectionPoint + hit.surfaceNormal * scene.shadow_ray_epsilon));
        Ray shadowRay = {hit.intersectionPoint + hit.surfaceNormal * scene.shadow_ray_epsilon, shadow_ray_dir, true, 0};

        Hit shadow_hit = find_the_closest_hit(shadowRay, scene);
        
        if(shadow_hit.hitHappened && shadow_hit.time < find_length(light_ray)){
            continue;
        }
        else{
            color = color + difuse_color_calculator(ray, hit, scene, current_light) + specular_color_calculator(ray, hit, scene, current_light);

        }
    }


    return color;
}


Vec3f pixel_color_calculator(const Ray &ray,const Scene &scene){
    //TODO: Calculate pixel color with shading and lights.
    Vec3f pixel_color = {0.0f,0.0f,0.0f};

    if(scene.max_recursion_depth < ray.ray_depth){
        return pixel_color;
    }

    Hit closest_hit = find_the_closest_hit(ray, scene);

    if(closest_hit.hitHappened){
        return compute_color(ray, scene, closest_hit);
    }
    else if(ray.ray_depth == 0){
        pixel_color.x = min(scene.background_color.x, 255);
        pixel_color.y = min(scene.background_color.y, 255);
        pixel_color.z = min(scene.background_color.z, 255);
    }

    return pixel_color;
}

void *render_image(void *args){
    int start_row = ((ThreadStruct *)args)->start_row;
    int end_row = ((ThreadStruct *)args)->end_row;
    int width = ((ThreadStruct *)args)->width;
    int height = ((ThreadStruct *)args)->height;
    unsigned char *image = ((ThreadStruct *)args)->image;
    Scene scene = *((ThreadStruct *)args)->scene;
    int camera_index = ((ThreadStruct *)args)->camera_index;

    for(int y = start_row; y < end_row; y++){
        for(int x = 0; x < width; x++){
            Ray ray = ray_generator(x, y, scene.cameras[camera_index]);
            Vec3f pixel_color = pixel_color_calculator(ray, scene);
            int pixel = (y * width + x) * 3;
            image[pixel] = clamp_function(pixel_color.x);
            image[pixel+1] = clamp_function(pixel_color.y);
            image[pixel+2] = clamp_function(pixel_color.z);
        }
    }

    pthread_exit(NULL);

}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    Scene scene;

    scene.loadFromXml(argv[1]);


    int number_of_cameras = scene.cameras.size();

    //This loop is for if there are multiple cameras. Doing the same calculation and
    //produces images as many as cameras.

    for(int i = 0; i < number_of_cameras; i++){
        int image_width = scene.cameras[i].image_width;
        int image_height = scene.cameras[i].image_height;

        unsigned char *image = new unsigned char[image_width * image_height * 3];

        int number_of_threads = 4;

        int row_per_thread = image_height / number_of_threads;

        pthread_t *threads = new pthread_t[number_of_threads];

        for(int j = 0; j < number_of_threads; j++){
            ThreadStruct *arg = new ThreadStruct;
            arg->start_row = j * row_per_thread;
            arg->end_row = j == number_of_threads-1 ? image_height : (j + 1) * row_per_thread;
            arg->width = image_width;
            arg->height = image_height;
            arg->image = image;
            arg->scene = &scene;
            arg->camera_index = i;
            pthread_create(&threads[j], NULL, render_image, (void *)arg);
        }

        pthread_join(threads[0], NULL);
        pthread_join(threads[1], NULL);
        pthread_join(threads[2], NULL);
        pthread_join(threads[3], NULL);

        /*
        int pixel = 0;

        for(int y = 0; y < image_height; y++){
            for(int x = 0; x < image_width; x++){
                Ray ray = ray_generator(x, y, scene.cameras[i]);
                Vec3f pixel_color = pixel_color_calculator(ray, scene);
                image[pixel++] = clamp_function(pixel_color.x);
                image[pixel++] = clamp_function(pixel_color.y);
                image[pixel++] = clamp_function(pixel_color.z);
            }
        }

        */
        //Implement threads.

        write_ppm(scene.cameras[i].image_name.c_str(), image, image_width, image_height);

        delete[] image;
    }

}


