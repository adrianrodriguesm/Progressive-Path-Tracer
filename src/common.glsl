/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979;
const float epsilon = 0.001;
const float displacementBias = 1.E-4;
struct Ray
{
    vec3 origin;        // origin
    vec3 direction;     // direction - always set with normalized vector
    float time;         // time, for motion blur
};

Ray createRay(vec3 o, vec3 d, float t)
{
    Ray r;
    r.origin = o;
    r.direction = d;
    r.time = t;
    return r;
}

Ray createRay(vec3 o, vec3 d)
{
    return createRay(o, d, 0.0);
}

vec3 pointOnRay(Ray r, float t)
{
    return r.origin + r.direction * t;
}

float gSeed = 0.0;
vec3 viewDirection;
float currRefractIndex;
uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

float rand(vec2 v)
{
    return fract(sin(dot(v.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 toLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

vec3 toGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

vec2 randomInUnitDisk(inout float seed) 
{
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 randomInUnitSphere(inout float seed)
{
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

struct Camera
{
    vec3 eye;
    vec3 u, v, n; // x y z
    float width, height;
    float lensRadius;
    float planeDist, focusDist;
    float time0, time1;
};

Camera createCamera(
    vec3 eye,
    vec3 at,
    vec3 worldUp,
    float fovy,
    float aspect,
    float aperture,  //diametro em multiplos do pixel size
    float focusDist,  //focal ratio
    float time0,
    float time1)
{
    Camera cam;
    if(aperture == 0.0) cam.focusDist = 1.0; //pinhole camera then focus in on vis plane
    else cam.focusDist = focusDist;
    vec3 w = eye - at;
    cam.planeDist = length(w);
    cam.height = 2.0 * cam.planeDist * tan(fovy * pi / 180.0 * 0.5);
    cam.width = aspect * cam.height;

    cam.lensRadius = aperture * 0.5 * cam.width / iResolution.x;  //aperture ratio * pixel size; (1 pixel=lente raio 0.5)
    cam.eye = eye;
    cam.n = normalize(w);
    cam.u = normalize(cross(worldUp, cam.n));
    cam.v = cross(cam.n, cam.u);
    cam.time0 = time0;
    cam.time1 = time1;
    return cam;
}

Ray getRay(Camera cam, vec2 pixelSample)  //rnd pixel_sample viewport coordinates
{
    vec2 ls = cam.lensRadius * randomInUnitDisk(gSeed);  //ls - lens sample for DOF
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);
    
    //Calculate eye_offset and ray direction
    vec3 rayOrigin = cam.eye + ls.x * cam.u + ls.y * cam.v;
    vec3 pointInFocalPlane;
	pointInFocalPlane.x = cam.width * (pixelSample.x / iResolution.x - 0.5f) * cam.focusDist;
    pointInFocalPlane.y = cam.height * (pixelSample.y / iResolution.y - 0.5f) * cam.focusDist;
    pointInFocalPlane.z = -cam.planeDist;

    vec3 rayDirection = (pointInFocalPlane.x - ls.x) * cam.u + 
                         (pointInFocalPlane.y - ls.y) * cam.v +    
                         pointInFocalPlane.z * cam.n;

    return createRay(rayOrigin, normalize(rayDirection), time);
}

// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIALECTRIC 2

struct Material
{
    int type;
    vec3 albedo;                // the color used for diffuse lighting
    float shininess;
    float specularChance;       // percentage chance of doing a specular reflection
    float specularRoughness;    // how rough the specular reflections are
    vec3  specularColor;        // the color tint of specular reflections
    float refractionIndex;      // index of refraction. used by fresnel and refraction.
    float refractionChance;     // percent chance of doing a refractive transmission
    float refractionRoughness;  // how rough the refractive transmissions are
    vec3  refractionColor;      // absorption for beer's law 
};

Material createDiffuseMaterial(vec3 albedo)
{
    Material material;
    material.type = MT_DIFFUSE;
    material.albedo = albedo;
    material.shininess = 10.f;        
    material.specularChance = 0.f;
    material.specularRoughness = 0.0f;
    material.specularColor = vec3(.1f); // Grey 
    material.refractionIndex = 0.f;
    material.refractionChance = 0.0f;
    return material;
}
Material createZeroedMaterial()
{
    Material ret;
    ret.albedo = vec3(0.0f, 0.0f, 0.0f);
    ret.specularChance = 0.0f;
    ret.specularRoughness = 0.0f;
    ret.shininess = 10.f;
    ret.specularColor = vec3(0.0f, 0.0f, 0.0f);
    ret.refractionIndex = 1.0f;
    ret.refractionChance = 0.0f;
    ret.refractionRoughness = 0.0f;
    return ret;
}
Material createMetalMaterial(vec3 specular, float roughness)
{
    Material material;
    material.type = MT_METAL;
    material.albedo = vec3(0);
    material.shininess = 220.f;      
    material.specularChance = 1.f;
    material.specularRoughness = roughness;
    material.specularColor = specular;
    material.refractionIndex = 0.f;
    material.refractionChance = 0.0f;
    return material;
}

Material createDialectricMaterial(vec3 albedo, float refIdx)
{
    Material material;
    material.type = MT_DIALECTRIC;
    material.albedo = albedo;
    material.shininess = 550.f;          
    material.specularChance = 0.02f;
    material.specularRoughness = 0.f;
    material.specularColor = vec3(1.0f) * 0.8f; // Grey 
    material.refractionIndex = refIdx;
    material.refractionChance = 1.0f;
    return material;
}

struct HitRecord
{
    vec3 pos;
    vec3 normal;
    float t;            // ray parameter
    bool hitFromInside;
    Material material;
};

// Calculate Reflection power (Shlicks Approximation)
float schlick(float cosine, float refractionIndex, float otherRefIndex)
{
    float R0 = refractionIndex * refractionIndex;
    return R0 + (1.0f - R0) * pow(max(1.0f - cosine, 0.0), 5.0f);
}

float fresnelReflectAmount(float n1, float n2, vec3 normal, vec3 incident, float f0, float f90)
{
        // Schlick aproximation
        float r0 = (n1-n2) / (n1+n2);
        r0 *= r0;
        float cosX = -dot(normal, incident);
        if (n1 > n2)
        {
            float n = n1/n2;
            float sinT2 = n*n*(1.0-cosX*cosX);
            // Total internal reflection
            if (sinT2 > 1.0)
                return f90;
            cosX = sqrt(1.0-sinT2);
        }
        float x = 1.0-cosX;
        float ret = r0+(1.0-r0)*x*x*x*x*x;

        // adjust reflect multiplier for object reflectivity
        return mix(f0, f90, ret);
}

bool scatter(Ray rayIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    if(rec.material.type == MT_DIFFUSE)
    {
        vec3 rayOrigin = rec.pos + rec.normal * displacementBias;
        vec3 s = rayOrigin + rec.normal + randomInUnitSphere(gSeed);
        vec3 rayDirection = normalize(s - rayOrigin);
        rScattered = createRay(rec.pos, rayDirection);
        atten = rec.material.albedo * max(dot(rScattered.direction, rec.normal), 0.0) / pi;
        return true;
    }
    if(rec.material.type == MT_METAL)
    {
        vec3 rayOrigin = rec.pos + rec.normal * displacementBias;
        vec3 rayDirection = normalize(reflect(rayIn.direction, rec.normal));
        rayDirection = normalize(rayDirection + randomInUnitSphere(gSeed) * rec.material.specularRoughness);
        rScattered = createRay(rayOrigin, rayDirection);
        atten = rec.material.specularColor;
        return true;
    }
    if(rec.material.type == MT_DIALECTRIC)
    {
        atten = rec.material.albedo;
        vec3 outwardNormal;
        float niOverNt;
        float cosine;
/**/
        // Hit inside
        if(dot(rayIn.direction, rec.normal) > 0.0) 
        {
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refractionIndex;
            cosine = rec.material.refractionIndex * dot(rayIn.direction, rec.normal); 
        }
        // Hit from outside
        else 
        {
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refractionIndex;
            cosine = -dot(rayIn.direction, rec.normal);
        }

        //Use probabilistic math to decide if scatter a reflected ray or a refracted ray
/**/
        float reflectProb = schlick(cosine, niOverNt, currRefractIndex);
        if( hash1(gSeed) < reflectProb)  //Reflection
        {

            vec3 rayOrigin = rec.pos + outwardNormal * displacementBias;
            vec3 rayDirection = normalize(reflect(rayIn.direction, outwardNormal));
            rScattered = createRay(rayOrigin, rayDirection);
            atten *= vec3(reflectProb); // not necessary since we are only scattering reflectProb rays and not all reflected rays 
        }
        else
        {
            vec3 rayOrigin = rec.pos - outwardNormal * displacementBias;
            vec3 inverViewDir = rayIn.direction * -1.f;
            vec3 rayDirection = refract(inverViewDir, outwardNormal, niOverNt);
            rScattered = createRay(rayOrigin, rayDirection);
            
            atten *= vec3(1.0 - reflectProb);// not necessary since we are only scattering 1-reflectProb rays and not all refracted rays
        }

        return true;
    }
    return false;
}

struct pointLight 
{
    vec3 pos;
    vec3 color;
};

pointLight createPointLight(vec3 pos, vec3 color) 
{
    pointLight l;
    l.pos = pos;
    l.color = color;
    return l;
}

struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle createTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}

bool hit_triangle(Triangle triangle, Ray ray, float tmin, float tmax, out HitRecord rec)
{
    vec3 edge1 = triangle.b - triangle.a;
    vec3 edge2 = triangle.c - triangle.a;
    vec3 normal = cross(edge1, edge2);
    // If the determinant is zero this means that 
    // the ray is paralle to triangle plane and if
    // it is bellow zero that means is behind
    float det = -dot(ray.direction, normal);
	float invDet = 1.f / det; // We divide to avoid used floating point division
	vec3 AO = ray.origin - triangle.a;
    // Using the scalar triple product and the Cramer's rule
    // to calculate the bary-centric coords
    vec3 DAO = cross(AO, ray.direction);
    float u = dot(edge2, DAO) * invDet;
    float v = -dot(edge1, DAO) * invDet;
    float t = dot(AO, normal) * invDet;
    //calculate a valid t and normal
    if(det >= 1e-6 && t >= 0.0 && t < tmax && t > tmin && u >= 0.0 && v >= 0.0 && (u + v) <= 1.0f)
    {
        rec.t = t;
        rec.hitFromInside = false;
        rec.normal = normal;
        rec.pos = pointOnRay(ray, rec.t);
        return true;
    }
    return false;
}

struct Sphere
{
    vec3 center;
    float radius;
};

Sphere createSphere(vec3 center, float radius)
{
    Sphere s;
    s.center = center;
    s.radius = radius;
    return s;
}


struct MovingSphere
{
    vec3 center0, center1;
    float radius;
    float time0, time1;
};

MovingSphere createMovingSphere(vec3 center0, vec3 center1, float radius, float time0, float time1)
{
    MovingSphere s;
    s.center0 = center0;
    s.center1 = center1;
    s.radius = radius;
    s.time0 = time0;
    s.time1 = time1;
    return s;
}

vec3 center(MovingSphere mvsphere, float time)
{
    return vec3(0.f);
    //return moving_center;
}

/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */

bool hit_sphere(Sphere sphere, Ray ray, float tmin, float tmax, out HitRecord rec)
{
    
    // Calculate a valid t and normal
	// Center and origin inversed and signs of b inversed for sphere optimization
    vec3 temp = sphere.center - ray.origin;

    float b = dot(temp, ray.direction);
    float c = dot(temp, temp) - sphere.radius * sphere.radius;
    // If origin outside and pointing away from sphere
    if (c > 0.f && b <= 0.f)
        return false;

    float a = 1.f;    // length of ray.Direction should be 1 
    float disc = b * b - a * c;

    if (disc <= 0.0)
        return false;

    float e = sqrt(disc);
    float t = (b - e); // root 1

    if (t > epsilon && t < tmax && t > tmin)
    {
        rec.t = t;
        rec.pos = pointOnRay(ray, rec.t);
        rec.hitFromInside = (dot((rec.pos - sphere.center), ray.direction)) > 0.f;
        rec.normal = normalize(rec.pos - sphere.center);
        return true;
    }
    
    t = (b + e); //root 2

    if (t > epsilon && t < tmax && t > tmin)
    {
        rec.t = t;
        rec.pos = pointOnRay(ray, rec.t);
        rec.hitFromInside = (dot((rec.pos - sphere.center), ray.direction)) > 0.f;
        rec.normal =  normalize(rec.pos - sphere.center);
        return true;
    }
    return  false;
}

bool hit_movingSphere(MovingSphere sphere, Ray ray, float tmin, float tmax, out HitRecord rec)
{
    float B, C, delta;
    bool outside;
    float t;


     //INSERT YOUR CODE HERE
     //Calculate the moving center
    //calculate a valid t and normal
	
    if(t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(ray, rec.t);
        rec.normal = normalize(rec.pos - sphere.center0);
        return true;
    }
    else return false;
    
}
