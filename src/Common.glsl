/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979f;
const float epsilon = 0.001f;
const float displacementBias = 1.E-4f;
struct Ray
{
    vec3 origin;        // origin
    vec3 direction;     // direction - always set with normalized vector
    float time;         // time, for motion blur
};

Ray CreateRay(vec3 o, vec3 d, float t)
{
    Ray r;
    r.origin = o;
    r.direction = d;
    r.time = t;
    return r;
}

Ray CreateRay(vec3 o, vec3 d)
{
    return CreateRay(o, d, 0.0);
}

vec3 PointOnRay(Ray r, float t)
{
    return r.origin + r.direction * t;
}

float gSeed = 0.0;
uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) 
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) 
{
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

vec3 ToLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

vec3 ToGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

vec3 ACESFilm(vec3 color, float exposure)
{
    // Apply exposure (how long the shutter is open)
    color = color * exposure;
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    return clamp((color*(a*color + b)) / (color*(c*color + d) + e), 0.0f, 1.0f);
}

vec2 RandomInUnitDisk(inout float seed) 
{
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 RandomInUnitSphere(inout float seed)
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

Camera CreateCamera(
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

    // Aperture ratio * pixel size; (1 pixel=lente raio 0.5)
    cam.lensRadius = aperture * 0.5f * cam.width / iResolution.x;  
    cam.eye = eye;
    cam.n = normalize(w);
    cam.u = normalize(cross(worldUp, cam.n));
    cam.v = cross(cam.n, cam.u);
    cam.time0 = time0;
    cam.time1 = time1;
    return cam;
}

Ray GetPrimaryRay(Camera cam, vec2 pixelSample)  //rnd pixel_sample viewport coordinates
{
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);
    
    // ls - lens sample for DOF
    vec2 ls = cam.lensRadius * RandomInUnitDisk(gSeed);  
    //Calculate eye_offset and ray direction
    vec3 rayOrigin = cam.eye + ls.x * cam.u + ls.y * cam.v;
    float focalRatio = (cam.lensRadius > 0.f) ? cam.focusDist/cam.planeDist : 1.f;
    vec3 pointInFocalPlane;
	pointInFocalPlane.x = cam.width * (pixelSample.x / iResolution.x - 0.5f) * focalRatio;
    pointInFocalPlane.y = cam.height * (pixelSample.y / iResolution.y - 0.5f) * focalRatio;
    pointInFocalPlane.z = -cam.planeDist;

    vec3 rayDirection = (pointInFocalPlane.x - ls.x) * cam.u + 
                         (pointInFocalPlane.y - ls.y) * cam.v +    
                         pointInFocalPlane.z * cam.n;

    return CreateRay(rayOrigin, normalize(rayDirection), time);
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
    vec3 emissive;              // how much the surface glows
    float specularPercent;      // percentage chance of doing a specular reflection
    float diffusePercent;       // percentage chance of doing a specular reflection
    float specularRoughness;    // how rough the specular reflections are
    vec3  specularColor;        // the color tint of specular reflections
    float refractionRoughness;   // how rough the refractive transmissions are
    float refractionChance;
    float specularChance;
    float refractionIndex;      // index of refraction. used by fresnel and refraction.
    vec3 refractionColor;       // absorption for beer's law (color is how much light to block)
};
Material createZeroedMaterial()
{
    Material ret;
    ret.albedo = vec3(0.0f, 0.0f, 0.0f);
    ret.specularPercent = 0.0f;
    ret.specularRoughness = 0.0f;
    ret.shininess = 10.f;
    ret.diffusePercent = 0.f;
    ret.specularColor = vec3(0.0f, 0.0f, 0.0f);
    ret.refractionIndex = 1.0f;
    ret.refractionColor = vec3(0.0f, 0.0f, 0.0f);
    return ret;
}
// Some of this values could be sent as a parameter (e.g shininess)
Material CreateDiffuseMaterial(vec3 albedo)
{
    Material material;
    material.type = MT_DIFFUSE;
    material.albedo = albedo;
    material.emissive = vec3(0.f);
    material.shininess = 10.f;
    material.diffusePercent = 1.f;
    material.specularPercent = 0.f;
    material.specularRoughness = 0.0f;
    material.specularColor = vec3(.1f); // Grey 
    material.refractionIndex = 0.f;
    material.refractionColor = vec3(0.0f, 0.0f, 0.0f);
    return material;
}

Material CreateMetalMaterial(vec3 specular, float roughness)
{
    Material material;
    material.type = MT_METAL;
    material.albedo = vec3(0);
    material.emissive = vec3(0.f);
    material.shininess = 220.f;
    material.diffusePercent = 0.f;     
    material.specularPercent = 1.f;
    material.specularRoughness = roughness;
    material.specularColor = specular;
    material.refractionIndex = 0.f;
    material.refractionColor = vec3(0.0f, 0.0f, 0.0f);
    return material;
}

Material CreateDialectricMaterial(vec3 albedo, float refIdx)
{
    Material material;
    material.type = MT_DIALECTRIC;
    material.albedo = albedo;
    material.emissive = vec3(0.f);
    material.shininess = 20.f;
    material.diffusePercent = 0.f;       
    material.specularPercent = 0.7f;
    material.specularRoughness = 0.5f;
    material.specularColor = vec3(0); 
    material.refractionIndex = refIdx;
    material.refractionColor = vec3(0.0f, 0.0f, 0.0f);
    return material;
}

struct HitRecord
{
    vec3 pos;
    vec3 normal;
    float t;        // Ray parameter
    bool hitFromInside;
    Material material;
};

// Calculate Reflection power (Shlicks Approximation)
float FresnelSchlick(float n1, float n2, float incidentAngle)
{
        // Schlick aproximation
        float r0 = (n1-n2) / (n1+n2);
        r0 *= r0;
        float cosX = incidentAngle;
        if(n1 >= n2)
        {
            float n = n1/n2;
            float sinT2 = n*n*(1.0-cosX*cosX);
            // Total internal reflection
            if (sinT2 > 1.0)
                return 1.f;

            cosX = sqrt(1.0-sinT2);
        }
        float x = 1.0 - cosX;
        return  r0 + (1.0-r0) *x*x*x*x*x;
}

bool Scatter(Ray rayIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    vec3 normal = rec.normal;

    if(rec.material.type == MT_DIFFUSE)
    {
        vec3 rayOrigin = rec.pos + normal * displacementBias;
        vec3 s = rayOrigin + normal + RandomInUnitSphere(gSeed);
        vec3 rayDirection = normalize(s - rayOrigin);
        rScattered = CreateRay(rec.pos, rayDirection, rayIn.time);
        atten = rec.material.albedo * max(dot(rScattered.direction, normal), 0.0) / pi;
        return true;
    }
    if(rec.material.type == MT_METAL)
    {
        vec3 rayOrigin = rec.pos + normal * displacementBias;
        vec3 rayDirection = normalize(reflect(rayIn.direction, normal));
        rayDirection = normalize(rayDirection + RandomInUnitSphere(gSeed) * rec.material.specularRoughness);
        rScattered = CreateRay(rayOrigin, rayDirection, rayIn.time);
        atten = rec.material.specularColor;
        return true;
    }
    if(rec.material.type == MT_DIALECTRIC)
    {
        vec3 rayDir = rayIn.direction;
        vec3 inveRayDir = rayDir * -1.f;

        float incidentAngle = (!rec.hitFromInside) ? -dot(normal, rayDir) : dot(normal, rayDir);
        // Take fresnel into account for specularChance and adjust other chances.
        float reflectiveWeight = FresnelSchlick(
                !rec.hitFromInside ? 1.f :  rec.material.refractionIndex,
                !rec.hitFromInside ? rec.material.refractionIndex : 1.f,
                incidentAngle);
        // Calculate whether we are going to do a reflective, or refractive ray
        float rayProbability = 1.0f;
        if (hash1(gSeed) < reflectiveWeight)
        {
            rayProbability = reflectiveWeight;
            rScattered.origin =  rec.pos + normal * displacementBias;
            vec3 specularRayDir = normalize(reflect(rayDir, normal));         
            rScattered.direction = normalize(specularRayDir + RandomInUnitSphere(gSeed) * rec.material.specularRoughness); 
        }
        else
        {
            rayProbability = 1.f - reflectiveWeight;
            rScattered.origin =  rec.pos - normal * displacementBias;
            vec3 rayDir = refract(rayDir, normal, rec.hitFromInside ? rec.material.refractionIndex : 1.0f / rec.material.refractionIndex);
            rScattered.direction = normalize(mix(normalize(rayDir), normalize(normal + RandomInUnitSphere(gSeed)), rec.material.refractionRoughness * rec.material.refractionRoughness));
        }
        rScattered.time = rayIn.time; 
        atten = rec.material.albedo * rayProbability;
        return true;
    }
    return false;
}

struct PointLight 
{
    vec3 pos;
    vec3 color;
};

PointLight CreatePointLight(vec3 pos, vec3 color) 
{
    PointLight l;
    l.pos = pos;
    l.color = color;
    return l;
}

struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle CreateTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}

bool HitTriangle(Triangle triangle, Ray ray, float tmin, float tmax, out HitRecord rec)
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
        rec.normal = normalize(normal);
        rec.pos = PointOnRay(ray, rec.t);
        return true;
    }
    return false;
}

struct Sphere
{
    vec3 center;
    float radius;
};

Sphere CreateSphere(vec3 center, float radius)
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

MovingSphere CreateMovingSphere(vec3 center0, vec3 center1, float radius)
{
    MovingSphere s;
    s.center0 = center0;
    s.center1 = center1;
    s.radius = radius;
    s.time0 = 0.f;
    s.time1 = 1.f;
    return s;
}

vec3 Center(MovingSphere mvsphere, float time)
{
    return mvsphere.center0 * (1.f - time) + mvsphere.center1 * time;   
}

/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */
bool HitSphere(Sphere sphere, Ray ray, float tmin, float tmax, out HitRecord rec)
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
        rec.pos = PointOnRay(ray, rec.t);
        vec3 hitToCenter = sphere.radius >= 0.f ? rec.pos - sphere.center : sphere.center - rec.pos;
        rec.normal = normalize(hitToCenter);
        rec.hitFromInside = (dot(rec.normal, ray.direction)) > 0.f;
        rec.normal *= (rec.hitFromInside ? -1.f : 1.f);
        return true;
    }
    
    t = (b + e); //root 2

    if (t > epsilon && t < tmax && t > tmin)
    {
        rec.t = t;
        rec.pos = PointOnRay(ray, rec.t);
        vec3 hitToCenter = sphere.radius >= 0.f ? rec.pos - sphere.center : sphere.center - rec.pos;
        rec.normal = normalize(hitToCenter);
        rec.hitFromInside = (dot(rec.normal, ray.direction)) > 0.f;
        rec.normal *= (rec.hitFromInside ? -1.f : 1.f);
        return true;
    }
    return  false;
}

bool HitMovingSphere(MovingSphere sphere, Ray ray, float tmin, float tmax, out HitRecord rec)
{
    // Calculate a valid t and normal
	// Center and origin inversed and signs of b inversed for sphere optimization
    vec3 temp = Center(sphere, ray.time) - ray.origin;

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
        rec.pos = PointOnRay(ray, rec.t);
        rec.hitFromInside = (dot((rec.pos - Center(sphere, ray.time)), ray.direction)) > 0.f;
        rec.normal = normalize(rec.pos - Center(sphere, ray.time));
        rec.normal *= (rec.hitFromInside ? -1.f : 1.f);
        return true;
    }
    
    t = (b + e); //root 2

    if (t > epsilon && t < tmax && t > tmin)
    {
        rec.t = t;
        rec.pos = PointOnRay(ray, rec.t);
        rec.hitFromInside = (dot((rec.pos - Center(sphere, ray.time)), ray.direction)) > 0.f;
        rec.normal =  normalize(rec.pos - Center(sphere, ray.time));
        rec.normal *= (rec.hitFromInside ? -1.f : 1.f);
        return true;
    }
    return  false;
    
}

struct Quad
{
    vec3 a,b,c,d;
};

Quad CreateQuad(vec3 a, vec3 b, vec3 c,vec3 d)
{
    Quad q;
    q.a = a;
    q.b = b;
    q.c = c;
    q.d = d;
    return q;
}

float ScalarTriple(vec3 u, vec3 v, vec3 w)
{
    return dot(cross(u, v), w);
}

bool HitQuad(in Quad quad, in Ray ray, float tmin, float tmax, inout HitRecord info)
{
    // calculate normal and flip vertices order if needed
    vec3 normal = normalize(cross(quad.c-quad.a, quad.c-quad.b));
    if (dot(normal, ray.direction) > 0.0f)
    {
        normal *= -1.0f;
        
		vec3 temp = quad.d;
        quad.d = quad.a;
        quad.a = temp;
        
        temp = quad.b;
        quad.b = quad.c;
        quad.c = temp;
    }
    
    vec3 p = ray.origin;
    vec3 q = ray.origin + ray.direction;
    vec3 pq = q - p;
    vec3 pa = quad.a - p;
    vec3 pb = quad.b - p;
    vec3 pc = quad.c - p;
    
    // determine which triangle to test against by testing against diagonal first
    vec3 m = cross(pc, pq);
    float v = dot(pa, m);
    vec3 intersectPos;
    if (v >= 0.0f)
    {
        // test against triangle a,b,c
        float u = -dot(pb, m);
        if (u < 0.0f) return false;
        float w = ScalarTriple(pq, pb, pa);
        if (w < 0.0f) return false;
        float denom = 1.0f / (u+v+w);
        u*=denom;
        v*=denom;
        w*=denom;
        intersectPos = u*quad.a+v*quad.b+w*quad.c;
    }
    else
    {
        vec3 pd = quad.d - p;
        float u = dot(pd, m);
        if (u < 0.0f) return false;
        float w = ScalarTriple(pq, pa, pd);
        if (w < 0.0f) return false;
        v = -v;
        float denom = 1.0f / (u+v+w);
        u*=denom;
        v*=denom;
        w*=denom;
        intersectPos = u*quad.a+v*quad.d+w*quad.c;
    }
    
    float dist;
    if (abs(ray.direction.x) > 0.1f)
    {
        dist = (intersectPos.x - ray.origin.x) / ray.direction.x;
    }
    else if (abs(ray.direction.y) > 0.1f)
    {
        dist = (intersectPos.y - ray.origin.y) / ray.direction.y;
    }
    else
    {
        dist = (intersectPos.z - ray.origin.z) / ray.direction.z;
    }
    
	if (dist > tmin && dist < tmax)
    {
        info.hitFromInside = false;
        info.t = dist;        
        info.normal = normal;        
        return true;
    }    
    
    return false;
}