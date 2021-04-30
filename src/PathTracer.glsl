/**
* ver hash functions em
* https://www.shadertoy.com/view/XlGcRh hash functions GPU
* http://www.jcgt.org/published/0009/03/02/
 */

 #include "./Common.glsl"
 #iChannel0 "self"

bool hit_world(Ray r, float tmin, float tmax, out HitRecord rec)
{
    bool hit = false;
    rec.t = tmax;
    rec.hitFromInside = false;
    if(hit_triangle(createTriangle(vec3(-10.0, -0.01, 10.0), vec3(10.0, -0.01, 10.0), vec3(-10.0, -0.01, -10.0)), r, tmin, rec.t, rec))
    {
        hit = true;
        rec.material = createDiffuseMaterial(vec3(0.2));
    }

    if(hit_triangle(createTriangle(vec3(-10.0, -0.01, -10.0), vec3(10.0, -0.01, 10), vec3(10.0, -0.01, -10.0)), r, tmin, rec.t, rec))
    {
        hit = true;
        rec.material = createDiffuseMaterial(vec3(0.2));
    }

    if(hit_sphere(
        createSphere(vec3(-4.0, 1.0, 0.0), 1.0),
        r,
        tmin,
        rec.t,
        rec))
    {
        hit = true;
        rec.material = createDiffuseMaterial(vec3(0.2, 0.95, 0.1));
        //rec.material = createDiffuseMaterial(vec3(0.4, 0.2, 0.1));
    }

    if(hit_sphere(
        createSphere(vec3(4.0, 1.0, 0.0), 1.0),
        r,
        tmin,
        rec.t,
        rec))
    {
        hit = true;
        rec.material = createMetalMaterial(vec3(0.7, 0.6, 0.5), 0.f);
    }
/**/
    if(hit_sphere(
        createSphere(vec3(0.0, 1.0, 0.0), 1.0),
        r,
        tmin,
        rec.t,
        rec))
    {
        hit = true;
        rec.material = createDialectricMaterial(vec3(1.0), 1.5);
    }
/**/
    if(hit_sphere(
        createSphere(vec3(0.0, 1.0, 0.0), -0.55),
        r,
        tmin,
        rec.t,
        rec))
    {
        hit = true;
        rec.material = createDialectricMaterial(vec3(1.0), 1.5);
    }
/**/
    int numxy = 5;
    
    for(int x = -numxy; x < numxy; ++x)
    {
        for(int y = -numxy; y < numxy; ++y)
        {
            float fx = float(x);
            float fy = float(y);
            float seed = fx + fy / 1000.0;
            vec3 rand1 = hash3(seed);
            vec3 center = vec3(fx + 0.9 * rand1.x, 0.2, fy + 0.9 * rand1.y);
            float chooseMaterial = rand1.z;
            if(distance(center, vec3(4.0, 0.2, 0.0)) > 0.9)
            {
                if(chooseMaterial < 0.3)
                {
                    vec3 center1 = center + vec3(0.0, hash1(gSeed) * 0.5, 0.0);
                    /** /
                    // diffuse
                    if(hit_movingSphere(
                        createMovingSphere(center, center1, 0.2, 0.0, 1.0),
                        r,
                        tmin,
                        rec.t,
                        rec))
                    {
                        hit = true;
                        rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                    }
                    /**/
                }
                else if(chooseMaterial < 0.5)
                {
                    // diffuse
                    if(hit_sphere(
                        createSphere(center, 0.2),
                        r,
                        tmin,
                        rec.t,
                        rec))
                    {
                        hit = true;
                        rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                    }
                }
                else if(chooseMaterial < 0.7)
                {
                    // metal
                    if(hit_sphere(
                        createSphere(center, 0.2),
                        r,
                        tmin,
                        rec.t,
                        rec))
                    {
                        hit = true;
                       // rec.material.type = MT_METAL;
                        rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, 0.0);
                    }
                }
                else if(chooseMaterial < 0.9)
                {
                    // metal
                    if(hit_sphere(
                        createSphere(center, 0.2),
                        r,
                        tmin,
                        rec.t,
                        rec))
                    {
                        hit = true;
                       // rec.material.type = MT_METAL;
                        rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, hash1(seed));
                    }
                }
                else
                {
                    /**/
                    // glass (dialectric)
                    if(hit_sphere(
                        createSphere(center, 0.2),
                        r,
                        tmin,
                        rec.t,
                        rec))
                    {
                        hit = true;
                        rec.material.type = MT_DIALECTRIC;
                        rec.material = createDialectricMaterial(hash3(seed), 1.5);
                    }
                    /**/
                }
            }
        }
    }
    /**/
    return hit;
}

vec3 directlighting(pointLight pl, Ray r, HitRecord rec)
{
    vec3 normal = normalize(rec.normal);
    vec3 viewDirection = normalize(r.direction);
    vec3 emissionPoint = rec.pos + normal * displacementBias;
    vec3 lightDirection = pl.pos - rec.pos;
    float tMax = length(lightDirection);
    lightDirection = normalize(lightDirection);
    // Light is bellow the surface
    float lightIntensity = dot(lightDirection, normal);
    if(lightIntensity < 0.f)
        return vec3(0.);

    HitRecord lightRecord;
    Ray lightRay = createRay(emissionPoint,lightDirection);
    if(hit_world(lightRay, 0.001, tMax, lightRecord))
        return vec3(0.f);
    
    // Diffuse
    vec3 diffuseColor = (pl.color *  rec.material.diffusePercent * rec.material.albedo) * lightIntensity;
    // Specular
    vec3 halfwayVector = normalize(-viewDirection + lightDirection);
    float specAngle = max(dot(halfwayVector, normal), 0.f);
    float ksSpecular = pow(specAngle, rec.material.shininess) * rec.material.specularPercent;
    vec3 specularColor = pl.color *  ksSpecular * rec.material.specularColor;
	return diffuseColor + specularColor;
}

#define MAX_BOUNCES 10

vec3 rayColor(Ray ray)
{
    HitRecord rec;

    vec3 color = vec3(0.0);
    vec3 throughput = vec3(1.0f, 1.0f, 1.0f);
    for(int i = 0; i < MAX_BOUNCES; ++i)
    {   
        if(hit_world(ray, 0.001, 10000.0, rec))
        {      
            // Calculate direct lighting with 3 white point lights:
            {
                pointLight pl0 = createPointLight(vec3(-10.0, 15.0, 0.0), vec3(1.0, 1.0, 1.0));
                pointLight pl1 = createPointLight(vec3(8.0, 15.0, 3.0), vec3(1.0, 1.0, 1.0));
                pointLight pl2 = createPointLight(vec3(1.0, 15.0, -9.0), vec3(1.0, 1.0, 1.0));

                color += directlighting(pl0, ray, rec) * throughput;
                color += directlighting(pl1, ray, rec) * throughput;
                color += directlighting(pl2, ray, rec) * throughput;    
            }
           
            // Calculate secondary ray and update throughput
            Ray scatterRay;
            vec3 atten;
            if(scatter(ray, rec, atten, scatterRay))
            {    
                throughput *= atten;
                ray = scatterRay;

                // Russian Roulette
                // As the throughput gets smaller, the ray is more likely to get terminated early.
                // Survivors have their value boosted to make up for fewer samples being in the average.
                {
                    float p = max(throughput.r, max(throughput.g, throughput.b));
                    if (hash1(gSeed) > p)
                        break;
                
                    // Add the energy we 'lose' by randomly terminating paths
                    throughput *= 1.0f / p;
                }
            }
        }
        else  // Background
        {
            float t = 0.8 * (ray.direction.y + 1.0);
            color += throughput * mix(vec3(1.0), vec3(0.5, 0.7, 1.0), t);
            break;
        }
    }
    return color;
}

#define MAX_SAMPLES 5000.0

void main()
{
    gSeed = float(baseHash(floatBitsToUint(gl_FragCoord.xy))) / float(0xffffffffU) + iTime;

    vec2 mouse = iMouse.xy / iResolution.xy;
    mouse.x = mouse.x * 2.0 - 1.0;

    vec3 camPos = vec3(mouse.x * 10.0, mouse.y * 5.0, 8.0);
    vec3 camTarget = vec3(0.0, 0.0, -1.0);
    float fovy = 60.0;
    float aperture = 0.0;
    float distToFocus = 2.5;
    float time0 = 0.0;
    float time1 = 1.0;
    Camera cam = createCamera(
        camPos,
        camTarget,
        vec3(0.0, 1.0, 0.0),    // world up vector
        fovy,
        iResolution.x / iResolution.y,
        aperture,
        distToFocus,
        time0,
        time1);

    

    vec4 prev = texture(iChannel0, gl_FragCoord.xy / iResolution.xy);
    vec3 prevLinear = toLinear(prev.xyz);  

    vec2 ps = gl_FragCoord.xy + hash2(gSeed);
    //vec2 ps = gl_FragCoord.xy;
    vec3 color = rayColor(getRay(cam, ps));
    if(iMouseButton.x != 0.0 || iMouseButton.y != 0.0)
    {
        gl_FragColor = vec4(toGamma(color), 1.0);  //samples number reset = 1
        return;
    }
    if(prev.w > MAX_SAMPLES)   
    {
        gl_FragColor = prev;
        return;
    }

    float w = prev.w + 1.f;
    color = mix(prevLinear, color, 1.0/w); 
    gl_FragColor = vec4(toGamma(color), w);
}
