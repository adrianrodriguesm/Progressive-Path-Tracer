#iChannel0 "./PathTracer.glsl"
#include "./Common.glsl"

void main()
{
    vec3 color = texture(iChannel0, gl_FragCoord.xy / iResolution.xy).rgb;

    // convert unbounded HDR color range to SDR color range
    float exposure = 1.f;
    color = ACESFilm(color, exposure);

    // convert from linear to sRGB for display
    //color = LinearToSRGB(color);
    
    gl_FragColor = vec4(toGamma(color), 1.0f);
}