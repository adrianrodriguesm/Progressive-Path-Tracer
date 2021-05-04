#iChannel0 "./PathTracer.glsl"
#include "./Common.glsl"

void main()
{
    vec3 color = texture(iChannel0, gl_FragCoord.xy / iResolution.xy).rgb;

    // Convert unbounded HDR color range to SDR color range
    float exposure = 0.9f;
    color = ACESFilm(color, exposure);
    gl_FragColor = vec4(ToGamma(color), 1.0f);
}