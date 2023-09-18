#version 410
// Fragment shader

uniform float uniIsoFreq;
uniform mat3 normalmatrix;

layout (location = 0) in vec3 vertcoords_camera_fs;
layout (location = 1) in vec3 vertnormal_world_fs;

out vec4 fColor;

void main() {
	vec3 vertnormal_camera_fs = normalmatrix * vertnormal_world_fs;

	// fixed normal to compute similarity with
        vec3 control_normal = normalize(vec3(1.0, 1.0, 1.0));

	// compute dot-product
	float dotprod = dot(control_normal, normalize(vertnormal_camera_fs));

	// default: black color
        fColor = vec4(0.0, 0.0, 0.6, 1.0);

	// For 3/4 of the values -> white color
        if(sin(uniIsoFreq*dotprod) > 0.0) { //> 0.0 // -0.5
                fColor = vec4(0.4, 1.0, 1.0, 1.0);
	}
}

