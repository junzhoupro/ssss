#version 410
// Fragment shader

layout (location = 1) in vec3 vertnormal_camera_fs;

out vec4 fColor;

void main() {

	fColor = vec4(vertnormal_camera_fs, 1.0);

}
