#version 410
// Fragment shader

uniform mat3 normalmatrix;

layout (location = 0) in vec3 vertcoords_camera_fs;
layout (location = 1) in vec3 vertnormal_world_fs;

uniform uint shadingType;

out vec4 fColor;

// Computes fixed ambient + reflection
void main() {

//  if (!gl_FrontFacing) {
//    discard;
//  }

  vec3 vertnormal_camera_fs = normalmatrix * vertnormal_world_fs;

  vec3 lightpos = vec3(3.0, 0.0, 2.0);
  vec3 lightcolour = vec3(1.0);

  vec3 matcolour = vec3(0.53, 0.80, 0.87);
  vec3 matspeccolour = vec3(1.0);

  float matambientcoeff = 0.2;
  float matdiffusecoeff = 0.6;
  float matspecularcoeff = 0.4;

  vec3 normal = normalize(vertnormal_camera_fs);

  vec3 surftolight = normalize(lightpos - vertcoords_camera_fs);
  float diffusecoeff = max(0.0, dot(surftolight, normal));
  vec3 camerapos = vec3(0.0);
  vec3 surftocamera = normalize(camerapos - vertcoords_camera_fs);
  vec3 reflected = 2 * diffusecoeff * normal - surftolight;
  float specularcoeff = max(0.0, dot(reflected, surftocamera));

  vec3 compcolour = min(1.0, matambientcoeff + matdiffusecoeff * diffusecoeff) * lightcolour * matcolour;
  compcolour += matspecularcoeff * specularcoeff * lightcolour * matspeccolour;

  // ---

  if (shadingType % 2 == 0) {
    // Wireframe
    gl_FragDepth = gl_FragCoord.z - 0.0005;
  }
  else {
    gl_FragDepth = gl_FragCoord.z;
  }

  switch (shadingType) {
  case 0:
  case 1:
    fColor = vec4(compcolour, 1.0);
    break;
  case 2:
    fColor = vec4(vec3(0.0), 1.0);
    break;
  case 3:
    fColor = vec4(1.0);
    break;
  }
  //fColor = vec4(normalize(vertnormal_world_fs), 1.0);
}

