#version 120

attribute vec4 aPos;
attribute vec3 aNor;
attribute vec2 aTex;
uniform mat4 P;
uniform mat4 MV;
varying vec3 vPos;
varying vec3 vNor;
varying vec2 vTex;

void main()
{
	vec4 posCam = MV * aPos;
	gl_Position = P * posCam;
	vPos = posCam.xyz;
	vNor = (MV * vec4(aNor, 0.0)).xyz;
	vTex = aTex;
}
