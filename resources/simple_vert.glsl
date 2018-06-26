#version 120

uniform mat4 P;
uniform mat4 MV;
varying vec3 color;

void main()
{
	gl_Position = P * MV * gl_Vertex;
	color = gl_Color.rgb;
}
