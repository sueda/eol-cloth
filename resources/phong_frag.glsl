#version 120

varying vec3 vPos; // in camera space
varying vec3 vNor; // in camera space
varying vec2 vTex;
uniform vec3 kdFront;
uniform vec3 kdBack;

void main()
{
	vec3 lightPos = vec3(0.0, 0.0, 0.0);
	vec3 n = normalize(vNor);
	vec3 l = normalize(lightPos - vPos);
	vec3 v = -normalize(vPos);
	vec3 h = normalize(l + v);
	vec3 kd = kdFront;
	float ln = dot(l, n);
	if(ln < 0.0) {
		kd = kdBack;
		ln = -ln;
	}
	vec3 diffuse = ln * kd;
	vec3 color = diffuse;
	gl_FragColor = vec4(color, 1.0);
}
