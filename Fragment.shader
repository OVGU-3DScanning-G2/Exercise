varying vec3 N;
varying vec3 v;
varying vec4 FrontColor;

void main(void)
{
	vec3 lightPosition = vec3(1.0, 1.0, 1.0);
	vec4 ambientColor = vec4(0.1, 0.1, 0.1, 1.0);

	vec4 specularColor = vec4(0.7, 0.7, 0.7, 1.0);
	float shininess = 100.0;
	vec3 N = normalize(N);
	vec3 L = normalize(lightPosition - v);

	vec3 E = normalize(-v);
	vec3 R = normalize(-reflect(L, N));

	vec4 Iamb = ambientColor;
	vec4 Idiff = FrontColor * max(abs(dot(N, L)), 0.0);

	Idiff = clamp(Idiff, 0.0, 1.0);

	vec4 Ispec = specularColor * pow(max(dot(R, E), 0.0), 0.03 * shininess);
	Ispec = clamp(Ispec, 0.0, 1.0);
	gl_FragColor = Iamb + Idiff + Ispec;
}