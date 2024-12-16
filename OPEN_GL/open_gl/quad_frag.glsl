#version 330 core

uniform float scale;
uniform float offset;
uniform float offsetZ;
uniform vec3 color1; // Color for one type of square (e.g., black)
uniform vec3 color2; // Color for the other type of square (e.g., white)
out vec4 FragColor;

in vec3 FragPos;

void main() {

    bool x = bool(int((FragPos.x + offset + .6f) * scale * 1.25) % 2);
    bool y = bool(int((FragPos.y + offset) * scale * 1.25) % 2);
    bool z = bool(int((FragPos.z + offsetZ) * scale) % 2);
    bool xorXY = x != z;

    if (xorXY)
        FragColor = vec4(color1, 1.0); // First color (e.g., black)
    else
        FragColor = vec4(color2, 1.0); // Second color (e.g., white)
}
