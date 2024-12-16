#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <GL/glew.h>
//#include <OpenGL/gl3.h>   // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp> // GL Math library header
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp> 
#include <ft2build.h>  //For text rendering
#include FT_FREETYPE_H  //For text rendering
#include <map> //For text rendering

#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;

int my_rand = 2;
float score = 1.f;

float elapsedTime = 0.0f;

float speed = 1.0f; // Adjust this for the desired speed
float bunnySpeed = 1.0f; // Adjust this for the desired speed
float groundZMovement = 0.0f; // This will be used to move the ground quad

float box_position = 0.0f;
float reset = 0.0f;
float groundPos = 0.0f;

GLuint gProgram[4];
int gWidth = 1000, gHeight = 800;

GLint modelingMatrixLoc[3];
GLint viewingMatrixLoc[3];
GLint projectionMatrixLoc[3];
GLint eyePosLoc[3];

glm::mat4 projectionMatrix;
glm::mat4 viewingMatrix;
glm::mat4 modelingMatrix;
glm::vec3 eyePos(0, 0, 0);

float bunny_x = 0;
float bunny_y = 0;

bool move_right = false;
bool move_left = false;

bool restart = false;

bool hit_yellow = false;
bool hit_red = false;
int hit_box = -1;
bool bunny_yay = false;
float rotation_y = 0.0f;

glm::vec3 score_color(1.f, 1.f, 0.f);


struct Vertex
{
	Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
	GLfloat x, y, z;
};

struct Texture
{
	Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) { }
	GLfloat u, v;
};

struct Normal
{
	Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
	GLfloat x, y, z;
};

struct Face
{
	Face(int v[], int t[], int n[]) {
		vIndex[0] = v[0];
		vIndex[1] = v[1];
		vIndex[2] = v[2];
		tIndex[0] = t[0];
		tIndex[1] = t[1];
		tIndex[2] = t[2];
		nIndex[0] = n[0];
		nIndex[1] = n[1];
		nIndex[2] = n[2];
	}
	GLuint vIndex[3], tIndex[3], nIndex[3];
};


// Objects
vector<Vertex> gVertices[3];
vector<Normal> gNormals[3];
vector<Texture> gTextures[3];
vector<Face> gFaces[3];

/// Holds all state information relevant to a character as loaded using FreeType
struct Character {
    GLuint TextureID;   // ID handle of the glyph texture
    glm::ivec2 Size;    // Size of glyph
    glm::ivec2 Bearing;  // Offset from baseline to left/top of glyph
    GLuint Advance;    // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;


GLuint gTextVBO, gTextVAO;
GLuint gVertexAttribBuffer[3], gIndexBuffer[3];
int gVertexDataSizeInBytes[3], gNormalDataSizeInBytes[3];
bool ParseObj(const string& fileName, int i)
{
	cout << "ParseObj(" << i << ")" << endl;
	fstream myfile;

	// Open the input 
	myfile.open(fileName.c_str(), std::ios::in);

	if (myfile.is_open())
	{
		string curLine;

		while (getline(myfile, curLine))
		{
			stringstream str(curLine);
			GLfloat c1, c2, c3;
			GLuint index[9];
			string tmp;

			if (curLine.length() >= 2)
			{
				if (curLine[0] == 'v')
				{
					if (curLine[1] == 't') // texture
					{
						str >> tmp; // consume "vt"
						str >> c1 >> c2;
						gTextures[i].push_back(Texture(c1, c2));
					}
					else if (curLine[1] == 'n') // normal
					{
						str >> tmp; // consume "vn"
						str >> c1 >> c2 >> c3;
						gNormals[i].push_back(Normal(c1, c2, c3));
					}
					else // vertex
					{
						str >> tmp; // consume "v"
						str >> c1 >> c2 >> c3;
						gVertices[i].push_back(Vertex(c1, c2, c3));
					}
				}
				else if (curLine[0] == 'f') // face
				{
					str >> tmp; // consume "f"
					char c;
					int vIndex[3], nIndex[3], tIndex[3];
					str >> vIndex[0]; str >> c >> c; // consume "//"
					str >> nIndex[0];
					str >> vIndex[1]; str >> c >> c; // consume "//"
					str >> nIndex[1];
					str >> vIndex[2]; str >> c >> c; // consume "//"
					str >> nIndex[2];

					assert(vIndex[0] == nIndex[0] &&
						vIndex[1] == nIndex[1] &&
						vIndex[2] == nIndex[2]); // a limitation for now

					// make indices start from 0
					for (int c = 0; c < 3; ++c)
					{
						vIndex[c] -= 1;
						nIndex[c] -= 1;
						tIndex[c] -= 1;
					}

					gFaces[i].push_back(Face(vIndex, tIndex, nIndex));
				}
				else
				{
					cout << "Ignoring unidentified line in obj file: " << curLine << endl;
				}
			}

			//data += curLine;
			if (!myfile.eof())
			{
				//data += "\n";
			}
		}

		myfile.close();
	}
	else
	{
		return false;
	}

	/*
	for (int i = 0; i < gVertices.size(); ++i)
	{
		Vector3 n;

		for (int j = 0; j < [i].size(); ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				if ([i][j].vIndex[k] == i)
				{
					// face j contains vertex i
					Vector3 a(gVertices[[i][j].vIndex[0]].x,
							  gVertices[[i][j].vIndex[0]].y,
							  gVertices[[i][j].vIndex[0]].z);

					Vector3 b(gVertices[[i][j].vIndex[1]].x,
							  gVertices[[i][j].vIndex[1]].y,
							  gVertices[[i][j].vIndex[1]].z);

					Vector3 c(gVertices[[i][j].vIndex[2]].x,
							  gVertices[[i][j].vIndex[2]].y,
							  gVertices[[i][j].vIndex[2]].z);

					Vector3 ab = b - a;
					Vector3 ac = c - a;
					Vector3 normalFromThisFace = (ab.cross(ac)).getNormalized();
					n += normalFromThisFace;
				}

			}
		}

		n.normalize();

		gNormals.push_back(Normal(n.x, n.y, n.z));
	}
	*/

	assert(gVertices[i].size() == gNormals[i].size());

	return true;
}

bool ReadDataFromFile(
	const string& fileName, ///< [in]  Name of the shader file
	string& data)     ///< [out] The contents of the file
{
	fstream myfile;

	// Open the input 
	myfile.open(fileName.c_str(), std::ios::in);

	if (myfile.is_open())
	{
		string curLine;

		while (getline(myfile, curLine))
		{
			data += curLine;
			if (!myfile.eof())
			{
				data += "\n";
			}
		}

		myfile.close();
	}
	else
	{
		return false;
	}

	return true;
}

void createVS(GLuint& program, const string& filename)
{
    string shaderSource;

    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &shader, &length);
    glCompileShader(vs);

    char output[1024] = {0};
    glGetShaderInfoLog(vs, 1024, &length, output);
    printf("VS compile log: %s\n", output);

    glAttachShader(program, vs);
}

void createFS(GLuint& program, const string& filename)
{
    string shaderSource;

    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &shader, &length);
    glCompileShader(fs);

    char output[1024] = {0};
    glGetShaderInfoLog(fs, 1024, &length, output);
    printf("FS compile log: %s\n", output);

    glAttachShader(program, fs);
}

void initShaders()
{
	// Create the programs

    gProgram[0] = glCreateProgram(); 
    gProgram[1] = glCreateProgram(); 
    gProgram[2] = glCreateProgram();
	gProgram[3] = glCreateProgram();

    // Create the shaders for each program

    createVS(gProgram[0], "bunny_vert.glsl");
	createFS(gProgram[0], "bunny_frag.glsl");

	createVS(gProgram[1], "quad_vert.glsl");
	createFS(gProgram[1], "quad_frag.glsl");

	createVS(gProgram[2], "cube_vert.glsl");
	createFS(gProgram[2], "cube_frag.glsl");

	createVS(gProgram[3], "text_vert.glsl");
	createFS(gProgram[3], "text_frag.glsl");
	

	glBindAttribLocation(gProgram[0], 0, "inVertex");
    glBindAttribLocation(gProgram[0], 1, "inNormal");
    glBindAttribLocation(gProgram[1], 0, "inVertex");
    glBindAttribLocation(gProgram[1], 1, "inNormal");
    glBindAttribLocation(gProgram[2], 0, "inVertex");
    glBindAttribLocation(gProgram[2], 1, "inNormal");
    glBindAttribLocation(gProgram[3], 2, "vertex");



    // Link the programs
	GLint status;

    glLinkProgram(gProgram[0]);
	

	glLinkProgram(gProgram[1]);
	

    glLinkProgram(gProgram[2]);
	

	glLinkProgram(gProgram[3]);
	

    // Get the locations of the uniform variables from each program

    for (int i = 0; i < 3; ++i)
    {
        modelingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "modelingMatrix");
        viewingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "viewingMatrix");
        projectionMatrixLoc[i] = glGetUniformLocation(gProgram[i], "projectionMatrix");
        eyePosLoc[i] = glGetUniformLocation(gProgram[i], "eyePos");
    }
}

void initVBO(int obj)
{
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	assert(glGetError() == GL_NONE);

	glGenBuffers(1, &gVertexAttribBuffer[obj]);
	glGenBuffers(1, &gIndexBuffer[obj]);

	assert(gVertexAttribBuffer[obj] > 0 && gIndexBuffer[obj] > 0);

	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer[obj]);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer[obj]);

	gVertexDataSizeInBytes[obj] = gVertices[obj].size() * 3 * sizeof(GLfloat);
	gNormalDataSizeInBytes[obj] = gNormals[obj].size() * 3 * sizeof(GLfloat);
	int indexDataSizeInBytes = gFaces[obj].size() * 3 * sizeof(GLuint);
	GLfloat* vertexData = new GLfloat[gVertices[obj].size() * 3];
	GLfloat* normalData = new GLfloat[gNormals[obj].size() * 3];
	GLuint* indexData = new GLuint[gFaces[obj].size() * 3];

	float minX = 1e6, maxX = -1e6;
	float minY = 1e6, maxY = -1e6;
	float minZ = 1e6, maxZ = -1e6;

	for (int i = 0; i < gVertices[obj].size(); ++i)
	{
		vertexData[3 * i] = gVertices[obj][i].x;
		vertexData[3 * i + 1] = gVertices[obj][i].y;
		vertexData[3 * i + 2] = gVertices[obj][i].z;

		minX = std::min(minX, gVertices[obj][i].x);
		maxX = std::max(maxX, gVertices[obj][i].x);
		minY = std::min(minY, gVertices[obj][i].y);
		maxY = std::max(maxY, gVertices[obj][i].y);
		minZ = std::min(minZ, gVertices[obj][i].z);
		maxZ = std::max(maxZ, gVertices[obj][i].z);
	}

	std::cout << "minX = " << minX << std::endl;
	std::cout << "maxX = " << maxX << std::endl;
	std::cout << "minY = " << minY << std::endl;
	std::cout << "maxY = " << maxY << std::endl;
	std::cout << "minZ = " << minZ << std::endl;
	std::cout << "maxZ = " << maxZ << std::endl;

	for (int i = 0; i < gNormals[obj].size(); ++i)
	{
		normalData[3 * i] = gNormals[obj][i].x;
		normalData[3 * i + 1] = gNormals[obj][i].y;
		normalData[3 * i + 2] = gNormals[obj][i].z;
	}

	for (int i = 0; i < gFaces[obj].size(); ++i)
	{
		indexData[3 * i] = gFaces[obj][i].vIndex[0];
		indexData[3 * i + 1] = gFaces[obj][i].vIndex[1];
		indexData[3 * i + 2] = gFaces[obj][i].vIndex[2];
	}


	glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes[obj] + gNormalDataSizeInBytes[obj], 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes[obj], vertexData);
	glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes[obj], gNormalDataSizeInBytes[obj], normalData);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

	// done copying to GPU memory; can free now from CPU memory
	delete[] vertexData;
	delete[] normalData;
	delete[] indexData;

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes[obj]));
}

void initFonts(int windowWidth, int windowHeight)
{
    // Set OpenGL options
    //glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(windowWidth), 0.0f, static_cast<GLfloat>(windowHeight));
    glUseProgram(gProgram[3]);
    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "projection"), 1, GL_FALSE, glm::value_ptr(projection));

    // FreeType
    FT_Library ft;
    // All functions return a value different than 0 whenever an error occurred
    if (FT_Init_FreeType(&ft))
    {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
    }

    // Load font as face
    FT_Face face;
    if (FT_New_Face(ft, "/usr/share/fonts/truetype/liberation/LiberationSerif-Italic.ttf", 0, &face))
    {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
    }

    // Set size to load glyphs as
    FT_Set_Pixel_Sizes(face, 0, 48);

    // Disable byte-alignment restriction
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 

    // Load first 128 characters of ASCII set
    for (GLubyte c = 0; c < 128; c++)
    {
        // Load character glyph 
        if (FT_Load_Char(face, c, FT_LOAD_RENDER))
        {
            std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
            continue;
        }
        // Generate texture
        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(
                GL_TEXTURE_2D,
                0,
                GL_RED,
                face->glyph->bitmap.width,
                face->glyph->bitmap.rows,
                0,
                GL_RED,
                GL_UNSIGNED_BYTE,
                face->glyph->bitmap.buffer
                );
        // Set texture options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // Now store character for later use
        Character character = {
            texture,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            face->glyph->advance.x
        };
        Characters.insert(std::pair<GLchar, Character>(c, character));
    }
	//glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glBindTexture(GL_TEXTURE_2D, 0);
    // Destroy FreeType once we're finished
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    //
    // Configure VBO for texture quads
    //
	//glGenVertexArrays(1, &gTextVAO);
	//glBindVertexArray(gTextVAO);
    glGenBuffers(1, &gTextVBO);
    glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void init()
{
	ParseObj("bunny.obj", 0);
	ParseObj("quad.obj", 1);
	ParseObj("cube.obj", 2);

	glEnable(GL_DEPTH_TEST);
	initShaders();
	initFonts(gWidth, gHeight);
	initVBO(0);
	initVBO(1);
	initVBO(2);
}

void drawModel(int i)
{
    glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer[i]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer[i]);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes[i]));

    glDrawElements(GL_TRIANGLES, gFaces[i].size() * 3, GL_UNSIGNED_INT, 0);
}

void renderText(const std::string& text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color)
{
    // Activate corresponding render state	
    glUseProgram(gProgram[3]);
    glUniform3f(glGetUniformLocation(gProgram[3], "textColor"), color.x, color.y, color.z);
    glActiveTexture(GL_TEXTURE0);
	glBindVertexArray(gTextVAO);

    // Iterate through all characters
    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++) 
    {
        Character ch = Characters[*c];

        GLfloat xpos = x + ch.Bearing.x * scale;
        GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

        GLfloat w = ch.Size.x * scale;
        GLfloat h = ch.Size.y * scale;

        // Update VBO for each character
        GLfloat vertices[6][4] = {
            { xpos,     ypos + h,   0.0, 0.0 },            
            { xpos,     ypos,       0.0, 1.0 },
            { xpos + w, ypos,       1.0, 1.0 },

            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos + w, ypos,       1.0, 1.0 },
            { xpos + w, ypos + h,   1.0, 0.0 }           
        };

        // Render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.TextureID);

        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

        //glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)

        x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }

    glBindTexture(GL_TEXTURE_2D, 0);
}


void display()
{
    glClearColor(0, 0, 0, 1);
	glClearDepth(1.0f);
	glClearStencil(0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);



	float knocked = glm::radians(0.0f);

	if(move_left && bunny_x > -2.5f){
		bunny_x -= 0.02f * speed;
	}
	if(move_right && bunny_x < 2.5f){
		bunny_x += 0.02f * speed;
	}

	if(hit_red){
		knocked = glm::radians(-90.0f);
		score_color = glm::vec3(1.f, 0.f, 0.f);


	}
	else{
		speed += 0.003f;
		bunnySpeed += 0.0005f;
		bunny_y += 0.01f;
		groundPos += 0.05f;
		box_position += 0.05f;
		score += 0.05 * speed;
	}

	// BUNNY IMPLEMENTATION
	float rotationAngle = glm::radians(270.0f);


	if(bunny_yay && rotation_y < 360.0f){
		rotation_y += 2.f * speed;
	}
	else{
		bunny_yay = false;
		rotation_y = 0.0f;
	}

	glm::mat4 transformMatrix = glm::mat4(1.0);
	// Compute the modeling matrix for the bunny
    transformMatrix = glm::translate(transformMatrix, glm::vec3(bunny_x, 0.3f * sin(10 * bunny_y * bunnySpeed) - 1.f, -3.f));
	// Add rotation to the bunny's modeling matrix with 270 degrees rotation around the z axis
	transformMatrix = glm::rotate(transformMatrix, glm::radians(rotation_y), glm::vec3(0, 1, 0));
	// Add rotation to the bunny's modeling matrix with 270 degrees rotation around the y axis
	transformMatrix = glm::rotate(transformMatrix, rotationAngle, glm::vec3(0, 1, 0));
	// Add rotation to the bunny's modeling matrix with 270 degrees rotation around the x axis
	transformMatrix = glm::rotate(transformMatrix, knocked, glm::vec3(1, 0, 0));
	// Scale the bunny to 1/10 of its original size
	transformMatrix = glm::scale(transformMatrix, glm::vec3(0.25f, 0.25f, 0.25f));

    // Set the active program and the values of its uniform variables for the bunny
    glUseProgram(gProgram[0]);
    glUniformMatrix4fv(projectionMatrixLoc[0], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
    glUniformMatrix4fv(viewingMatrixLoc[0], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
    glUniformMatrix4fv(modelingMatrixLoc[0], 1, GL_FALSE, glm::value_ptr(transformMatrix));
    glUniform3fv(eyePosLoc[0], 1, glm::value_ptr(eyePos));

    // Draw Bunny
    drawModel(0);


	// GROUND IMPLEMENTATION
	// Compute the modeling matrix for the ground

	float groundRotationAngle = glm::radians(-90.0f);
	glm::vec3 groundTranslation(0.0f, 0.0f, -2.0f);
	glm::mat4 groundTransformMatrix = glm::mat4(1.0);
	groundTransformMatrix = glm::rotate(groundTransformMatrix, groundRotationAngle, glm::vec3(1, 0, 0));
	groundTransformMatrix = glm::translate(groundTransformMatrix, groundTranslation);
	groundTransformMatrix = glm::scale(groundTransformMatrix, glm::vec3(5.f, 1000.0f, 1.0f)); // Adjust scale as needed


	// Set the active program and the values of its uniform variables for the ground quad
	glUseProgram(gProgram[1]); // Use the program for the ground (change the index if needed)

	// Set uniform values for the checkerboard shader
	glUniform3f(glGetUniformLocation(gProgram[1], "color1"), 1.f, 1.f, 1.f); // White
	glUniform3f(glGetUniformLocation(gProgram[1], "color2"), 0.f, 0.f, 0.f); // white
	glUniform1f(glGetUniformLocation(gProgram[1], "scale"), .5f); // Adjust scale as needed
	glUniform1f(glGetUniformLocation(gProgram[1], "offset"), 5.f); // Adjust width as needed
	glUniform1f(glGetUniformLocation(gProgram[1], "offsetZ"), groundPos * speed); // Adjust width as needed

	// Set the modeling matrix, viewing matrix, projection matrix, and eye position as needed
	glUniformMatrix4fv(projectionMatrixLoc[1], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
	glUniformMatrix4fv(viewingMatrixLoc[1], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
	glUniformMatrix4fv(modelingMatrixLoc[1], 1, GL_FALSE, glm::value_ptr(groundTransformMatrix));
	glUniform3fv(eyePosLoc[1], 1, glm::value_ptr(eyePos));

	// Draw the ground quad
	drawModel(1); // Use the appropriate index for the ground quad



	for(int j = 0; j < 3; j++){
		if(-40.f * reset + box_position * speed >= eyePos.z){
			my_rand = int(elapsedTime * 10) % 3;
			hit_yellow = false;
			hit_red = false;
			reset += 1;
			j = -1;
		}

		if((hit_red || hit_yellow) && j == hit_box){
			continue;
		}

		// Compare bunny position with box position z
		if(-3.f < -40.f * reset + box_position * speed + 1.5f && -3.f > -40.f * reset + box_position * speed - 1.5f){
			if(bunny_x < -2.f + j * 2.f + 0.5f && bunny_x > -2.f + j * 2.f - 0.5f){
				if(j == my_rand){
					hit_yellow = true;
					bunny_yay = true;
					score += 1000;					// INCREASE SCORE HERE
				}
				else{
					move_left = false;
					move_right = false;
					hit_red = true;
				}
				
				hit_box = j;
			}
			
		}

		

		glm::mat4 quadTransformationMatrix = glm::mat4(1.0);
		// Add rotation to the bunny's modeling matrix with 270 degrees rotation around the y axis
		quadTransformationMatrix = glm::rotate(quadTransformationMatrix, glm::radians(0.0f), glm::vec3(1, 0, 0));
		// Compute the modeling matrix for the bunny
		quadTransformationMatrix = glm::translate(quadTransformationMatrix, glm::vec3(-2.5f + j * 2.5f, -1.f, -40.f * reset + box_position * speed));
		// Scale the bunny to 1/2 of its original size
		quadTransformationMatrix = glm::scale(quadTransformationMatrix, glm::vec3(0.5f, 1.f, 0.5f));

		glUseProgram(gProgram[2]);
		glUniformMatrix4fv(projectionMatrixLoc[2], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
		glUniformMatrix4fv(viewingMatrixLoc[2], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
		glUniformMatrix4fv(modelingMatrixLoc[2], 1, GL_FALSE, glm::value_ptr(quadTransformationMatrix));
		glUniform3fv(eyePosLoc[2], 1, glm::value_ptr(eyePos));
		 if(j == my_rand)
		 	glUniform3f(glGetUniformLocation(gProgram[2], "kd"), 1.f, 1.f, 0.f); 
		 else
		 	glUniform3f(glGetUniformLocation(gProgram[2], "kd"), 1.f, 0.f, 0.f); 
		 //Draw Cube
		drawModel(2);
	}
	// CUBE IMPLEMENTATION
	
	//assert(glGetError() == GL_NO_ERROR);

	renderText("SCORE: " + std::to_string(int(score)), 0.f, 760.0f, 1.f, score_color);


	//assert(glGetError() == GL_NO_ERROR);


}

void reshape(GLFWwindow* window, int w, int h)
{
	w = w < 1 ? 1 : w;
	h = h < 1 ? 1 : h;

	gWidth = w;
	gHeight = h;

	glViewport(0, 0, w, h);

	// Use perspective projection
	float fovyRad = (float)(90.0 / 180.0) * M_PI;
	projectionMatrix = glm::perspective(fovyRad, w / (float)h, 1.0f, 100.0f);

	// Assume default camera position and orientation (camera is at
	// (0, 0, 0) with looking at -z direction and its up vector pointing
	// at +y direction)
	// 
	//viewingMatrix = glm::mat4(1);
	viewingMatrix = glm::lookAt(glm::vec3(0, 0, 0), glm::vec3(0, 0, 0) + glm::vec3(0, 0, -1), glm::vec3(0, 1, 0));

}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_Q && action == GLFW_PRESS)
	{
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	}
	else if (key == GLFW_KEY_A && action == GLFW_PRESS)
	{	if(!hit_red){
			move_left = true;
		}
	}
	else if (key == GLFW_KEY_A && action == GLFW_RELEASE)
	{	
		move_left = false;
	}
	else if (key == GLFW_KEY_D && action == GLFW_PRESS)
	{	if(!hit_red){
			move_right = true;
		}
	}
	else if (key == GLFW_KEY_D && action == GLFW_RELEASE)
	{
		move_right = false;
	}
	else if (key == GLFW_KEY_R && (action == GLFW_PRESS || action == GLFW_REPEAT))
	{
		restart = true;
	}
}

void mainLoop(GLFWwindow* window)
{
	while (!glfwWindowShouldClose(window))
	{
		if(restart){
			restart = false;

			score = 0;

			elapsedTime = 0.0f;

			speed = 1.0f; // Adjust this for the desired speed
			bunnySpeed = 1.0f; // Adjust this for the desired speed
			groundZMovement = 0.0f; // This will be used to move the ground quad

			box_position = 0.0f;
			reset = 0.0f;
			groundPos = 0.0f;

			bunny_x = 0;
			bunny_y = 0;

			move_right = false;
			move_left = false;

			hit_yellow = false;
			hit_red = false;
			hit_box = -1;
			bunny_yay = false;
			rotation_y = 0.0f;

			score_color = glm::vec3(1.f, 1.f, 0.f);
			
		}
		float currentTime = glfwGetTime();
        float deltaTime = currentTime - elapsedTime;
        elapsedTime = currentTime;
		display();
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
}

int main(int argc, char** argv)   // Create Main Function For Bringing It All Together
{
	GLFWwindow* window;
	if (!glfwInit())
	{
		exit(-1);
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this if on MacOS


	window = glfwCreateWindow(gWidth, gHeight, "Simple Example", NULL, NULL);

	if (!window)
	{
		glfwTerminate();
		exit(-1);
	}

	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	// Initialize GLEW to setup the OpenGL Function pointers
	if (GLEW_OK != glewInit())
	{
		std::cout << "Failed to initialize GLEW" << std::endl;
		return EXIT_FAILURE;
	}

	char rendererInfo[512] = { 0 };
	strcpy(rendererInfo, (const char*)glGetString(GL_RENDERER)); // Use strcpy_s on Windows, strcpy on Linux
	strcat(rendererInfo, " - "); // Use strcpy_s on Windows, strcpy on Linux
	strcat(rendererInfo, (const char*)glGetString(GL_VERSION)); // Use strcpy_s on Windows, strcpy on Linux
	glfwSetWindowTitle(window, rendererInfo);

	init();

	glfwSetKeyCallback(window, keyboard);
	glfwSetWindowSizeCallback(window, reshape);

	reshape(window, gWidth, gHeight); // need to call this once ourselves
	mainLoop(window); // this does not return unless the window is closed

	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
