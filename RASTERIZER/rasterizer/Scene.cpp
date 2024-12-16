#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace std;

using namespace tinyxml2;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	command = "convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}


/*
	Helper functions.
*/

/*
	This function calculates the camera transformation matrix.
*/
Matrix4 camera_transform_matrix(Camera *camera){
	double T_matrix[4][4] = 
	{
		{1, 0, 0, -camera->position.x},
		{0, 1, 0, -camera->position.y},
		{0, 0, 1, -camera->position.z},
		{0, 0, 0, 1}
	};
	Matrix4 T(T_matrix);

	double R_matrix[4][4] = 
	{
		{camera->u.x, camera->u.y, camera->u.z, 0},
		{camera->v.x, camera->v.y, camera->v.z, 0},
		{camera->w.x, camera->w.y, camera->w.z, 0},
		{0, 0, 0, 1}
	};
	Matrix4 R(R_matrix);

	return multiplyMatrixWithMatrix(R, T);
}

/*
	This function calculates the perspective projection matrix. It gets called in the forwardRendereingPipeline
	depending on the projectionType variable.
*/
Matrix4 perspective_projection_matrix(Camera *camera){
	double PP_matrix[4][4] = 
	{
		{2 * camera->near/(camera->right - camera->left), 0, (camera->right + camera->left) / (camera->right - camera->left), 0},
		{0, 2 * camera->near/(camera->top - camera->bottom), (camera->top + camera->bottom) / (camera->top - camera->bottom), 0},
		{0, 0, -(camera->far + camera->near) / (camera->far - camera->near), -(2 * camera->far * camera->near) / (camera->far - camera->near)},
		{0, 0, -1, 0}
	};
	Matrix4 PP(PP_matrix);

	return PP;
}

/*
	This function calculates the ortographic projection matrix. It gets called in the forwardRendereingPipeline
	depending on the projectionType variable.
*/
Matrix4 orthographic_projection_matrix(Camera *camera){
	double OP_matrix[4][4] = 
	{
		{2 / (camera->right - camera->left), 0, 0, -(camera->right + camera->left) / (camera->right - camera->left)},
		{0, 2 / (camera->top - camera->bottom), 0, -(camera->top + camera->bottom) / (camera->top - camera->bottom)},
		{0, 0, -2 / (camera->far - camera->near), -(camera->far + camera->near) / (camera->far - camera->near)},
		{0, 0, 0, 1}
	};
	Matrix4 OP(OP_matrix);
	return OP;
}


/*
	This function calculates the viewport transformation matrix .
*/
Matrix4 viewport_transform_matrix(Camera *camera){
	double VT_matrix[4][4] = 
	{
		{camera->horRes / 2.0, 0, 0, (camera->horRes - 1) / 2.0},
		{0, camera->verRes / 2.0, 0, (camera->verRes - 1) / 2.0},
		{0, 0, 0.5, 0.5},
		{0, 0, 0, 1}
	};
	Matrix4 VT(VT_matrix);
	return VT;
}

/*
	This function calculates the modeling transformation matrix.
*/
Matrix4 modeling_transform_matrix(Mesh *mesh, const vector<Scaling *> &scalings, const vector<Rotation *> &rotations, const vector<Translation *> &translations){
	Matrix4 MT = getIdentityMatrix();

	for (int i = 0; i < mesh->numberOfTransformations; i++)
	{
		if(mesh->transformationTypes[i] == 't'){
			Translation *translation = translations[mesh->transformationIds[i] - 1];
			double fransformation_matrix[4][4] = 
			{
				{1, 0, 0, translation->tx},
				{0, 1, 0, translation->ty},
				{0, 0, 1, translation->tz},
				{0, 0, 0, 1}
			};
			Matrix4 T(fransformation_matrix);
			MT = multiplyMatrixWithMatrix(T, MT);
		}
			
		if(mesh->transformationTypes[i] == 's'){
			Scaling *scaling = scalings[mesh->transformationIds[i] - 1];
			double scaling_matrix[4][4] = 
			{
				{scaling->sx, 0, 0, 0},
				{0, scaling->sy, 0, 0},
				{0, 0, scaling->sz, 0},
				{0, 0, 0, 1}
			};
			Matrix4 S(scaling_matrix);
			MT = multiplyMatrixWithMatrix(S, MT);
		}

		if(mesh->transformationTypes[i] == 'r'){
			Rotation *rotation = rotations[mesh->transformationIds[i] - 1];

			Vec3 u = Vec3(rotation->ux, rotation->uy, rotation->uz, -1);
			Vec3 v, w;

			if(std::min(std::min(abs(rotation->ux), abs(rotation->uy)), abs(rotation->uz)) == abs(rotation->ux)){
				v = Vec3(0, -u.z, u.y, -1);
			}
			else if(std::min(std::min(abs(rotation->ux), abs(rotation->uy)), abs(rotation->uz)) == abs(rotation->uy)){
				v = Vec3(-u.z, 0, u.x, -1);
			}
			else if(std::min(std::min(abs(rotation->ux), abs(rotation->uy)), abs(rotation->uz)) == abs(rotation->uz)){
				v = Vec3(-u.y, u.x, 0, -1);
			}

			w = crossProductVec3(u, v);
			v = normalizeVec3(v);
			w = normalizeVec3(w);

			//Construct model matrix for rotation.
			double model_matrix[4][4] = 
			{
				{u.x, u.y, u.z, 0},
				{v.x, v.y, v.z, 0},
				{w.x, w.y, w.z, 0},
				{0, 0, 0, 1}
			};
			double transpose_model_matrix[4][4] = 
			{
				{u.x, v.x, w.x, 0},
				{u.y, v.y, w.y, 0},
				{u.z, v.z, w.z, 0},
				{0, 0, 0, 1}
			};

			//Construct the x_rotation_matrix.
			double x_rotation_matrix[4][4] = 
			{
				{1, 0, 0, 0},
				{0, cos(rotation->angle * M_PI / 180), -sin(rotation->angle * M_PI / 180), 0},
				{0, sin(rotation->angle * M_PI / 180), cos(rotation->angle * M_PI / 180), 0},
				{0, 0, 0, 1}
			};

			Matrix4 first_rotation_matrix = multiplyMatrixWithMatrix(Matrix4(x_rotation_matrix), Matrix4(model_matrix));
			Matrix4 second_rotation_matrix = multiplyMatrixWithMatrix(Matrix4(transpose_model_matrix), first_rotation_matrix);
			MT = multiplyMatrixWithMatrix(second_rotation_matrix, MT);
			
		}
	}

	return MT;
}



bool cull_triangle(Vec4 &vertex1, Vec4 &vertex2, Vec4 &vertex3){
	//Take the vertices and calculate the normal vector.
	Vec3 edge2_1 = Vec3(vertex2.x - vertex1.x, vertex2.y - vertex1.y, vertex2.z - vertex1.z, -1);
	Vec3 edge3_1 = Vec3(vertex3.x - vertex1.x, vertex3.y - vertex1.y, vertex3.z - vertex1.z, -1);

	Vec3 normal_vector = normalizeVec3(crossProductVec3(edge2_1, edge3_1));

	double dot_product = dotProductVec3(normal_vector, Vec3(vertex1.x, vertex1.y, vertex1.z, -1));

	return dot_product < 0;
}


double f_helper(double x_a, double y_a, double x_b, double y_b, double x, double y){
    return (x * (y_a - y_b)) + (y * (x_b - x_a)) + (x_a * y_b) - (y_a * x_b);
}

double fit_in(double num, double res){
	num = num >= 0 ? num : 0;
	num = num < res - 1 ? num : res - 1;
	return num;
}

void rasterize_triangle(Vec4 &viewport_v1, Vec4 &viewport_v2, Vec4 &viewport_v3, Color *color_v1, Color *color_v2, Color *color_v3, vector<vector<Color>> &image, int horRes, int verRes, vector<vector<double>> &depth){
	
	double x_max = max(viewport_v1.x, max(viewport_v2.x, viewport_v3.x));
	x_max = fit_in(x_max, horRes);
	double x_min = min(viewport_v1.x, min(viewport_v2.x, viewport_v3.x));
	x_min = fit_in(x_min, horRes);

	double y_max = max(viewport_v1.y, max(viewport_v2.y, viewport_v3.y));
	y_max = fit_in(y_max, verRes);
	double y_min = min(viewport_v1.y, min(viewport_v2.y, viewport_v3.y));
	y_min = fit_in(y_min, verRes);

	double f_23 = f_helper(viewport_v2.x, viewport_v2.y, viewport_v3.x, viewport_v3.y, viewport_v1.x, viewport_v1.y);
	double f_31 = f_helper(viewport_v3.x, viewport_v3.y, viewport_v1.x, viewport_v1.y, viewport_v2.x, viewport_v2.y);
	double f_12 = f_helper(viewport_v1.x, viewport_v1.y, viewport_v2.x, viewport_v2.y, viewport_v3.x, viewport_v3.y);
	double f_23_xy, f_31_xy, f_12_xy;
	double alpha, beta, gamma;

	for(int x = x_min; x <= x_max; x++){
		for(int y = y_min; y <= y_max; y++){
			f_23_xy = f_helper(viewport_v2.x, viewport_v2.y, viewport_v3.x, viewport_v3.y, x, y);
			f_31_xy = f_helper(viewport_v3.x, viewport_v3.y, viewport_v1.x, viewport_v1.y, x, y);
			f_12_xy = f_helper(viewport_v1.x, viewport_v1.y, viewport_v2.x, viewport_v2.y, x, y);

			alpha = f_23_xy / f_23;
			beta = f_31_xy / f_31;
			gamma = f_12_xy / f_12;

			if(alpha >= 0 && beta >= 0 && gamma >= 0){
				double depth_z = alpha * viewport_v1.z + beta * viewport_v2.z + gamma * viewport_v3.z;
				if(depth_z <= depth[x][y]){
					depth[x][y] = depth_z;
					Color color = (*color_v1) * alpha + (*color_v2) * beta + (*color_v3) * gamma;
					image[x][y] = color.round();			// round it.
				}
				
			}

		}
	}
}


void rasterize_line(Vec4 &v1, Vec4 &v2, Color &c1, Color &c2, vector< vector<Color> > &image){
	

	Color color_v1 = c1;
	Color color_v2 = c2;
	double dx = abs(v2.x - v1.x);
	double dy = abs(v2.y - v1.y);
	
	int step;
	Color color, dc;
	double d, x, y;

	// slope m = (0,1]
	if(dx >= dy){

		if(v2.x < v1.x){
			swap(v1, v2);
			swap(color_v1, color_v2);
		
		}

		if(v2.y >= v1.y){
			step = 1;
		}
		else{
			step = -1;
		}

		color = color_v1;
		d = (v1.y - v2.y) + 0.5 * step * (v2.x - v1.x);
		dc = (color_v2 - color_v1) / (v2.x - v1.x);
		y = v1.y;

		for(x = v1.x; x <= v2.x; x++){
			image[x][y] = color.round();			// round it.

			if(d * step < 0){
				y = y + step;
				d += (v1.y - v2.y) + step * (v2.x - v1.x);
			}
			
			else{
				d += (v1.y - v2.y);
			}
			color += dc;
		}
	}

	// slope m = (1, inf)
	else{
		if(v2.y < v1.y){
			swap(v1, v2);
			swap(color_v1, color_v2);
		}

		if(v2.x >= v1.x){
			step = 1;
		}
		else{
			step = -1;
		}

		color = color_v1;
		d = (v2.x - v1.x) + 0.5 * step * (v1.y - v2.y);
		dc = (color_v2 - color_v1) / (v2.y - v1.y);
		x = v1.x;

		for(y = v1.y; y <= v2.y; y++){
			image[x][y] = color.round();				// round it.

			if(d * step > 0){
				x = x + step;
				d += (v2.x - v1.x) + step * (v1.y - v2.y);
			}
			
			else{
				d += (v2.x - v1.x);
			}
			color += dc;
		}
	}

}


bool clip_line(Vec4 &vertex1, Vec4 &vertex2, Color &color1, Color &color2){
	//In this function we will use Liang-Barsky algorithm to clip the line and return true if the line is visible.
	//We used Liang-Barsky algorithm because it is more efficient(it is approximately %40 faster) than Cohen-Sutherland algorithm.

	double x_min = -1;
	double x_max = 1;
	double y_min = -1;
	double y_max = 1;
	double z_min = -1;
	double z_max = 1;

	double t_E = 0;
	double t_L = 1;

	double dx = vertex2.x - vertex1.x;
	double dy = vertex2.y - vertex1.y;
	double dz = vertex2.z - vertex1.z;

	Color dc = color2 - color1;

	double p[6] = {dx, -dx, dy, -dy, dz, -dz}; //Den

	double q[6] = {x_min - vertex1.x, vertex1.x - x_max, y_min - vertex1.y, vertex1.y - y_max, z_min - vertex1.z, vertex1.z - z_max}; //Num

	for(int i = 0; i < 6; i++){
		if(p[i] > 0){
			double t = q[i] / p[i];

			if(t > t_L){
				return false;
			}
			else if(t > t_E){
				t_E = t;
			}
		}
		else if(p[i] < 0){
			double t = q[i] / p[i];

			if(t < t_E){
				return false;
			}
			else if(t < t_L){
				t_L = t;
			}
		}
		else if(q[i] > 0){
			return false;
		}
	}

	if(t_L < 1){
		vertex2.x = vertex1.x + t_L * dx;
		vertex2.y = vertex1.y + t_L * dy;
		vertex2.z = vertex1.z + t_L * dz;
		color2 = color1 + dc * t_L;
	}

	if(t_E > 0){
		vertex1.x = vertex1.x + t_E * dx;
		vertex1.y = vertex1.y + t_E * dy;
		vertex1.z = vertex1.z + t_E * dz;
		color1 = color1 + dc * t_E;
	}
	

	return true;
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	Matrix4 camera_trans_matrix = camera_transform_matrix(camera); // Calculate Camera translation matrix.
	Matrix4 viewport_trans_matrix = viewport_transform_matrix(camera); // Calculate viewport translation matrix.

	Matrix4 projection_matrix;  // Calculate projection translation matrix.
	if(camera->projectionType){
		projection_matrix = perspective_projection_matrix(camera);
	}
	else{
		projection_matrix = orthographic_projection_matrix(camera);
	}

	// For each model apply translations, clipping, culling, rasterization.
	for(int i = 0; i < this->meshes.size(); i++){
		Mesh *mesh = this->meshes[i];

		Matrix4 modeling_matrix = modeling_transform_matrix(mesh, scalings, rotations, translations);
		Matrix4 camera_modeling_matrix = multiplyMatrixWithMatrix(camera_trans_matrix, modeling_matrix);
		Matrix4 camera_modeling_projection_matrix = multiplyMatrixWithMatrix(projection_matrix, camera_modeling_matrix);

		// For each triangle apply translations, clipping, culling, rasterization.
		for(int j = 0; j < mesh->triangles.size(); j++){
			Triangle triangle = mesh->triangles[j];

			Vec4 proj_cam_v1 = multiplyMatrixWithVec4(camera_modeling_projection_matrix, Vec4(vertices[triangle.vertexIds[0] - 1]->x,
																							 vertices[triangle.vertexIds[0] - 1]->y,
																							 vertices[triangle.vertexIds[0] - 1]->z, 1,
																							 vertices[triangle.vertexIds[0] - 1]->colorId));

			Vec4 proj_cam_v2 = multiplyMatrixWithVec4(camera_modeling_projection_matrix, Vec4(vertices[triangle.vertexIds[1] - 1]->x,
																							 vertices[triangle.vertexIds[1] - 1]->y,
																							 vertices[triangle.vertexIds[1] - 1]->z, 1,
																							 vertices[triangle.vertexIds[1] - 1]->colorId));

			Vec4 proj_cam_v3 = multiplyMatrixWithVec4(camera_modeling_projection_matrix, Vec4(vertices[triangle.vertexIds[2] - 1]->x,
																							 vertices[triangle.vertexIds[2] - 1]->y,
																							 vertices[triangle.vertexIds[2] - 1]->z, 1,
																							 vertices[triangle.vertexIds[2] - 1]->colorId));

			//Check culling if the culling is enabled.
			if(cullingEnabled){
				if(cull_triangle(proj_cam_v1, proj_cam_v2, proj_cam_v3)){
					continue;
				}
			}

			//Get the colors.
			Color *color_v1 = colorsOfVertices[vertices[triangle.vertexIds[0] - 1]->colorId - 1];
			Color *color_v2 = colorsOfVertices[vertices[triangle.vertexIds[1] - 1]->colorId - 1];
			Color *color_v3 = colorsOfVertices[vertices[triangle.vertexIds[2] - 1]->colorId - 1];

			//Check the mesh type.(Wireframe or solid)
			if(mesh->type){
				//This is the solid mesh.

				proj_cam_v1 = proj_cam_v1 / proj_cam_v1.t;
				proj_cam_v2 = proj_cam_v2 / proj_cam_v2.t;
				proj_cam_v3 = proj_cam_v3 / proj_cam_v3.t;

				Vec4 viewport_v1 = multiplyMatrixWithVec4(viewport_trans_matrix, proj_cam_v1);
				Vec4 viewport_v2 = multiplyMatrixWithVec4(viewport_trans_matrix, proj_cam_v2);
				Vec4 viewport_v3 = multiplyMatrixWithVec4(viewport_trans_matrix, proj_cam_v3);

				rasterize_triangle(viewport_v1, viewport_v2, viewport_v3, color_v1, color_v2, color_v3, this->image, camera->horRes, camera->verRes, depth);

			}

			else{
				//This is the wireframe mesh.
				proj_cam_v1 = proj_cam_v1 / proj_cam_v1.t;
				proj_cam_v2 = proj_cam_v2 / proj_cam_v2.t;
				proj_cam_v3 = proj_cam_v3 / proj_cam_v3.t;

				Vec4 v1_c1 = proj_cam_v1;
				Vec4 v2_c1 = proj_cam_v2;
				Vec4 v3_c1 = proj_cam_v3;
				Vec4 v1_c2 = proj_cam_v1;
				Vec4 v2_c2 = proj_cam_v2;
				Vec4 v3_c2 = proj_cam_v3;

				Color color1_c1 = Color(*color_v1);
				Color color2_c1 = Color(*color_v2);
				Color color3_c1 = Color(*color_v3);
				Color color1_c2 = Color(*color_v1);
				Color color2_c2 = Color(*color_v2);
				Color color3_c2 = Color(*color_v3);

				bool edge1_is_visible = clip_line(v1_c1, v2_c1, color1_c1, color2_c1);
				bool edge2_is_visible = clip_line(v2_c2, v3_c1, color2_c2, color3_c1);
				bool edge3_is_visible = clip_line(v3_c2, v1_c2, color3_c2, color1_c2);

				//Do viewport transformation to all vertices.

				v1_c1 = multiplyMatrixWithVec4(viewport_trans_matrix, v1_c1);
				v2_c1 = multiplyMatrixWithVec4(viewport_trans_matrix, v2_c1);
				v3_c1 = multiplyMatrixWithVec4(viewport_trans_matrix, v3_c1);

				v1_c2 = multiplyMatrixWithVec4(viewport_trans_matrix, v1_c2);
				v2_c2 = multiplyMatrixWithVec4(viewport_trans_matrix, v2_c2);
				v3_c2 = multiplyMatrixWithVec4(viewport_trans_matrix, v3_c2);

				//Do rasterization to all edges.
				if(edge1_is_visible){
					rasterize_line(v1_c1, v2_c1, color1_c1, color2_c1, this->image);
				}
				if(edge2_is_visible){
					rasterize_line(v2_c2, v3_c1, color2_c2, color3_c1, this->image);
				}
				if(edge3_is_visible){
					rasterize_line(v3_c2, v1_c2, color3_c2, color1_c2, this->image);
				}
			}
		}
	}
}
