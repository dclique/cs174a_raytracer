//
// raytrace.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
using namespace std;

int g_width;
int g_height;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

// TODO: add structs for spheres, lights and anything else you may need.

vector<vec4> g_colors;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;

vec4 g_back;
vec4 g_ambient;
char* g_output;

struct Light
{
	string name;
	vec4 position;
	vec4 intensity;
};

vector<Light> g_lights;

struct Sphere
{
	string name;
	vec4 position;
	vec3 scale;
	vec4 color;
	float ka, kd, ks, kr;
	float n;

	mat4 sphereMatrix;
};

vector<Sphere> g_spheres;

mat4 getMatrixSphere(Sphere s);


// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    //TODO: add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.
    if (vs[0] == "RES")
    {
        g_width = (int)toFloat(vs[1]);
        g_height = (int)toFloat(vs[2]);
        g_colors.resize(g_width * g_height);
    }
	else if(vs[0] == "NEAR")
	{
		g_near = toFloat(vs[1]);
	}
	else if(vs[0] == "LEFT")
	{
		g_left = toFloat(vs[1]);
	}
	else if(vs[0] == "RIGHT")
	{
		g_right = toFloat(vs[1]);
	}
	else if(vs[0] == "BOTTOM")
	{
		g_bottom = toFloat(vs[1]);
	}
	else if(vs[0] == "TOP")
	{
		g_top = toFloat(vs[1]);
	}
	else if(vs[0] == "SPHERE")
	{
		Sphere sphere;
		string spherename = vs[1];
		float x = toFloat(vs[2]);
		float y = toFloat(vs[3]);
		float z = toFloat(vs[4]);
		float scalex = toFloat(vs[5]);
		float scaley = toFloat(vs[6]);
		float scalez = toFloat(vs[7]);
		float r = toFloat(vs[8]);
		float g = toFloat(vs[9]);
		float b = toFloat(vs[10]);
		float ka = toFloat(vs[11]);
		float kd = toFloat(vs[12]);
		float ks = toFloat(vs[13]);
		float kr = toFloat(vs[14]);
		float n = toFloat(vs[15]);
		sphere.name = spherename;
		sphere.position = vec4(x,y,z,1.0f);
		sphere.scale = vec3(scalex,scaley,scalez);
		sphere.color = vec4(r,g,b,0.0f);
		sphere.ka = ka;
		sphere.kd = kd;
		sphere.ks = ks;
		sphere.kr = kr;
		sphere.n = n;

		sphere.sphereMatrix = getMatrixSphere(sphere);
		g_spheres.push_back(sphere);
	}
	else if(vs[0] == "LIGHT")
	{
		Light light;
		string lightname = vs[1];
		float x = toFloat(vs[2]);
		float y = toFloat(vs[3]);
		float z = toFloat(vs[4]);
		float r = toFloat(vs[5]);
		float g = toFloat(vs[6]);
		float b = toFloat(vs[7]);
		light.name = lightname;
		light.position = vec4(x,y,z,1.0f);
		light.intensity = (r,g,b,1.0f);
		g_lights.push_back(light);
	}
	else if(vs[0] == "BACK")
	{
		float r = toFloat(vs[1]);
		float g = toFloat(vs[2]);
		float b = toFloat(vs[3]);
		g_back = vec4(r,g,b,0.0f);
	}
	else if(vs[0] == "AMBIENT")
	{
		float r = toFloat(vs[1]);
		float g = toFloat(vs[2]);
		float b = toFloat(vs[3]);
		g_ambient = vec4(r,g,b,0.0f);
	}
	else if(vs[0] == "OUTPUT")
	{
		string output = vs[1];
		char* output2 = (char *)malloc(output.size() + 1);
		memcpy(output2, output.c_str(), output.size() + 1);
		g_output = output2;
	}
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}

//calculates the transformation matrix for a sphere
mat4 getMatrixSphere(Sphere s)
{
    mat4 m = mat4();                     // Create matrix
    m *= Translate(s.position.x, s.position.y, s.position.z); // Translate to correct position
    m *= Scale(s.scale.x, s.scale.y, s.scale.z);     // Scale to appropriate size
    
    return m;
}

//calculates the normalized vector between two points
vec4 vectorBetweenPoints(vec4 a, vec4 b)
{
    return normalize(vec4(b.x - a.x, b.y - a.y, b.z - a.z, 0.0f));
}


// -------------------------------------------------------------------
// Intersection routine

// TODO: add your ray-sphere intersection routine here.

int Intersection(Ray ray, Sphere sphere)
{


	vec4 origin = ray.origin;
	vec4 dir = ray.dir;

	mat4 M = sphere.sphereMatrix;
	mat4 M_inverse = mat4(1.0f);
	InvertMatrix(M,M_inverse);

	Ray inverseRay;
	inverseRay.origin = M_inverse * origin;
	inverseRay.dir = M_inverse * dir;
	vec3 origin2 = vec3(inverseRay.origin.x,inverseRay.origin.y,inverseRay.origin.z);
	vec3 dir2 = vec3(inverseRay.dir.x,inverseRay.dir.y,inverseRay.dir.z);

	//quadratic calculations
	float a = pow(length(dir2), 2.0f);
    float b = dot(origin2, dir2);
    float c = pow(length(origin2), 2.0f) - 1;
    float det = b*b - a*c;

	if(det < 0) // no intersection
		return -1;
	else    // intersection exists
    {
        float d1 = -(b / a) + sqrt(det) / a;
        float d2 = -(b / a) - sqrt(det) / a;
        
        return (d1 < d2) ? d1 : d2;
    }

}

vec3 Illuminate(Ray ray, Sphere sphere, float t)
{
	mat4 m_inverse = mat4(1.0f);
	InvertMatrix(sphere.sphereMatrix, m_inverse);
    vec4 hit = ray.origin + ray.dir * t;
    vec4 origin = m_inverse * ray.origin;
    vec4 dir = m_inverse * ray.dir;
    
    vec4 hit2 = origin + dir * t;
    hit2.w = 0;
    vec4 norm_hit = normalize(hit2);
    
    vec4 norm = normalize(m_inverse * norm_hit);
    
    float specular_r = 0;
    float diffuse_r  = 0;
    float specular_g = 0;
    float diffuse_g  = 0;
    float specular_b = 0;
    float diffuse_b  = 0;
    
	//for each light
    for (int i = 0; i < g_lights.size(); i++)
    {
        Light light = g_lights[i];
        vec4 L = vectorBetweenPoints(hit, vec4(light.position.x, light.position.y, light.position.z, 1.0f));
        vec4 R = normalize(2 * norm * dot(norm, L) - L);
        
        float spec = pow(dot(R, ray.dir), sphere.n);
        float diff = dot(norm, L);
        
        specular_r += light.intensity.x * spec;
        diffuse_r  += light.intensity.x * diff;
        specular_g += light.intensity.y * spec;
        diffuse_g  += light.intensity.y * diff;
        specular_b += light.intensity.z * spec;
        diffuse_b  += light.intensity.z * diff;
    }
    
    diffuse_r  *= sphere.kd;
    diffuse_g  *= sphere.kd;
    diffuse_b  *= sphere.kd;
    specular_r *= sphere.ks;
    specular_g *= sphere.ks;
    specular_b *= sphere.ks;
    
    return vec3(diffuse_r + specular_r + g_ambient.x * sphere.ka,
                diffuse_g + specular_g + g_ambient.y * sphere.ka,
                diffuse_b + specular_b + g_ambient.z * sphere.ka);
}


// -------------------------------------------------------------------
// Ray tracing


vec4 trace(const Ray& ray)
{
    // TODO: implement your ray tracing routine here.
    

	//find the first point of intersection for the ray
    vector<float> floats;

    for (int i = 0; i < g_spheres.size(); i++)
        floats.push_back(Intersection(ray, g_spheres[i]));
    
    int lowestfloat = -1;
    float lowest = numeric_limits<float>::max();
    for (int i = 0; i < g_spheres.size(); i++)
    {
        if (floats[i] > 1 && floats[i] < lowest)
        {
            lowest = floats[i];
            lowestfloat = i;
        }
    }
    
	//if an intersection exists
    if (lowestfloat >= 0)
    {
        Sphere temp = g_spheres[lowestfloat];
        vec3 intensity = Illuminate(ray, temp, lowest);
        return vec4(temp.color.x * intensity.x, temp.color.y * intensity.y, temp.color.z * intensity.z, 1.0f);
    }
    
	//if nothing intersects with ray, return color of background
    return vec4(g_back.x, g_back.y, g_back.z, 1.0f); 
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    
    vec4 dir;
    
    float x = g_left + ((float)ix / g_width) * (g_right - g_left);
    float y = g_bottom + ((float)iy / g_height) * (g_top - g_bottom);
    float z = -g_near;
    
	dir = vec4(x,y,z,0.0f);
	normalize(dir);

    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // TODO: clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
            {
                float a = ((float*)g_colors[y * g_width+x])[i];
                if (a > 1)
                    a = 1;
                buf[y * g_width * 3+x* 3+i] = (unsigned char)(a * 255.0f);
            }
    // TODO: change file name based on input file name.
    savePPM(g_width, g_height, g_output, buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
	free(g_output);
	return 0;
}

