#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

// Image resolution
int g_width;
int g_height;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

// TODO: add structs for spheres, lights and anything else you may need.
// DONE

struct Sphere
{
    string name;
    
    // Position components
    float p_x;
    float p_y;
    float p_z;
    
    // Scaling components
    float s_x;
    float s_y;
    float s_z;
    
    // Color components
    float c_r;
    float c_g;
    float c_b;
    
    // K - components
    float k_a;
    float k_d;
    float k_s;
    float k_r;
    
    // Specular exponent
    float n;
    
    mat4 M;
    mat4 M_inverse;
    
    bool invertible;
    
};

struct Light
{
    string name;
    
    // Position components
    float p_x;
    float p_y;
    float p_z;
    
    // Intensity components
    float I_r;
    float I_g;
    float I_b;
};

vector<vec4> g_colors;
vector<Sphere> spheres;
vector<Light> lights;

string output_file;  // output file name

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;

// Backgroung color components
float bg_r;
float bg_g;
float bg_b;

// Ambient intensity components
float amb_r;
float amb_g;
float amb_b;

inline mat4 getSphereMatrix(Sphere s);

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
    //DONE
    int i;
    for (i = 0; i < vs.size(); i++)
    {
        if (vs[i] == "RES")
        {
            i++;
            g_width = (int)toFloat(vs[i++]);
            g_height = (int)toFloat(vs[i]);
            g_colors.resize(g_width * g_height);
        }
        else if (vs[i] == "NEAR")
        {
            i++;
            g_near = toFloat(vs[i]);
        }
        else if (vs[i] == "LEFT")
        {
            i++;
            g_left = toFloat(vs[i]);
        }
        else if (vs[i] == "RIGHT")
        {
            i++;
            g_right = toFloat(vs[i]);
        }
        else if (vs[i] == "BOTTOM")
        {
            i++;
            g_bottom = toFloat(vs[i]);
        }
        else if (vs[i] == "TOP")
        {
            i++;
            g_top = toFloat(vs[i]);
        }
        else if (vs[i] == "SPHERE")
        {
            i++;
            Sphere temp;
            
            temp.name = vs[i++];
            
            temp.p_x = toFloat(vs[i++]);
            temp.p_y = toFloat(vs[i++]);
            temp.p_z = toFloat(vs[i++]);
            
            temp.s_x = toFloat(vs[i++]);
            temp.s_y = toFloat(vs[i++]);
            temp.s_z = toFloat(vs[i++]);
            
            temp.c_r = toFloat(vs[i++]);
            temp.c_g = toFloat(vs[i++]);
            temp.c_b = toFloat(vs[i++]);
            
            temp.k_a = toFloat(vs[i++]);
            temp.k_d = toFloat(vs[i++]);
            temp.k_s = toFloat(vs[i++]);
            temp.k_r = toFloat(vs[i++]);
            
            temp.n   = toFloat(vs[i]);
            
            temp.M   = getSphereMatrix(temp);
            
            temp.invertible = InvertMatrix(temp.M, temp.M_inverse);
            
            spheres.push_back(temp);
        }
        else if (vs[i] == "LIGHT")
        {
            i++;
            Light temp;
            
            temp.name = vs[i++];
            
            temp.p_x = toFloat(vs[i++]);
            temp.p_y = toFloat(vs[i++]);
            temp.p_z = toFloat(vs[i++]);
            
            temp.I_r = toFloat(vs[i++]);
            temp.I_g = toFloat(vs[i++]);
            temp.I_b = toFloat(vs[i]);
            
            lights.push_back(temp);
        }
        else if (vs[i] == "BACK")
        {
            i++;
            bg_r = toFloat(vs[i++]);
            bg_g = toFloat(vs[i++]);
            bg_b = toFloat(vs[i]);
        }
        else if (vs[i] == "AMBIENT")
        {
            i++;
            amb_r = toFloat(vs[i++]);
            amb_g = toFloat(vs[i++]);
            amb_b = toFloat(vs[i]);
        }
        else if (vs[i] == "OUTPUT")
        {
            i++;
            output_file = vs[i];
        }
    }
}

void checkParser()
{
    cout << "PARSER CHECK" << endl;
    
    // Print planes / values
    cout << "Planes:" << endl;
    cout << "\tNear: " << g_near << endl;
    cout << "\tLeft: " << g_left << endl;
    cout << "\tRight: " << g_right << endl;
    cout << "\tBottom: " << g_bottom << endl;
    cout << "\tTop: " << g_top << endl;
    cout << endl;
    
    // Print res
    cout << "Res: " << g_width << " x " << g_height << endl;
    
    // Print Spheres
    cout << "Spheres (" << spheres.size() << "):" << endl;
    for (int i = 0; i < spheres.size(); i++)
    {
        Sphere temp = spheres[i];
        cout << temp.name << ": " << endl;
        cout << "\tPos(" << temp.p_x << ", " << temp.p_y << ", " << temp.p_z << ")" << endl;
        cout << "\tScale(" << temp.s_x << ", " << temp.s_y << ", " << temp.s_z << ")" << endl;
        cout << "\tColor(" << temp.c_r << ", " << temp.c_g << ", " << temp.c_b << ")" << endl;
        cout << "\tK(" << temp.k_a << ", " << temp.k_d << ", " << temp.k_s << ", " << temp.k_r << ")" << endl;
        cout << "\tn = " << temp.n << endl;
    }
    
    cout << "Lights (" << lights.size() << "):" << endl;
    for (int i = 0; i < lights.size(); i++)
    {
        Light temp = lights[i];
        cout << temp.name << ": " << endl;
        cout << "\tPos(" << temp.p_x << ", " << temp.p_y << ", " << temp.p_z << ")" << endl;
        cout << "\tIntenisty(" << temp.I_r << ", " << temp.I_g << ", " << temp.I_b << ")" << endl;
    }
    
    cout << "Background color: (" << bg_r << ", " << bg_g << ", " << bg_b << ")" << endl;
    cout << "Ambient intensity: (" << amb_r << ", " << amb_g << ", " << amb_b << ")" << endl;
    cout << "Output filename: " << output_file << endl;
    
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
//    checkParser();
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}

// Obtains transform matrix of sphere
inline
mat4 getSphereMatrix(Sphere s)
{
    mat4 s_matrix = mat4();                     // Create matrix
    s_matrix *= Translate(s.p_x, s.p_y, s.p_z); // Translate to correct position
    s_matrix *= Scale(s.s_x, s.s_y, s.s_z);     // Scale to appropriate size
    
    return s_matrix;
}

// Strips the 0 or 1 component of points and vectors respectively, returning
// a length=3 vector
inline
vec3 vec4ToVec3(vec4 v)
{
    return vec3(v.x, v.y, v.z);
}

// returns B^2 - AC, the determinate of a quadratic formula
inline
float quadraticDeterminate(float a, float b, float c)
{
    return (b * b) - (a * c);
}

inline
float magnitude(vec3 v)
{
    return sqrt((v.x * v.x) + (v.y * v.y) + (v.z + v.z));
}

inline
void printVec3(vec3 v)
{
    cout << "Vector: (" << v.x << ", " << v.y << ", " << v.z << ")" << endl;
}

vec4 vectorBetweenPoints(vec4 a, vec4 b)
{
    return normalize(vec4(b.x - a.x, b.y - a.y, b.z - a.z, 0.0f));
}

// -------------------------------------------------------------------
// Intersection routine

// TODO: add your ray-sphere intersection routine here.
// DONE

// intersect(ray, sphere)
// This function takes a ray and a sphere and returns the lowest t-value
// of their intersection.
// If there are none, will return -1;
float intersect(Ray ray, Sphere sphere)
{
    // Save ray components
    vec4 S = ray.origin;
    vec4 c = ray.dir;
    
    // Generate the M and M_inverse matrices for each sphere
    mat4 M = sphere.M;
    mat4 M_inverse = sphere.M_inverse;
    if (sphere.invertible == false)
        printf("Matrix for sphere %s not invertible", sphere.name.c_str());
    
    // create inverse transform ray
    Ray inv_trans_ray;
    inv_trans_ray.origin = (M_inverse * S);
    inv_trans_ray.dir = (M_inverse * c);
    vec3 S_prime = vec4ToVec3(inv_trans_ray.origin);
    vec3 c_prime = vec4ToVec3(inv_trans_ray.dir);
    
    // set up quadratic equation
    float q_a = pow(length(c_prime), 2.0f);
    float q_b = dot(S_prime, c_prime);
    float q_c = pow(length(S_prime), 2.0f) - 1;
    float q_d = quadraticDeterminate(q_a, q_b, q_c);
    
    
    if (q_d < 0) // no solution
        return -1;
    else    // solution exists
    {
        float t_1 = -(q_b / q_a) + sqrt(q_d) / q_a;
        float t_2 = -(q_b / q_a) - sqrt(q_d) / q_a;
        
        return (t_1 < t_2) ? t_1 : t_2;
    }
    
    return -1;
}

// -------------------------------------------------------------------
// Ray tracing

vec3 Illuminate(Ray ray, Sphere sphere, float t)
{
    vec4 hit_point = ray.origin + ray.dir * t;
    vec4 S_prime = sphere.M_inverse * ray.origin;
    vec4 c_prime = sphere.M_inverse * ray.dir;
    
    vec4 trans_hit_point = S_prime + c_prime * t;
    trans_hit_point.w = 0;  // turn the point into a vector
    vec4 trans_norm = normalize(trans_hit_point);   // Normalize to a unit vector
    
    vec4 norm = normalize(sphere.M_inverse * trans_norm);
    
    float specular_sum_r = 0;
    float diffuse_sum_r  = 0;
    float specular_sum_g = 0;
    float diffuse_sum_g  = 0;
    float specular_sum_b = 0;
    float diffuse_sum_b  = 0;
//    float 
    
    for (int i = 0; i < lights.size(); i++)
    {
        Light t_light = lights[i];
        vec4 L = vectorBetweenPoints(hit_point, vec4(t_light.p_x, t_light.p_y, t_light.p_z, 1.0f));
        vec4 R = normalize(2 * norm * dot(norm, L) - L);
        
        float spec = pow(dot(R, ray.dir), sphere.n);
        float diff = dot(norm, L);
        
        specular_sum_r += t_light.I_r * spec;
        diffuse_sum_r  += t_light.I_r * diff;
        specular_sum_g += t_light.I_g * spec;
        diffuse_sum_g  += t_light.I_g * diff;
        specular_sum_b += t_light.I_b * spec;
        diffuse_sum_b  += t_light.I_b * diff;
    }
    
    diffuse_sum_r  *= sphere.k_d;
    diffuse_sum_g  *= sphere.k_d;
    diffuse_sum_b  *= sphere.k_d;
    specular_sum_r *= sphere.k_s;
    specular_sum_g *= sphere.k_s;
    specular_sum_b *= sphere.k_s;
    
    return vec3(diffuse_sum_r + specular_sum_r + amb_r * sphere.k_a,
                diffuse_sum_g + specular_sum_g + amb_g * sphere.k_a,
                diffuse_sum_b + specular_sum_b + amb_b * sphere.k_a);
}

vec4 trace(const Ray& ray)
{
    // TODO: implement your ray tracing routine here.
    // DONE (No illumnination)
    
    vector<float> t_values;
    
    // For each sphere...
    for (int i = 0; i < spheres.size(); i++)
        t_values.push_back(intersect(ray, spheres[i]));
    
    // find the value and index of lowest t value
    int lowest_index = -1;
    float lowest_t = numeric_limits<float>::max();
    for (int i = 0; i < spheres.size(); i++)
    {
        if (t_values[i] > 1 && t_values[i] < lowest_t)
        {
            lowest_t = t_values[i];
            lowest_index = i;
        }
    }
    
    if (lowest_index >= 0)
    {
        Sphere temp = spheres[lowest_index];
        vec3 intensity = Illuminate(ray, temp, lowest_t);
        return vec4(temp.c_r * intensity.x, temp.c_g * intensity.y, temp.c_b * intensity.z, 1.0f);
    }
        
    return vec4(bg_r, bg_g, bg_b, 1.0f); // default: set the color to the color of the background
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    // DONE
    
    vec4 dir;
    
    // Calculate individual components
    float x = g_left + ((float)ix / g_width) * (g_right - g_left);
    float y = g_bottom + ((float)iy / g_height) * (g_top - g_bottom);
    float z = -g_near;
    
    // Calculate resulting vector's magnitude
    float mag = sqrt((x * x) + (y * y) + (z * z));
    
    // Normalize during creation
    dir = vec4(x / mag, y / mag, z / mag, 0.0f);
    
    // Testing getDir()
//    cout << "Pixel: (" << ix << ", " << iy << "); Dir(" << dir.x << ", " << dir.y << ", " << dir.z << ")" << endl;
    
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

void savePPM(int Width, int Height, const char* fname, unsigned char* pixels)
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
    // DONE
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
            {
                float c = ((float*)g_colors[y*g_width+x])[i];
                if (c > 1)
                    c = 1;
                buf[y*g_width*3+x*3+i] = (unsigned char)(c * 255.0f);
            }
    // TODO: change file name based on input file name.
    // DONE
    savePPM(g_width, g_height, output_file.c_str(), buf);
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
	return 0;
}