
//Set the global scene parameter variables
//TODO: Set the scene parameters based on the values in the scene file

#ifndef PARSE_PGA_H
#define PARSE_PGA_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include "PGA_3D.h"
//Camera & Scene Parmaters (Global Variables)
class Ray {
public:
    Ray(Point3D p, Line3D l, Dir3D dir)
    {
        start = p;
        shoot = l.normalized();
        this->dir = dir;
    }
    Ray(Point3D p, Dir3D dir)
    {
        start = p;
        this->dir = dir;
        shoot = vee(p, dir).normalized();
    }
    Point3D start;
    Line3D shoot;
    Dir3D dir;
};
class Hit {
public:
    Dir3D N;
    Dir3D V;
    Dir3D R;
    Dir3D L;
    float intensity;
    float t; // distance
    Point3D interacts;
    bool hit;
    Hit()
    {
        hit = false;
    }
    Hit(Dir3D N, Dir3D V, Dir3D R, Dir3D L,float intensity,Point3D interacts)
    {
        this->N = N.normalized();
        this->V = V.normalized();
        this->R = R.normalized();
        this->L = L.normalized();
        this->intensity = intensity;
        this->interacts = interacts;
    }
};
class PointLights {
public:
    Color light;
    Point3D center;
    PointLights(Color color, Point3D position)
    {
        light = color;
        center = position;
    }
};
// material class
class Material {
public:
    float ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, ior;
    Material()
    {
        ar = ag = ab = dr = dg = db = sr = sg = sb = ns = tr = tg = tb = ior = 1;
    }
    Material(float a, float b, float c, float d, float e, float f, float g, float h,float i, float j, float k, float l, float m, float n)
    {
        ar = a;
        ag = b;
        ab = c;
        dr = d;
        dg = e;
        db = f;
        sr = g;
        sg = h;
        sb = i;
        ns = j;
        tr = k;
        tg = l;
        tb = m;
        ior = n;
    }
    
};
class Triangles
{
public:
    Point3D a;
    Point3D b;
    Point3D c;
    Dir3D norm;
    Material material;
    Triangles(Point3D a, Point3D b, Point3D c, Dir3D norm,Material material)
    {
        this->a = a;
        this->b = b;
        this->c = c;
        this->norm = norm;
        this->material = material;
    }
};
class Boxes
{

};
class Planes
{

};
class DirectLights
{
public:
    Color light;
    Dir3D dir;
    DirectLights(Color color, Dir3D position)
    {
        light = color;
        dir = position;
    }
};
class Spheres {
public:
    Spheres(Point3D center, float radius, Material mat)
    {
        pos = center;
        this->radius = radius;
        material = mat;
    }
    Point3D pos = Point3D(0, 0, 2);
    float radius = 1;
    Material material;

};
class Scene
{
public:
    std::vector<Spheres> spheres;
    std::vector<Triangles> triangles;
    std::vector<Boxes> boxes;
    std::vector<Planes> planes;
    std::vector<PointLights> pointLights;
    std::vector<DirectLights> directLights;

    //Scene Geometry
    int max_vextices;
    int max_normals;
    std::vector<Point3D> vertex;
    std::vector<Dir3D> normal;
};

//Image Parmaters
int img_width = 800, img_height = 600;
std::string imgName = "raytraced.png";

//Camera Parmaters
Point3D eye = Point3D(0,0,0); 
Dir3D forward = Dir3D(0,0,-1).normalized();
Dir3D up = Dir3D(0,1,0).normalized();
Dir3D right = Dir3D(-1,0,0).normalized();
float halfAngleVFOV = 35; 

//Scene (Sphere) Parmaters
Scene scene;
Point3D spherePos = Point3D(0,0,2);
float sphereRadius = 1; 
std::vector<Spheres> spheres;


Point3D origin = Point3D(0, 0, 0);

//background
Color backgroundcolor = Color(1,1,1);

//material parameters
Material CurrentMaterial = Material();


//lighting parameters
Color ambient = Color(0, 0, 0);

// sampling
int sampling = -1;

// focus parameter
float disk_radius = 0;  //default no focus 
float focus_length = 500;// img_height / 2 / tanf(halfAngleVFOV * (M_PI / 180.0f));

// max depth
int max_depth = 0; // 0 for debug
void parseSceneFile(std::string fileName){
  //TODO: Override the default values with new data from the file "fileName"
    // open the file containing the scene description
    std::ifstream input(fileName.c_str());
    std::string line;
    //scene.vertex.push_back(origin);
    // check for errors in opening the file
    if (input.fail()) {
        std::cout << "Can't open file '" << fileName << "'" << std::endl;
        return;
    }

    // determine the file size (this is optional -- feel free to delete the 6 lines below)
    std::streampos begin, end;
    begin = input.tellg();
    input.seekg(0, std::ios::end);
    end = input.tellg();
    std::cout << "File '" << fileName << "' is: " << (end - begin) << " bytes long.\n\n";
    input.seekg(0, std::ios::beg);


    //Loop through reading each line
    std::string command;
    while (input >> command) { //Read first word in the line (i.e., the command type)

        if (command[0] == '#') {
            getline(input, line); //skip rest of line
            std::cout << "Skipping comment: " << command << line << std::endl;
            continue;
        }
        if (command == "sphere:") { //If the command is a sphere command
            float x = 0, y = 0, z = 0, r = 2;
            input >> x >> y >> z >> r;
            printf("Sphere as position (%f,%f,%f) with radius %f\n", x, y, z, r);
            spherePos = Point3D(x, y, z);
            sphereRadius = r;
            Spheres newsphere = Spheres(spherePos, sphereRadius, CurrentMaterial);
            scene.spheres.push_back(newsphere);
        }
        else if (command == "image_resolution:") { //If the command is image_resolution
            int w = 800, h = 600;
            input >> w >> h;
            printf("Imgae resolution of (%d,%d)\n", w, h);
            img_width = w;
            img_height = h;
        }
        else if (command == "output_image:") { //If the command is output_image: 
            std::string outputname = "raytraced.png";
            input >> outputname;
            imgName = outputname;
            std::cout << "Output file name: " << outputname << std::endl;
        }
        else if (command == "camera_pos:") { //If the command is camera_pos:
            float x, y, z;
            x = y = z = 0;
            input >> x >> y >> z;
            printf("Camera position: (%f,%f,%f)\n", x, y, z);
            eye = Point3D(x, y, z);
        }
        else if (command == "camera_fwd:") { //If the command is camera_pos:
            float fx, fy;
            fx = fy = 0;
            float fz = -1;
            input >> fx >> fy >> fz;
            printf("camera fwd: (%f,%f,%f)\n", fx, fy, fz);
            forward = Dir3D(fx, fy, fz);
        }
        else if (command == "camera_up:") { //If the command is an camera_up: ux uy uz
            float ux, uy;
            ux = uy = 0;
            float uz = -1;
            input >> ux >> uy >> uz;
            printf("camera_up: (%f,%f,%f)\n", ux, uy, uz);
            up = Dir3D(ux, uy, uz);
        }
        else if (command == "camera_fov_ha:") { //If the command is an camera_fov_ha: ha
            float hf = 35;
            input >> hf;
            printf("camera_fov_ha: %f\n", hf);
            halfAngleVFOV = hf;
        }
        else if (command == "max_vertices:") { //If the command is an max vextices
            input >> scene.max_vextices;
            printf("max_vextices: %f\n", scene.max_vextices);
        }
        else if (command == "max_normals:") { //If the command is an max normals
            input >> scene.max_normals;
            printf("max_normals: %f\n", scene.max_normals);
        }
        else if (command == "vertex:") { //If the command is a vertex
            float x, y, z;
            input >> x>>y>>z;
            scene.vertex.push_back(Point3D(x, y, z));
            printf("Add vertex: %f %f %f\n", x,y,z);
        }
        else if (command == "normal:") { //If the command is a normal
            float x, y, z;
            input >> x >> y >> z;
            Dir3D temp = Dir3D(x, y, z).normalized();
            scene.normal.push_back(temp);
            printf("Add normal: %f %f %f\n", x, y, z);
        }
        else if (command == "triangle:") { //If the command is a triangle
            int x, y, z;
            input >> x >> y >> z;
            Dir3D norm = cross(scene.vertex[z] - scene.vertex[y], scene.vertex[x] - scene.vertex[y]).normalized(); // direction of the norm???
            Triangles temp(scene.vertex[x], scene.vertex[y], scene.vertex[z], norm, CurrentMaterial);
            scene.triangles.push_back(temp);
            printf("Add triangle: %d %d %d\n", x, y, z);
        }
        else if (command == "background:") { //If the command is a background
            float x, y, z;
            input >> x >> y >> z;
            backgroundcolor =  Color(x, y, z);
            printf("Set background: %f %f %f\n", x, y, z);
        }
        else if (command == "material:") { //If the command is a material
            float ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, ior;
            input >> ar >> ag >> ab >> dr >> dg >> db >> sr >> sg >> sb >> ns >> tr >> tg >> tb >> ior;
            CurrentMaterial = Material(ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, ior);
            printf("New material: .....\n");
        }
        else if (command == "point_light:") { //If the command is a light
            float r, g, b, x, y, z;
            input >> r >> g >> b >> x >> y >> z;
            Color color = Color(r, g, b);
            Point3D center = Point3D(x, y, z);
            scene.pointLights.push_back(PointLights(color, center));
            printf("New Point light: .....\n");
        }
        else if (command == "directional_light:") { //If the command is a light
            float r, g, b, x, y, z;
            input >> r >> g >> b >> x >> y >> z;
            Color color = Color(r, g, b);
            Dir3D center = Dir3D(x, y, z);
            center = center.normalized();
            scene.directLights.push_back(DirectLights(color, center));
            printf("New Direction light: .....\n");
        }
        else if (command == "ambient_light:") { //If the command is a light
            float r, g, b;
            input >> r >> g >> b;
            ambient = Color(r, g, b);
            printf("New ambient light: .....\n");
        }
        else if (command == "sampling:") { //If the command is an sampling method
            std::string samplings = "basic";
            input >> samplings;
            if (samplings == "basic")
            {
                sampling = 2;
                printf("Basic sampling\n");
            }
            else if (samplings == "jittered")
            {
                sampling = 0;
                printf("Jittered sampling\n");
            }
            else if (samplings == "adaptive")
            {
                sampling = 1;
                printf("Adaptive sampling\n");
            }
        }
        else if (command == "max_depth:") { //If the command is max_depth
            input >> max_depth;
            printf("Max depth %d\n", max_depth);
        }
        else if (command == "focus:") {
            input >> disk_radius >> focus_length;
            printf("Disk radius %f  focus length  %f", disk_radius, focus_length);
        }
        else {
            getline(input, line); //skip rest of line
            std::cout << "WARNING. Unknow command: " << command << std::endl;
        }
    }
    up = up.normalized();
    forward = forward.normalized();
    right = cross(up,forward);// = forward.cross(up);
  //TODO: Create an orthagonal camera basis, based on the provided up and right vectors
  printf("Orthagonal Camera Basis:\n");
  forward.print("forward");
  right.print("right");
  up.print("up");
}

#endif