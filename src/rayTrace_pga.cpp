//CSCI 5607 HW3 - Rays & Files
//This HW has three steps:
// 1. Compile and run the program (the program takes a single command line argument)
// 2. Understand the code in this file (rayTrace_pga.cpp), in particular be sure to understand:
//     -How ray-sphere intersection works
//     -How the rays are being generated
//     -The pipeline from rays, to intersection, to pixel color
//    After you finish this step, and understand the math, take the HW quiz on canvas
// 3. Update the file parse_pga.h so that the function parseSceneFile() reads the passed in file
//     and sets the relevant global variables for the rest of the code to product to correct image

//To Compile: g++ -fsanitize=address -std=c++11 rayTrace_pga.cpp

//For Visual Studios
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//Images Lib includes:
#define STB_IMAGE_IMPLEMENTATION //only place once in one .cpp file
#define STB_IMAGE_WRITE_IMPLEMENTATION //only place once in one .cpp files
#define M_PI 3.1415926
#include "image_lib.h" //Defines an image class and a color class

//#3D PGA
#include "PGA_3D.h"

//High resolution timer
#include <chrono>

//Scene file parser
#include "parse_pga.h"


Hit rayTriangleIntersect(Ray ray, Triangles triangle);
bool raySphereIntersect_fast(Point3D rayStart, Line3D rayLine, Point3D sphereCenter, float sphereRadius);
Hit raySphereIntersect(Point3D rayStart, Line3D rayLine, Point3D sphereCenter, float sphereRadius);
Hit lightSphereIntersect(Point3D rayStart, Line3D rayLine, Point3D sphereCenter, float sphereRadius, Point3D lightSource);
bool sameside(Point3D p1, Point3D p2, Point3D a, Point3D b);
bool inside_triangle(Point3D interact_point, Triangles triangle);
Color GetColor(Material material, Dir3D N, Dir3D L, Dir3D V, int Sl, float intensity, Dir3D R); // sum over point light term
Color AddColor(Color color, Color emissive, Material mat, Color ambient); // add on ambient term
Color RenderingTriangle(Scene scene, Ray ray, float& shortest, Color color,int depth);
Color Focus(Point3D eye, Point3D focus,int depth);
Color ColorOnThisPixel(Point3D eye, Ray views,int depth);
Color RenderingSphere(Scene scene, Ray views, float& shortest, Color color,int depth);
Color ReflectAndRefractContribtion(Color reflect, Color refract, Material mat);
Ray Refract(Ray ray, Hit hit, Material mat);
Ray Reflect(Ray ray, Dir3D norm, Point3D start);

int main(int argc, char** argv){
  
  //Read command line paramaters to get scene file
  if (argc != 2){
     std::cout << "Usage: ./a.out scenefile\n";
     return(0);
  }
  std::string secenFileName = argv[1];

  //Parse Scene File
  parseSceneFile(secenFileName);

  float imgW = img_width, imgH = img_height;
  float halfW = imgW/2, halfH = imgH/2;
  float d = halfH / tanf(halfAngleVFOV * (M_PI / 180.0f));
  // defualt d = 428.444
  // d = halfW / tanf(halfAngleVFOV * (M_PI / 180.0f));
  Color back = backgroundcolor;
  Image outputImg = Image(img_width,img_height);
  auto t_start = std::chrono::high_resolution_clock::now();


  for (int i = 0; i < img_width; i++){
      //std::cout << i << " "<<img_width<<std::endl;
    for (int j = 0; j < img_height; j++){
      //TODO: In what way does this assumes the basis is orthonormal?
        float u,v;
        if (sampling == 0) // jittered supersampling
        {
            float rand1 = (float)rand() / (RAND_MAX);
            float rand2 = (float)rand() / (RAND_MAX);
            u = (halfW - (imgW) * ((i + rand1) / imgW));
            v = (halfH - (imgH) * ((j + rand2) / imgH));
        }
        else if (sampling == 1) //adaptive supersampling
        {

        }
        else// basic sampling
        {
            u = (halfW - (imgW) * ((i + 0.5) / imgW));
            v = (halfH - (imgH) * ((j + 0.5) / imgH));
        }
      
      Point3D p = eye - d*forward + u*right + v*up;
      Dir3D rayDir = (p - eye); 
      Line3D rayLine = vee(eye,rayDir).normalized();  //Normalizing here is optional

      Point3D focus = p + (focus_length-d) * rayDir.normalized();
      Ray view(p, rayLine, rayDir.normalized());
      Color temp;
      if (disk_radius == 0)
      {
          temp = ColorOnThisPixel(eye, view, max_depth);
      } 
      else
      {
          temp = Focus(eye, focus, max_depth);
      }
      outputImg.setPixel(i, j, temp);
    }
  }
  auto t_end = std::chrono::high_resolution_clock::now();
  printf("Rendering took %.2f ms\n",std::chrono::duration<double, std::milli>(t_end-t_start).count());

  outputImg.write(imgName.c_str());
  return 0;
}

Color ColorOnThisPixel(Point3D eye, Ray views,int depth) // eye is the point on disk????, view is the ray connect eye and focus
{
    Color color = backgroundcolor;
    float shortest = INT_MAX;
    
    color = RenderingSphere(scene, views, shortest, color,depth);
    color = RenderingTriangle(scene, views, shortest, color,depth);
    
    return color;
}
Color Focus(Point3D eye,Point3D focus,int depth) // calculate the average color pixel 
{
    float x, y, z;
    float stepsize = 0.5;
    int count = disk_radius/stepsize*180; // how many average points
    
    Color average(0, 0, 0);
    for (float k = 0; k < disk_radius; k += stepsize)
    {
        for (float i = 0; i < M_PI; i += M_PI / 180) // i has step of 2 degree range from 0 to 2 PI
        {
            x = cos(i) * k;
            y = sin(i) * k;
            Point3D new_eye(eye.x + x, eye.y + y, eye.z); // this is the new point of eye( wondering around the original one) ( based on spherical coordinate)
            Dir3D dir = (focus - new_eye).normalized();
            Ray view(new_eye, dir);
            // add depth here
            Color result = ColorOnThisPixel(new_eye, view,depth);
            //std::cout << result.r << result.g << result.b << std::endl;
            average = average + result *(1.0/ count);
        }
        
    }
    
    //std::cout << average.r << average.g << average.b << std::endl;
    return average;
}
Color RenderingSphere(Scene scene, Ray views, float& shortest, Color color,int depth)
{
    for (int k = 0; k < scene.spheres.size(); k++) // traverse all spheres in the scene
    {
        Hit hits = raySphereIntersect(eye, views.shoot, scene.spheres[k].pos, scene.spheres[k].radius);
        if (hits.hit)
        {
            Dir3D view = eye - scene.spheres[k].pos;
            float distance = view.magnitude();
            if (distance < shortest) //find the spheres closest to viewer
            {
                shortest = distance;
                // add ambient response
                color = AddColor(Color(0, 0, 0), Color(0, 0, 0), scene.spheres[k].material, ambient); // add ambient light

                // add point light
                for (int l = 0; l < scene.pointLights.size(); l++) // sum over all point light
                {
                    int shadow = 1;
                    Dir3D L = scene.pointLights[l].center - hits.interacts;
                    float intensity = L.magnitude();
                    intensity = 1 / (0.25 * intensity);
                    L = L.normalized();
                    Dir3D R = L - 2 * (dot(L, hits.N) * hits.N);

                    Line3D rayLine1 = vee(hits.interacts, L).normalized();
                    for (int ki = 0; ki < scene.spheres.size(); ki++) // this loop check the shadow, if the light is visible, not block by others
                    {
                        if (ki != k) // don't check itself
                        {
                            if (raySphereIntersect(hits.interacts, rayLine1, scene.spheres[ki].pos, scene.spheres[ki].radius).hit) // if hit, generate shadow
                            {
                                shadow = 0;
                            }
                        }
                    }
                    color = color + GetColor(scene.spheres[k].material, hits.N, L, hits.V, shadow, intensity, R); //Sl ??                     
                }
                

                // add direct light
                for (int l = 0; l < scene.directLights.size(); l++) // sum over all direction light
                {
                    int shadow = 1;
                    for (int ki = 0; ki < scene.spheres.size(); ki++) // this loop check the shadow, if the light is visible, not block by others
                    {
                        if (ki != k) // don't check itself
                        {
                            Dir3D rayDir1 = scene.directLights[l].dir;
                            Line3D rayLine1 = vee(hits.interacts, rayDir1).normalized();
                            if (raySphereIntersect(hits.interacts, rayLine1, scene.spheres[ki].pos, scene.spheres[ki].radius).hit) // if hit, generate shadow
                            {
                                shadow = 0;
                            }
                        }
                    }
                    Dir3D R = scene.directLights[l].dir - 2 * dot(scene.directLights[l].dir, hits.N) * hits.N;
                    color = color + GetColor(scene.spheres[k].material, hits.N, scene.directLights[l].dir, hits.V, shadow, 3, R);
                }
                
                // add reflect and refract
                if (depth > 0)
                {
                    Ray mirror = Reflect(views, hits.N, hits.interacts); 
                    Ray glass = Refract(views, hits, scene.spheres[k].material);
                    
                    Color reflect = ColorOnThisPixel(hits.interacts + mirror.dir*0.01, mirror, depth - 1); // add a small factor so it will not interact with current sphere
                    Color refract = ColorOnThisPixel(hits.interacts + glass.dir * 0.01, glass, depth - 1);
                    color = color + ReflectAndRefractContribtion(reflect, refract, scene.spheres[k].material);
                }
            }
        }

    }
    return color;
}
Color RenderingTriangle(Scene scene, Ray ray, float& shortest,Color color,int depth)
{
    for (int k = 0; k < scene.triangles.size(); k++)
    {
        Hit hits = rayTriangleIntersect(ray, scene.triangles[k]);
        
        if (hits.hit && hits.t < shortest)
        {
            
            shortest = hits.t; 
            // rendering triangle
            
            // add ambient response
            color = AddColor(Color(0, 0, 0), Color(0, 0, 0), scene.triangles[k].material, ambient); // add ambient light

            //add light
            
            // sum over all point light
            for (int l = 0; l < scene.pointLights.size(); l++)
            {
                int shadow = 1;
                for (int ki = 0; ki < scene.spheres.size(); ki++) // this loop check the shadow, if the light is visible, not block by others
                {
                    Dir3D rayDir1 = (scene.pointLights[l].center - hits.interacts);
                    Line3D rayLine1 = vee(hits.interacts, rayDir1).normalized();
                    if (raySphereIntersect(hits.interacts, rayLine1, scene.spheres[ki].pos, scene.spheres[ki].radius).hit) // if hit, generate shadow
                    {
                        shadow = 0;
                    }
                }
                for (int ki = 0; ki < scene.triangles.size(); ki++) // this loop check the shadow, if the light is visible, not block by others
                {
                    if (ki != k)
                    {
                        Dir3D rayDir1 = (scene.pointLights[l].center - hits.interacts);
                        Line3D rayLine1 = vee(hits.interacts, rayDir1).normalized();
                        Ray shadows(hits.interacts, rayDir1.normalized());
                        if (rayTriangleIntersect(shadows, scene.triangles[ki]).hit) // if hit, generate shadow
                        {
                            shadow = 0;
                        }
                    }
                }
                Dir3D L = (scene.pointLights[l].center - hits.interacts).normalized();
                Dir3D R = L - 2 * dot(L, hits.N) * hits.N;
                color = color + GetColor(scene.triangles[k].material, hits.N, L, hits.V, shadow, hits.intensity, R); //Sl ??                     
            }
            
            // sum over all direction light
            for (int l = 0; l < scene.directLights.size(); l++) 
            {
                int shadow = 1;
                for (int ki = 0; ki < scene.spheres.size(); ki++) // this loop check the shadow, if the light is visible, not block by others
                {
                    Dir3D rayDir1 = scene.directLights[l].dir;
                    Line3D rayLine1 = vee(hits.interacts, rayDir1).normalized();
                    if (raySphereIntersect(hits.interacts, rayLine1, scene.spheres[ki].pos, scene.spheres[ki].radius).hit) // if hit, generate shadow
                    {
                        shadow = 0;
                    }
                }
                for (int ki = 0; ki < scene.triangles.size(); ki++) // this loop check the shadow, if the light is visible, not block by others
                {
                    
                    if (ki != k) // for triangle loop, don't check itself
                    {
                        
                        Dir3D rayDir1 = scene.directLights[l].dir;
                        Line3D rayLine1 = vee(hits.interacts, rayDir1).normalized();
                        Ray shadows(hits.interacts, rayDir1.normalized());
                        if (rayTriangleIntersect(shadows, scene.triangles[ki]).hit) // if hit, generate shadow
                        {
                            shadow = 0;
                        }
                    }  
                }
                Dir3D R = scene.directLights[l].dir - 2 * dot(scene.directLights[l].dir, hits.N) * hits.N;
                color = color + GetColor(scene.triangles[k].material, hits.N, scene.directLights[l].dir, hits.V, shadow, 3, R);
            }

            // add reflection and refraction
            if (depth > 0)
            {
                Ray mirror = Reflect(ray, hits.N, hits.interacts); // only reflect half light intensity
                Ray glass = Refract(ray, hits, scene.triangles[k].material);
                int depth1 = depth;
                Color reflect = ColorOnThisPixel(hits.interacts, mirror, depth - 1);
                Color refract = ColorOnThisPixel(hits.interacts, glass, depth1 - 1);
                color = color + ReflectAndRefractContribtion(reflect, refract, scene.triangles[k].material);
            }
        }
    }
    
    return color;
}
Color GetColor(Material material, Dir3D N, Dir3D L, Dir3D V, int Sl, float intensity, Dir3D R)
{
    float parameter1 = dot(N, L);
    if (parameter1 < 0)
    {
        parameter1 = 0;
    }
    float parameter2 = pow(dot(V, R), material.ns);
    float r = (material.dr * parameter1 + material.sr * parameter2) * Sl * intensity;
    float g = (material.dg * parameter1 + material.sg * parameter2) * Sl * intensity;
    float b = (material.db * parameter1 + material.sb * parameter2) * Sl * intensity;
    return Color(r, g, b);
}
Color AddColor(Color color, Color emissive, Material mat, Color ambient)
{
    float r = color.r + emissive.r + mat.ar*ambient.r;
    float g = color.g + emissive.g + mat.ag*ambient.g;
    float b = color.b + emissive.b + mat.ab*ambient.b;
    return Color(r, g, b);
}


bool sameside(Point3D p1, Point3D p2, Point3D a, Point3D b)
{
    Dir3D cp1 = cross(b - a, p1 - a);
    Dir3D cp2 = cross(b - a, p2 - a);
    return dot(cp1, cp2) >= 0;
}
bool inside_triangle(Point3D interact_point, Triangles triangle)
{
    return sameside(interact_point, triangle.a, triangle.b, triangle.c) && sameside(interact_point, triangle.b, triangle.a, triangle.c) && sameside(interact_point, triangle.c, triangle.a, triangle.b);
}
Hit rayTriangleIntersect(Ray ray, Triangles triangle) { // if ray and sphere intersect
    //assume triangle normal is pointing outside
    //Plane3D tri = wedge(triangle.a, triangle.b, triangle.c);
    Hit result;
    if (dot(triangle.norm, ray.dir) == 0)
    {
        result.hit = false;
        return result;
    }
    if (dot(triangle.norm, ray.dir) > 0)
    {
        triangle.norm = triangle.norm * (-1); // use the other norm
    }
    float D = -1*dot(triangle.norm, triangle.a);

    float t = -(dot(triangle.norm, ray.start) + D) / dot(triangle.norm, ray.dir);
    
    Point3D interact_point = ray.start + t * ray.dir;
    //interact_point.print();
    //std::cout << D << " " << t << std::endl;
    //interact_point.print();
    if (inside_triangle(interact_point, triangle))
    {
        //std::cout << "hit\n";
        result.interacts = interact_point;
        result.N = triangle.norm;
        result.t = t;
        result.hit = true;
        result.V = (interact_point - ray.start).normalized();
        result.intensity = 1 / (0.5 * t);
        return result;
    }
    else
    {
        //std::cout << "hit!!!!!!!!!!!!!!\n";
        result.hit = false;
        return result;
    }
}

Hit raySphereIntersect(Point3D rayStart, Line3D rayLine, Point3D sphereCenter, float sphereRadius) {
    Hit result;
    result.hit = false;
    Point3D projPoint = proj(sphereCenter, rayLine);
    float d1 = sphereCenter.distTo(projPoint);
    if (d1 > sphereRadius)
    {
        result.hit = false;
        return result;
    }
    float w = sqrt(sphereRadius * sphereRadius - d1 * d1);
    Point3D p1 = projPoint + rayLine.dir() * w;                   //Add/subtract above distance to find hit points
    Point3D p2 = projPoint - rayLine.dir() * w;
    if (dot((p2 - rayStart), rayLine.dir()) >= 0)
    {
        result.interacts = p2;
        result.hit = true;
        result.V = (p2 - rayStart).normalized();
        result.N = (p2 - sphereCenter).normalized();
        return result;     //Is the second point in same direction as the ray line?
    }
    if (dot((p1 - rayStart), rayLine.dir()) >= 0)
    {

        result.interacts = p1;
        result.hit = true;
        result.V = (p1 - rayStart).normalized();
        result.N = (p1 - sphereCenter).normalized();
        return result;     //Is the first point in same direction as the ray line?
    }
    

    //std::cout << "wrong";
    return result;
}
Color ReflectAndRefractContribtion(Color reflect, Color refract, Material mat)
{
    float red = reflect.r * mat.sr + refract.r * mat.tr;
    float green = reflect.g * mat.sg + refract.g * mat.tg;
    float blue = reflect.b * mat.sb + refract.b * mat.tb;
    return Color(red, green, blue);
}
Ray Refract(Ray ray, Hit hit,Material mat)
{
    float parameter1 = dot(ray.dir, hit.N);
    float ni, nr;
    if (parameter1 < 0)//go into the material
    {
        nr = mat.ior;
        ni = 1;
    }
    else
    {
        ni = mat.ior;
        nr = 1;
    }
    
    //anglei = acos(dot(hit.N, ray.dir) / ray.dir.magnitude());
    //angler = asin(ni * sin(anglei) / nr);

    Dir3D T = ni / nr * (ray.dir - hit.N * parameter1) - hit.N * sqrt(1 - ni * ni / nr / nr * (1 - parameter1 * parameter1));
    return Ray(hit.interacts, T.normalized());
}
Ray Reflect(Ray ray, Dir3D norm, Point3D start)
{
    Dir3D R = ray.dir - 2 * dot(ray.dir, norm) * norm;
    R = R.normalized();
    return Ray(start, R);
}