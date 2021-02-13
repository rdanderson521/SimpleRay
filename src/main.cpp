
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "CImg.h"
using namespace cimg_library;

template <class T>
struct Vec3X
{
  T x,y,z;

  T dot(Vec3X<T> b);

  T abs();

  Vec3X<T> norm();

  Vec3X<T> operator+(Vec3X<T> b);

  Vec3X<T> operator-(Vec3X<T> b);

  Vec3X<T> operator-();

  Vec3X<T> operator*(T b);

  Vec3X<T> operator/(T b);

  T operator[](size_t idx);
};

template <class T>
Vec3X<T> operator*(T a, Vec3X<T> b)
{
  Vec3X<T> c = {a*b.x, a*b.y, a*b.z};
  return c;
}

template <class T>
T Vec3X<T>::dot(Vec3X<T> b)
{
  return this->x*b.x + this->y*b.y + this->z*b.z;
}

template <class T>
T Vec3X<T>::abs()
{
  return static_cast<T>(sqrt(x*x + y*y + z*z));
}

template <class T>
Vec3X<T> Vec3X<T>::norm()
{
  T abs = this->abs();
  return *this/abs;
}

template <class T>
Vec3X<T> Vec3X<T>::operator+(Vec3X<T> b)
{
  Vec3X<T> c = {this->x+b.x, this->y+b.y, this->z+b.z};
  return c;
}

template <class T>
Vec3X<T> Vec3X<T>::operator-(Vec3X<T> b)
{
  Vec3X<T> c = {this->x-b.x, this->y-b.y, this->z-b.z};
  return c;
}

template <class T>
Vec3X<T> Vec3X<T>::operator-()
{
  Vec3X<T> c = {-this->x, -this->y, -this->z};
  return c;
}

template <class T>
Vec3X<T> Vec3X<T>::operator*(T b)
{
  Vec3X<T> c = {this->x*b, this->y*b, this->z*b};
  return c;
}

template <class T>
Vec3X<T> Vec3X<T>::operator/(T b)
{
  Vec3X<T> c = {this->x/b, this->y/b, this->z/b};
  return c;
}

template <class T>
T Vec3X<T>::operator[](size_t idx)
{
  switch (idx)
  {
    case (0):
      return x;
    case (1):
      return y;
    case (2):
      return z;
  }
  return 0;
}

typedef Vec3X<double> Vec3d;
typedef Vec3X<int> Vec3i;


template <class T>
struct Vec2X
{
  T x,y;
};

typedef Vec2X<double> Vec2d;
typedef Vec2X<int> Vec2i;


struct BoundingBox
{
  Vec3d origin, size;
};


struct Ray
{
  Vec3d origin, direction;
  int depth;
};


struct RayHit
{
  Vec3i colour;
  Vec3d hitPos;
  Vec3d hitNorm;
  double dist;

  RayHit();
};

RayHit::RayHit()
{
  colour = {0,0,0};
  hitPos = {0,0,0};
  hitNorm = {0,0,0};
  dist = 0;
}


class Object
{
public:
  virtual bool intersect(Ray ray, RayHit& hit) {return false;};

private:
    BoundingBox bounds;
};


class Sphere: public Object
{
public:
  Sphere();

  Sphere(Vec3d o, double r, Vec3i c);

  virtual bool intersect(Ray ray, RayHit& hit);

private:
  Vec3d origin;
  double radius;

  Vec3i colour;
};

Sphere::Sphere()
{
  this->origin = {0,0,0};
  this->radius = 1;
  this->colour = {255,0,0}; // red
}

Sphere::Sphere(Vec3d o, double r, Vec3i c)
{
  this->origin = o;
  this->radius = r;
  this->colour = c;
}

bool Sphere::intersect(Ray ray, RayHit& hit)
{
  // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

  double delta = pow((ray.direction.norm().dot(ray.origin - this->origin)),2) - (pow((ray.origin - this->origin).abs(),2) - pow(this->radius,2) );

  if (delta < 0)
    return false;
  else if (delta >= 0)
  {
    double dist = -(ray.direction.norm().dot(ray.origin - this->origin)) - sqrt(delta);
    hit.hitPos = ray.origin + (ray.direction.norm()*dist);
    hit.dist = dist;
    hit.colour = this->colour;
    hit.hitNorm = (hit.hitPos - this->origin).norm();
    return true;
  }
  return false;
}


class Camera
{
public:
  //https://gabrielgambetta.com/computer-graphics-from-scratch/09-perspective-projection.html

  Camera(double d, double fov, Vec2i size, Vec3d o = {0,0,0}, Vec3d dir = {0,0,1});

  Vec3d origin() {return o;}

  Vec2i project(Vec3d &p);

  void getRay(Vec2i pos, Ray &ray);

private:
  Vec3d o, dir;
  double d; // depth of near plane
  Vec2d sizeV; // viewport size (camera coords)
  Vec2i sizeC; // canvas size (pixel coords)
};

Camera::Camera(double d, double fov, Vec2i size, Vec3d o, Vec3d dir)
{
  this->d = d;
  this->sizeC = size;
  this->o = o;
  this->dir = dir;

  double xMax = tan(fov/2.)/d;

  double fovY = (fov/size.x) * size.y;
  double yMax = tan(fovY/2.)/d;

  this->sizeV.x = (xMax*size.x)/(size.x-1 - (size.x/2));
  this->sizeV.y = (yMax*size.y)/(size.y-1 - (size.y/2));

}

Vec2i Camera::project(Vec3d &p)
{
  Vec2i projected;
  projected.x = (((p.x*d)/p.z)*sizeC.x)/sizeV.x + (sizeC.x/2);
  projected.y = (((p.y*d)/p.z)*sizeC.y)/sizeV.y + (sizeC.y/2);
  return projected;
}

void Camera::getRay(Vec2i pos, Ray &ray)
{
  Vec3d dir;
  dir.x = ((pos.x - (sizeC.x/2)) * sizeV.x) / sizeC.x;
  dir.y = ((pos.y - (sizeC.y/2)) * sizeV.y) / sizeC.y;
  dir.z = d;

  ray.direction = dir.norm();

  ray.origin = o;
}

const int WIDTH = 800;
const int HEIGHT = 600;
const Vec2i C_SIZE = {WIDTH,HEIGHT};
const double FOV = 0.698132;//40 deg
const double NEAR = 0.5;

int main (int argc, char *argv[])
{
  CImg<unsigned char> img(WIDTH,HEIGHT,1,3,0);

  Camera cam(NEAR,FOV,C_SIZE,{0,-2,-5},{0,0,1});

  Object *objArr[2];
  objArr[0] = new Sphere({-1,0,5},3,{255,0,0});
  objArr[1] = new Sphere({4,0,2},1,{0,255,0});

  Vec3d light = {8,-2, -3};

  for (int i = 0; i < HEIGHT; i++)
  {
    for (int j = 0; j < WIDTH; j++)
    {
      Ray ray;
      cam.getRay({j,i},ray);

      RayHit closestHit;
      closestHit.dist = DBL_MAX;

      bool objectHit = false;

      for (Object *obj: objArr)
      {
        RayHit hit;
        if(obj->intersect(ray, hit))
        {
          objectHit = true;
          if (hit.dist < closestHit.dist)
          {
            closestHit = hit;
          }
        }
      }

      Vec3d vecToGround;
      if (ray.direction.y > 0)
      {
        double yFactor = 0.75 - cam.origin().y;
        double rayDirY = ray.direction.y;
        vecToGround =  (ray.direction/ rayDirY) * yFactor;// + cam.origin();

        double dist = vecToGround.abs();

        if (dist < closestHit.dist)
        {
          closestHit.dist = dist;
          closestHit.hitNorm = Vec3d{0,-1,0}.norm();
          closestHit.hitPos = vecToGround + cam.origin();

          Vec3d comp = vecToGround;

          if (comp.x < 0)
          {
            comp.x = -comp.x + 1;
          }
          if (comp.z < 0)
          {
            comp.z = -comp.z + 1;
          }
          if ((int)(comp.x) % 2 == 0)
          {
            closestHit.colour =  ((int)(comp.z) % 2 == 0) ? Vec3i{0,0,255} : Vec3i{255,255,255};
          }
          else
          {
            closestHit.colour =  ((int)(comp.z) % 2 != 0) ? Vec3i{0,0,255} : Vec3i{255,255,255};
          }

        }
      }

      Ray shadowRay;
      shadowRay.depth = 1;
      shadowRay.origin = closestHit.hitPos;
      shadowRay.direction = (light - closestHit.hitPos).norm();

      bool shadowRayHit = false;

      for (Object *obj: objArr)
      {
        RayHit hit;
        if(obj->intersect(shadowRay, hit))
        {
          if (hit.dist > 0.)
            shadowRayHit = true;
        }
      }

      Vec3d lightDir = (light - closestHit.hitPos).norm();

      double diffuse = std::max(closestHit.hitNorm.dot(lightDir),0.);

      Vec3d reflect = (lightDir - 2. * closestHit.hitNorm.norm().dot(lightDir) * closestHit.hitNorm.norm()).norm();

      double specular = pow( std::max(reflect.dot(ray.direction.norm()), 0.),8);

      for (int c = 0; c < 3; c++)
      {
        if (!shadowRayHit)
        {
          img.atXY(j,i,0,c) = std::min( (std::min(diffuse  + 0.1, 1.0) * closestHit.colour[c]) + ( specular * 255), 255.);
        }
        else
        {
          img.atXY(j,i,0,c) = std::min(0.1 * closestHit.colour[c] , 255.);
        }

      }

    }
  }

  CImgDisplay main_disp(img,"disp");
  while (1)
  main_disp.wait();

}
