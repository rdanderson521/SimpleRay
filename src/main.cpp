
#include <vector>
#include <cmath>
#include <iostream>

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

  Vec3X<T> operator*(T b);

  Vec3X<T> operator/(T b);
};

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
};


class Object
{
public:
  virtual bool intersect(Ray ray, RayHit& hit)=0;

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
    return true;
  }
  return false;
}


class Camera
{
public:
  //https://gabrielgambetta.com/computer-graphics-from-scratch/09-perspective-projection.html

  Camera(double d, double fov, Vec2i size, Vec3d o = {0,0,0}, Vec3d dir = {0,0,1});

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
const double FOV = 1.0472;//60 deg
const double NEAR = 0.5;

int main (int argc, char *argv[])
{
  CImg<unsigned char> img(WIDTH,HEIGHT,1,3,0);

  Camera cam(NEAR,FOV,C_SIZE,{0,0,0},{0,0,1});

  Sphere obj({0,0,5},2,{255,0,0});
  Sphere obj2({0,0,3},1,{255,0,0});

  for (int i = 0; i < HEIGHT; i++)
  {
    for (int j = 0; j < WIDTH; j++)
    {
      Ray ray;
      cam.getRay({j,i},ray);

      RayHit hit;

      if(obj.intersect(ray, hit))
      {
        img.atXY(j,i) = 255;
      }
      if(obj2.intersect(ray, hit))
      {
        img.atXY(j,i,0,1) = 255;
      }
    }
  }

  CImgDisplay main_disp(img,"disp");
  while (1)
  main_disp.wait();

}
