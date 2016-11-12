//
//  lz1.cpp
//  smallpt
//
//  Created by Pan, Zimeng on 11/11/16.
//
//

#include <stdlib.h>
#include <cmath>
#include <iostream>

struct Vec {
    double x, y, z;

    Vec(double x_ = 0, double y_ = 0, double z_ = 0) {
        x = x_;
        y = y_;
        z = z_;
    }

    Vec operator+(const Vec &t) const {
        return Vec(x + t.x, y + t.y, z + t.z);
    }

    Vec operator-(const Vec &t) const {
        return Vec(x - t.x, y - t.y, z - t.z);
    }

    Vec operator*(double b) const {
        return Vec(x * b, y * b, z * b);
    }

    Vec mult(const Vec &b) const {
        return Vec(x * b.x, y * b.y, z * b.z);
    }

    Vec &norm() {
        return *this = *this * (1 / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    }

    double dot(const Vec &b) const {
        return x * b.x + y * b.y + z * b.z;
    }

    Vec operator%(Vec &b) {
        return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
};

struct Ray{
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};
enum Refl_t{
    DIFF, SPEC, REFR
};

struct Sphere{
    double rad;
    Vec p, e, c;
    Refl_t refl;

    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_): //position, emission, color
            rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}

    double intersect(const Ray &r) const{
        Vec op = p - r.o;
        double t, eps = 1e-4, b=op.dot(r.d), det=b*b - op.dot(op) + pow(rad, 2);
        if(det < 0) return 0; else det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

Sphere spheres[];

inline double clamp(double x){
    return x > 1 ? 1 : x < 0 ? 0 : x;
}

inline int toInt(double x){
    return int(pow(clamp(x), 1/2.2) * 255 + 0.5);
}

inline bool intersect(const Ray &r, double &t, int &id){
    double n = sizeof(spheres)/ sizeof(Sphere), d, inf = t = 1e20;
    for(int i = int(n); i--;){
        if((d=spheres[i].intersect(r)) && d < t){
            t = d;
            id = i;
        }
        return t < inf;
    }
}

Vec radiance(const Ray &r, int depth, unsigned short *Xi){
    double dist;
    int id = 0;
    if(!intersect(r, dist, id)){
        return Vec();
    }
    const Sphere &obj = spheres[id];
    Vec x = r.o + r.d*dist, n = (x - obj.p).norm(), nl = n.dot(r.d)<0 ? n : n*-1, f = obj.c;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max reflectance
    if (++depth>5)
        if(erand48(Xi) < p) f=f*(1/p); else return obj.e; //R.R.
    if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
        double r1=2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > 0.1 ? Vec(0,1) : Vec(1))%w).norm(), v = w%u;
        Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
        return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));
    } else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
        return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));
    Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
    bool into = n.dot(nl)>0;                // Ray from outside going in?
    double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
    if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
        return obj.e + f.mult(radiance(reflRay,depth,Xi));
    Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
    double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
    double Re=R0+(1-R0) * pow(c, 5), Tr = 1 - Re, P = 0.25 + 0.5 * Re, RP = Re/P, TP = Tr/(1 - P);
    return obj.e + f.mult(depth>2 ? (erand48(Xi) < P ?   // Russian roulette
                                     radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
                          radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}

