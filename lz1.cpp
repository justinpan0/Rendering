//
//  lz1.cpp
//  smallpt
//
//  Created by Pan, Zimeng on 11/11/16.
//
//

#include "lz1.hpp"

Struct Vec{
    double x,y,z;
    Vec (double x_=0; double y_=0, double z_=0){
        x=x_;
        y=y_;
        z=z_;
    }
    Vec operator +(const Vec &t) const{
        return Vec(x + t.x, y + z.y, z + t.z);
    }
    Vec operator -(const Vec &t) const{
        return Vec(x - t.x, y - t.y, z - t.z);
    }
    Vec operator *(double b) const{
        return Vec(x * b, y * b, z * b);
    }
    Vec mult(const Vec&b) const{
        return Vec(x * b.x, y * b.y, z * b.z);
    }
    Vec &norm(){
        return *this = *this * (1/sqrt(pow(x,2)+pow(y,2)+pow(z,2)));
    }
    double dot(const Vec &b) const{
        return Vec(x * b.x + y * b.y + z * b.z)
    }
    Vec operator%(Vec &b){
        return Vec(y * b.z - z*b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
    
    struct Ray{
        Vec o, d;
        Ray(Vec o_, d_) : o(o_), d(d_) {}
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
            double t, eps=1e-4, b=opdot(r.d), det=b*b - op.dot(op) + pow(rad, 2);
            
            if(det < 0){
                return 0;
            } else {
                return det = sqrt(det);
            }
            
            return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
        }
    }
}
