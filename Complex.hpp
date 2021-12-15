//
//  Complex.hpp
//  DampedHarmonics
//
//  Created by Connor Blake on 12/14/21.
//  Copyright © 2021 Connor Blake. All rights reserved.
//

#ifndef Complex_hpp
#define Complex_hpp
#include <math.h>
#include <stdio.h>
class Complex {
public:
    float r, i;
    Complex(float rin, float iin) {r=rin; i=iin;}
    Complex() {}
    float magz() {return sqrt(r*r+i*i);}
    void addz(Complex * other) {r+=other->r; i+=other->i;}
    void addz(float rin, float iin) {Complex c(rin,iin); this->addz(&c);}
    void subz(Complex * other) {r-=other->r; i-=other->i;}
    void subz(float rin, float iin) {Complex c(rin,iin); this->subz(&c);}
    void multz(float s) {r*=s;i*=s;}
    void powz(float p) {
        float rad = this->magz();
        float theta = atan2(i,r);
        rad = pow(rad,p);
        theta *= p;
        r = rad*cos(theta);
        i = rad*sin(theta);}
    void multz(Complex * other) {
        float io = i;
        float ro = r;
        r = ro*other->r - io*other->i;
        i = ro*other->i + io*other->r;}
    void multz(float rin, float iin) {
        Complex c(rin,iin);
        this->multz(&c);}
    void divz(Complex * other) {
        float io = i;
        float ro = r;
        float iot = other->i;
        float rot = other->r;
        r = (ro*rot+io*iot)/(rot*rot+iot*iot);
        i = (io*rot-ro*iot)/(rot*rot*iot*iot);}
    void divz(float rin, float iin) {
        Complex c(rin,iin);
        this->divz(&c);}
};
#endif /* Complex_hpp */
