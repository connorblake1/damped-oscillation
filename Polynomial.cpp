//
//  Polynomial.cpp
//  RampCurves
//
//  Created by Connor Blake on 9/27/21.
//  Copyright © 2021 Connor Blake. All rights reserved.
//
#include <iostream>
#include <math.h>
#include "Polynomial.hpp"
#include "Complex.hpp"
void Polynomial::init(int deg, float x0in, float y0in, float * incoeffs, int w, int h, bool split) {
    this->degree = deg;
    if (split) {this->setCoeffsSplit(incoeffs);}
    else {
        this->setCoeffs(incoeffs);
        this->kscale=1;
    }
    this->x0 = x0in;
    this->y0 = y0in;
    this->w = w;
    this->h = h;
    //printout();
}
void Polynomial::graph(sf::RenderWindow * win, int y1, int y2, int nin, sf::Color col) {

    sf::Vertex line1[] = { sf::Vertex(sf::Vector2f(.05*w,0),sf::Color::White), sf::Vertex(sf::Vector2f(.05*w,h),sf::Color::White)};
    sf::Vertex line2[] = { sf::Vertex(sf::Vector2f(0,map(0,y1,y2,.95*h, .05*h)),sf::Color::White), sf::Vertex(sf::Vector2f(w,map(0,y1,y2,.95*h, .05*h)),sf::Color::White)};
    (*win).draw(line1, 2, sf::Lines);
    (*win).draw(line2, 2, sf::Lines);
    sf::VertexArray graph(sf::Points, nin);
    sf::CircleShape c1;
    c1.setRadius(10);
    
    c1.setFillColor(sf::Color::Green);
    c1.setPosition(.05*w,map(y0,y1,y2,.95*h, .05*h));
    win->draw(c1);
    sf::CircleShape c2;
    c2.setRadius(10);
    c2.setFillColor(sf::Color::Green);
    c2.setPosition(.95*w,map(0,y1,y2,.95*h, .05*h));
    win->draw(c2);
    for (int n = 0; n < nin; n++) {
        graph[n].color = col;
        float xmap = (float)(n)/nin*x0;
        graph[n].position = sf::Vector2f( map(xmap,0.0,x0,.05*w,.95*w), map(evaluate(xmap), y1, y2, .95*h, .05*h) );}
        (*win).draw(graph);}

sf::Vector2f * Polynomial::mapCoords(float fx, float y1, float y2) {
    sf::Vector2f * out = new sf::Vector2f(map(fx,0.0,x0,.05*w,.95*w), map(evaluate(fx), y1, y2, .95*h, .05*h));
    return out;
}
void Polynomial::init(int deg) {
    this->degree = deg;
    coeffs = new float[degree+1];
}
void Polynomial::setCoeffs(float * incoeffs) {
    //coeffs = new float[degree+1];
    this->coeffs = incoeffs;
}
void Polynomial::setCoeffsSplit(float * incoeffs) {
//    for (int i = 0; i < degree-1; i++) {
//        std::cout<<"incoeff " << i << "  " << incoeffs[i] << std::endl;
//    }
    coeffs = new float[degree+1];
    this->kscale = -1.0*y0/x0/incoeffs[0];
    //std::cout << kscale <<  " " << degree << std::endl;
    for (int i = degree; i >= 0; i--) {
        if (i == degree) {
            coeffs[degree] = kscale;}
        else if (i == 0) {
            coeffs[i] = y0;
        }
        else if (degree-1 == i) {
            coeffs[i] =kscale*(incoeffs[i-1]-x0);
        }
        else {
            coeffs[i] = kscale*(incoeffs[i-1]-x0*incoeffs[i]);}}
}
float Polynomial::evaluate(float x) {
    float r = 0;
    for (int i = 0; i <= degree; i++) {
        r += coeffs[i]*pow(x, i);}
    return r;
}
Polynomial * Polynomial::derivative() {
    Polynomial * p1 = new Polynomial();
    float * newcoeffs = new float[degree];
    for (int i = 0; i < degree; i++) {
        newcoeffs[i] = coeffs[i+1]*(i+1);
    }
    p1->init(degree-1, x0, y0, newcoeffs, w, h, false);
    return p1;
}

float Polynomial::map(float in, float l1, float u1, float l2, float u2) {
    float r = (in-l1)/(u1-l1)*(u2-l2)+l2;
    //map(plotter[vals*p+v]/(bounds[p*2+1]-bounds[p*2]), 0.0,1.0,y+.05*h,y+.95*h)
    return r;}

void Polynomial::printout() {
    for (int i = degree; i >= 0; i--) {
        std::cout << coeffs[i] << "x^" << i << " + " ;}
    std::cout << std::endl;
}
float * Polynomial::getRoots() {
    if (degree == 4) {
        float A = coeffs[0];
        float B = coeffs[1];
        float C = coeffs[2];
        float D = coeffs[3];
        float E = coeffs[4];
        double alpha = (-3.*B*B/8./A/A) + (C/A);
        double beta = (B*B*B/8./A/A/A) - (B*C/2./A/A) + (D/A);
        double gamma = (-3.*B*B*B*B/256./A/A/A/A) + (C*B*B/16./A/A/A) - (B*D/4./A/A) + E/A;
        double P = -alpha*alpha/12. -gamma;
        float Q = -alpha*alpha*alpha/108. + alpha*gamma/3. - beta*beta/8.;
        float det1 =Q*Q/4.+P*P*P/27.;
        Complex R(-Q/2,0);
        if (det1<0) {
            R.i = sqrt(-det1);}
        else {
            R.r += sqrt(det1);}
        R.powz(.3333);
        Complex y(-5./6.*alpha,0);
        Complex PZ(P,0);
        PZ.divz(R.r*3,R.i*3);
        Complex PZ1(R.r-PZ.r,R.i-PZ.i);
        y.addz(&PZ1);
        Complex W(2*y.r+alpha,2*y.i);
        W.powz(.5);
        Complex R1(-B/4./A,0);
        Complex R2(-B/4./A,0);
        Complex R3(-B/4./A,0);
        Complex R4(-B/4./A,0);
        Complex a3a2y(-3*alpha-2*y.r,-2*y.i);
        Complex b2w(2*beta,0);
        b2w.divz(&W);
        W.multz(.5);
        Complex d1(a3a2y.r-b2w.r,a3a2y.i-b2w.i);
        Complex d2(a3a2y.r-b2w.r,a3a2y.i-b2w.i);
        Complex d3(a3a2y.r+b2w.r,a3a2y.i+b2w.i);
        Complex d4(a3a2y.r+b2w.r,a3a2y.i+b2w.i);
        R1.addz(&W);
        R2.addz(&W);
        R3.subz(&W);
        R4.subz(&W);
        d1.powz(.5);
        d1.multz(.5);
        d2.powz(.5);
        d2.multz(.5);
        d3.powz(.5);
        d3.multz(.5);
        d4.powz(.5);
        d4.multz(.5);
        R1.addz(&d1);
        R2.subz(&d2);
        R3.addz(&d3);
        R4.subz(&d4);
        float * roots1 = new float[8];
        roots1[0] = R1.r;
        roots1[1] = R1.i;
        roots1[2] = R2.r;
        roots1[3] = R2.i;
        roots1[4] = R3.r;
        roots1[5] = R3.i;
        roots1[6] = R4.r;
        roots1[7] = R4.i;
        return roots1;
    }
}
