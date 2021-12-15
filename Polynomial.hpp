//
//  Polynomial.hpp
//  RampCurves
//
//  Created by Connor Blake on 9/27/21.
//  Copyright © 2021 Connor Blake. All rights reserved.
//

#ifndef Polynomial_hpp
#define Polynomial_hpp
#include <SFML/Graphics.hpp>
#include <iostream>
#include <stdio.h>
class Polynomial {
private:
    float * coeffs;
    float kscale, x0, y0;
    static constexpr float g = 5.0;
    int degree, w, h;
    
public:
    Polynomial() {}
    void printout(); //done
    void init(int deg, float x0in, float y0in, float * incoeffs, int w, int h, bool split); //done
    void init(int deg);
    void graph(sf::RenderWindow * win, int y1, int y0, int nin, sf::Color col);
    void setCoeffs(float * coeffs); //done
    void setCoeffsSplit(float * coeffs);
    float evaluate(float x); //done
    float map(float in, float l1, float u1, float l2, float u2); //done
    Polynomial * derivative(); //done
    sf::Vector2f * mapCoords(float fx, float y1, float y2);
    int getDeg() {return this->degree;}
    float getCoeff(int index) {return coeffs[index];}
    float * getCoeffs() {return this->coeffs;}
    float * getRoots();
};

#endif /* Polynomial_hpp */
