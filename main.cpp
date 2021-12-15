
//
// Disclaimer:
// ----------
//
// This code will work only if you selected window, graphics and audio.
//
// Note that the "Run Script" build phase will copy the required frameworks
// or dylibs to your application bundle so you can execute it on any OS X
// computer.
//
// Your resource files (images, sounds, fonts, ...) are also copied to your
// application bundle. To get the path to these resources, use the helper
// function `resourcePath()` from ResourcePath.hpp
//

#include <SFML/Graphics.hpp>
#include "Polynomial.hpp"
#include "SliderSFML.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
// Here is a small helper for you! Have a look.
#include "ResourcePath.hpp"
int STATE = 0; //0 is home screen, 1 is animation
float m1 = 1; //kg
float m2 = 1;
float k1 = 200;
float k2 = 200;
float b1 = 1;
float b2 = 1;
float A, B, C, D, E;
float * roots = NULL;
std::vector<double> x(8);
std::vector<double> gauss(std::vector< std::vector<double> > A) {
    int n = A.size();
    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    std::vector<double> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}
int binomialCoeff(int n, int k)
{
    int res = 1;
    if (k > n - k)
    k = n - k;
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }
     
    return res;
}

double nDeriv(int ni, double a, double b, bool isCos) {
    double c = 0;
    for (int i = 0; i < ni+1; i++) {
        if (isCos && i % 2 == 1) {continue;}
        else if (!isCos && i % 2 == 0) {continue;}
        int c1 = 0;
        if (isCos) {
            if (((i+1)/2) % 2 == 0) {c1 = 1;}
            else {c1 = -1;}}
        else {
            if ((i/2) % 2 == 0) {c1 = 1;}
            else {c1 = -1;}}
        c += binomialCoeff(ni,i)*pow(a,ni-i)*c1*pow(b,i);
    }
    return c;
}
void bashRoots() {
    int n = 8;
    std::vector<double> line(n+1,0);
    std::vector< std::vector<double> > A(n,line);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            bool isC;
            if (j%2 ==0) {isC=true;}
            else{isC=false;}
            A[i][j] = nDeriv(i, roots[(j/2)*2], roots[(j/2)*2+1], isC);
        }
    }
    float s[8] = {1,0,0,0,0,0,0,0};
    for (int i=0; i<n; i++) {
        A[i][n] = s[i];}
    x = gauss(A);
    for (int i=0; i<n; i++) {
        std::cout << x[i] << " ";
    }
}

int main(int, char const**)
{
    sf::RenderWindow window(sf::VideoMode(800, 600), "SFML window");
    sf::Font fontMain;
    fontMain.loadFromFile("/Users/connorblake/Documents/SORTOld/RampCurves/RampCurves/arial.ttf");
   
    int fontsize = 24;
    int sFont = 20;
    float lalign = 50;
    float ualign = 100;
    float vspacer = 65;
    float optionOffset = 300;
    sf::Text poly1("Solution: ", fontMain,sFont);
    poly1.setFillColor(sf::Color::Yellow);
    poly1.setPosition(lalign,ualign-60);
    
    sf::Text m1t("M1",fontMain,fontsize);
    m1t.setFillColor(sf::Color::White);
    m1t.setPosition(lalign,ualign+1*vspacer);
    SliderSFML sliderM1(lalign+optionOffset,ualign+1*vspacer);
    sliderM1.create(1, 10,false);
    
    sf::Text m2t("M2",fontMain,fontsize);
    m2t.setFillColor(sf::Color::White);
    m2t.setPosition(lalign,ualign+2*vspacer);
    SliderSFML sliderM2(lalign+optionOffset,ualign+2*vspacer);
    sliderM2.create(1, 10,false);
    
    sf::Text k1t("K1",fontMain,fontsize);
    k1t.setFillColor(sf::Color::White);
    k1t.setPosition(lalign,ualign+3*vspacer);
    SliderSFML sliderK1(lalign+optionOffset,ualign+3*vspacer);
    sliderK1.create(1, 10,false);
    
    sf::Text k2t("K2",fontMain,fontsize);
    k2t.setFillColor(sf::Color::White);
    k2t.setPosition(lalign,ualign+4*vspacer);
    SliderSFML sliderK2(lalign+optionOffset,ualign+4*vspacer);
    sliderK2.create(1, 10,false);
    
    sf::Text b1t("B1",fontMain,fontsize);
    b1t.setFillColor(sf::Color::White);
    b1t.setPosition(lalign,ualign+5*vspacer);
    SliderSFML sliderB1(lalign+optionOffset,ualign+5*vspacer);
    sliderB1.create(1, 10,false);
    
    sf::Text b2t("B2",fontMain,fontsize);
    b2t.setFillColor(sf::Color::White);
    b2t.setPosition(lalign,ualign+6*vspacer);
    SliderSFML sliderB2(lalign+optionOffset,ualign+6*vspacer);
    sliderB2.create(1, 10,false);
    
    
    
    
    
    Polynomial * p = new Polynomial();
    A = m1/k2;
    B = b1/k2 + b2*m1/m2/k2;
    C = 1 + k1/k2 + b1*b2/m2/k2 + m1/m2;
    D = b2/m2*(1+k1/k2) + b1/m2;
    E = k1/m2;
    //float * coeffs = new float[5] {A,B,C,D,E};
    //float * coeffs = new float[5] {1,2,-41,-42,360};
    float coeffs [5] = {1,2,-41,-42,360};

    p->init(4, 20, 10, coeffs, 1000, 1000, false);
    roots = p->getRoots();
    std::cout << "Roots" << std::endl;
    for (int i = 0; i < 8; i ++) {
        std::cout << roots[i] << std::endl;}
    
    while (window.isOpen())
    {
        
        sf::Event event;
        while (window.pollEvent(event))
        {
            // Close window: exit
            if (event.type == sf::Event::Closed) {
                window.close();
            }

            // Escape pressed: exit
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape) {
                window.close();
            }
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::R) {
                m1 = sliderM1.getSliderValue();
                m2 = sliderM2.getSliderValue();
                k1 = sliderK1.getSliderValue();
                k2 = sliderK2.getSliderValue();
                b1 = sliderB1.getSliderValue();
                b2 = sliderB2.getSliderValue();
                A = m1/k2;
                B = b1/k2 + b2*m1/m2/k2;
                C = 1 + k1/k2 + b1*b2/m2/k2 + m1/m2;
                D = b2/m2*(1+k1/k2) + b1/m2;
                E = k1/m2;
                float coeffs [5] = {A,B,C,D,E};
                p->setCoeffs(coeffs);
                roots = p->getRoots();
                bashRoots();
            }
        }

        window.clear();
        
        std::stringstream stream3;
        for (int i = 0; i < 4 ;i++) {
            stream3 << std::fixed << std::setprecision(2) << x[2*i];
            //stream3 << "c" << 2*i << "e^(";
            stream3 << "e^(";
            stream3 << std::fixed << std::setprecision(2) << roots[2*i];
            stream3 << "t)cos(";
            stream3 << std::fixed << std::setprecision(2) << roots[2*i+1];
            stream3 << "t) + ";
            stream3 << std::fixed << std::setprecision(2) << x[2*i+1];
            stream3 << "e^(";
            stream3 << std::fixed << std::setprecision(2) << roots[2*i];
            stream3 << "t)sin(";
            stream3 << std::fixed << std::setprecision(2) << roots[2*i+1];
            stream3 << "t)";
            if (i != 3) {
                stream3 << " + \n";
            }
        }
        poly1.setString(stream3.str());
        window.draw(poly1);
        window.draw(m1t);
        window.draw(m2t);
        window.draw(k1t);
        window.draw(k2t);
        window.draw(b1t);
        window.draw(b2t);
        sliderM1.draw(window);
        sliderM2.draw(window);
        sliderK1.draw(window);
        sliderK2.draw(window);
        sliderB1.draw(window);
        sliderB2.draw(window);
        window.display();
    }

    return EXIT_SUCCESS;
}
