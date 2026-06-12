/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <math.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    auto partial = [](const float &x) {
        if(x < 0.5)
        {
            return -1 + sqrt(2*x);
        }
        else
        {
            return 1 - sqrt(2-2*x);
        }
    };
    return {partial(sample[0]), partial(sample[1])};
}

float Warp::squareToTentPdf(const Point2f &p) {

    return (1 - abs(p[0])) * (1 - abs(p[1]));
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    auto partialR = [](const float &r)
    {
        return sqrt(r);
    };
    float r = partialR(sample[0]);
    float theta = sample[1]*2*M_PI;
    float x = cos(theta)*r;
    float y = sin(theta)*r;
    return {x, y};
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    return ( p[0]*p[0]+p[1]*p[1] < 1 ) ? 1.0/(M_PI) : 0;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float wz = 2*sample[0]-1;
    float r = sqrt(1-wz*wz);
    float theta = 2*M_PI*sample[1];
    float wx = r*cos(theta);
    float wy = r*sin(theta);
    return {wx, wy, wz};
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {

    return 1.0/(4*M_PI);
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float wz = sample[0];
    float r = sqrt(1-wz*wz);
    float theta = 2*M_PI*sample[1];
    float wx = r*cos(theta);
    float wy = r*sin(theta);
    return {wx, wy, wz};
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {

    return v[2] > 0 ? 1.0/(2*M_PI) : 0;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    Point2f disc = squareToUniformDisk(sample);
    float z = sqrt(1-(disc[0]*disc[0] + disc[1] * disc[1]));
    return {disc[0], disc[1], z};
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {

    return v[2] > 0 ? v[2]/M_PI : 0;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {

    auto invCFD = [&alpha](const float &r)
    {
        return atan(sqrt(-log(1-r)*alpha*alpha));
    };
    float phi = sample[1]*2*M_PI;
    float theta = invCFD(sample[0]);
    float x = sin(theta)*cos(phi);
    float y = sin(theta)*sin(phi);
    float z = cos(theta);
    return {x, y, z};
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {

    float cosTheta = Frame::cosTheta(m);
    float tanTheta = Frame::sinTheta(m)/cosTheta;

    if( m[2] <= 0) return 0;

    return pow(M_E, - tanTheta*tanTheta / (alpha*alpha) ) /
        (alpha*alpha*pow(cosTheta,3)*M_PI);
}

NORI_NAMESPACE_END
