/*
 * File:   Bead.cpp
 * Author: guowei
 *
 * Created on September 21, 2012, 4:47 PM
 */

#include "Bead.h"

#include <iostream>
#include <iomanip>

Bead::Bead() {
}

Bead::Bead(double x, double y, double z, double r) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->r = r;
}


Bead::Bead(const Bead& orig) {
}

Bead::~Bead() {
}

double Bead::getZ() const {
    return z;
}

double Bead::getY() const {
    return y;
}

double Bead::getX() const {
    return x;
}

double Bead::getR() const {
    return r;
}


void Bead::translate(double offsetx, double offsety, double offsetz)
{
    this->x = offsetx + this->x;
    this->y = offsety + this->y;
    this->z = offsetz + this->z;
}

void Bead::scale(double factor)
{
    this->x = factor * this->x;
    this->y = factor * this->y;
    this->z = factor * this->z;
    this->r = factor * this->r;
}


bool Bead::neighbour(Bead * bead, double eps, double &dx, double &dy, double &dz) const{

    double x1 = this->x;
    double x2 = bead->getX();

    double y1 = this->y;
    double y2 = bead->getY();

    double z1 = this->z;
    double z2 = bead->getZ();

    double r1 = this->r;
    double r2 = bead->getR();

    dx = x2-x1;
    dy = y2-y1;
    dz = z2-z1;

    double dist = sqrt( pow(dx, 2) + pow(dy, 2) + pow(dz, 2) );

    return (dist < r1 + r2 + eps);

}

void Bead::setTag(int tag)
{
    this->tag = tag;
}

void Bead::setCTag(int ctag)
{
    this->ctag = ctag;
}

int Bead::getTag()
{
    return this->tag;
}

int Bead::getCTag()
{
    return this->ctag;
}

