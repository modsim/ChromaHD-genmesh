/*
 * File:   Bead.h
 * Author: guowei
 *
 * Created on September 21, 2012, 4:47 PM
 */

#ifndef BEAD_H
#define	BEAD_H

#include <math.h>

class Bead {
public:
    Bead();
    Bead(const Bead& orig);
    Bead(double x, double y, double z, double r);
    virtual ~Bead();
    double getZ() const;
    double getY() const;
    double getX() const;
    double getR() const;
    void translate(double offsetx, double offsety, double offsetz);
    void scale(double factor);
    void scaleRadius(double factor);
    bool neighbour(Bead* bead, double eps, double &dx, double &dy, double &dz) const;
    void setTag(int tag);
    void setCTag(int ctag);
    int getTag();
    int getCTag();
private:
    double x, y, z, r;
    int tag;
    int ctag;
};

#endif	/* BEAD_H */

