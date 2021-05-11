//
//  point.hh
//  hw3_code
//
//  Created by Elaine Cunha on 5/7/21.
//

#ifndef point_hh
#define point_hh

// Data structure of (x,y,z) positions with theta and the location of the point in the main array.
struct point {
    /** The x, y, and z coordinates. */
    double x;
    double y;
    double z;
    /** theta. */
    double theta;
    /** Index of point in main array. */
    int idx;
    point() {}
    point(double x_,double y_,double z_, double theta_, int idx_) : x(x_), y(y_), z(z_), theta(theta_), idx(idx_) {}
    point(const point &o) : x(o.x), y(o.y), z(o.z), theta(o.theta), idx(o.idx) {}
    void set(double x_,double y_,double z_, double theta_, int idx_) {
        x=x_;y=y_;z=z_;theta=theta_;idx=idx_;
    }
    void clear() {
        x=999999.;y=999999.;z=999999.;theta=999999.;idx=999999.;
    }
//    const bool operator<(const point &o) {
//        return z<o.z;
//    }
//    const bool operator==(const point &o) {
//        return z==o.z;
//    }
};

#endif /* point_hh */
