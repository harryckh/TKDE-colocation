#ifndef smallestenclosingcircle_h
#define smallestenclosingcircle_h

#include <vector>
#include "data_struct.h"
#include "costenum.h"

struct Point final {
	
	public: double x;
	public: double y;
	
	
	public: Point subtract(const Point &p) const;
	
	public: double distance(const Point &p) const;
	
	// Signed area / determinant thing
	public: double cross(const Point &p) const;
	
};


struct Circle final {
	
	public: static const Circle INVALID;
	
	private: static const double MULTIPLICATIVE_EPSILON;
	
	
	public: Point c;   // Center
	public: double r;  // Radius
	
	
	public: bool contains(const Point &p) const;
	
	public: bool contains(const std::vector<Point> &ps) const;
	
};


/* 
 * (Main function) Returns the smallest circle that encloses all the given points.
 * Runs in expected O(n) time, randomized. Note: If 0 points are given, a circle of
 * negative radius is returned. If 1 point is given, a circle of radius 0 is returned.
 */
Circle makeSmallestEnclosingCircle(const std::vector<Point> &points);

// For unit tests
Circle makeDiameter(const Point &a, const Point &b);
Circle makeCircumcircle(const Point &a, const Point &b, const Point &c);

disk_t* MCC(obj_set_t* obj_set_v);
#endif 

