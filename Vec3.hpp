#ifndef _VEC3_HPP
#define _VEC3_HPP
#include <iostream>
#include <cmath>

class Vec3{
	public:
	
	Vec3(double x = 0.0, double y = 0.0, double z = 0.0);
	void setR(const Vec3& copy);
	
	
	double getX() const { return x; }
	double getY() const { return y; }
	double getZ() const{ return z; }
	double getLength() const 
	{
		return sqrt(x*x + y*y + z*z);
	}
	double getLengthSquared() const {
		return (x*x + y*y + z*z);
	}
	
	friend Vec3 operator+(const Vec3& left, const Vec3& right);
	friend Vec3 operator-(const Vec3& left, const Vec3& right);
	friend Vec3 operator*(const Vec3& left, float right);
    friend Vec3 operator*(const Vec3& left, const Vec3& right);
    friend Vec3 operator*(float left, const Vec3& right);
    friend Vec3 operator/(const Vec3& left, float right);
    friend std::ostream& operator<<(std::ostream& stream, const Vec3& obj);
	private:
		double x;
		double y;
		double z;	
};

#endif
