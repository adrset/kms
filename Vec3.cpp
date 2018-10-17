#include "Vec3.hpp"

Vec3::Vec3(double xx, double yy, double zz) : x(xx), y(yy), z(zz)
{

}


void Vec3::setR(const Vec3& copy){
	
	this->x = copy.getX();
	this->y = copy.getY();
	this->z = copy.getZ();
	
}

Vec3 operator+(const Vec3& left, const Vec3& right) {
         Vec3 vec;
         vec.setR(left);
         vec.x += right.getX();  
         vec.y += right.getY();  
         vec.z += right.getZ(); 
         return vec;
} 

Vec3 operator*(const Vec3& left, const Vec3& right) {
         Vec3 vec;
         vec.x = right.getZ() * left.getY() - right.getY() * left.getZ(); 
         vec.y = -right.getZ() * left.getX() + right.getX() * left.getZ(); 
         vec.z = right.getY() * left.getX() - right.getX() * left.getY(); 
         return vec;
}

Vec3 operator*(const Vec3& left, float right) {
         Vec3 vec;
         vec.setR(left);
         vec.x *= right; 
         vec.y *= right; 
         vec.z *= right; 
         return vec;
}

Vec3 operator-(const Vec3& left, const Vec3& right){
	
	 Vec3 vec;
     vec.setR(left);
     vec.x -= right.getX();  
     vec.y -= right.getY();  
     vec.z -= right.getZ(); 
     return vec;
	
}

void Vec3::zero(){
	this->x = 0;
	this->y = 0;
	this->z = 0;
}


Vec3 operator/(const Vec3& left, float right) {
         Vec3 vec;
         vec.setR(left);
         vec.x /= right; 
         vec.y /= right; 
         vec.z /= right; 
         return vec;
}

Vec3 operator*(float left, const Vec3& right) {
         Vec3 vec(right);
         vec = vec * left;
         return vec;
}

std::ostream& operator<<(std::ostream& stream, const Vec3& obj){
	stream << obj.x << " " << obj.y << " " << obj.z;
	return stream;
	
}
