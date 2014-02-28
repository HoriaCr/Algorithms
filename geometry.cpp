#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES 1
#include <math.h>
#include <vector>

using namespace std;

template<class T>
class point { 
private :
	T x , y;
public:
	point() {
		this->x = 0;
		this->y = 0;
	}
	point(const T &a,const T &b) {
		this->x = a;
		this->y = b;
	}
	int X() {
		return this->x;
	}
	int Y() {
		return this->y;
	}
	void operator += (point<T> &P) {
		this->x+=P.X();
		this->y+=P.Y();
	}
	void operator -= (point<T> &P){
		this->x-=P.X();
		this->y-=P.Y();
	}
	T sqrDistanceTo(point<T> &P);
  friend ostream &operator << (ostream &stream, point ob) {
	  stream<<ob.X()<<' '<<ob.Y()<<'\n';
	  return stream;
  }
  template<class T>
  friend istream &operator >> (istream &stream, point<T> &ob) {
	  T a , b;
	  stream>>a>>b;
	  ob = point<T>(a,b);
	  return stream;
  }
};

template<class T> T point<T>::sqrDistanceTo(point<T> &P) {
	return (this->x - P.X())*(this->x - P.X()) + (this->y - P.Y())*(this->y - P.Y());
}

template<class T> 
class circle { 
private:
	T ray;
	point<T> center;
public:
	circle() {
		this->ray = 1;
		this->center = point<int>(0,0);
	}
	circle(const T &r,const point<T> &c) {
		this->ray = r;
		this->center = c;
	}
	point<T> C() {
		return this->center;
	}
	T R() {
		return this->ray;
	}
	T D() {
		return 2*this->ray;
	}
	double area() {
		return ray*ray*M_PI;
	}
	double circumference() {
		return 2.0*M_PI*ray;
	}
	bool overlapsCircle(circle<T> &C) {
		return this->center.sqrDistanceTo(C.C()) < (this->ray*this->ray + C.R()*C.R());
	}
	bool touchesCircle(circle<T> &C) {
		return this->center.sqrDistanceTo(C.C()) ==  (this->ray*this->ray + C.R()*C.R());
	}
  template<class T>
  friend istream &operator >> (istream &stream, circle<T> &ob) {
	  T a , b , r;
	  stream>>a>>b>>r;
	  ob = circle<T>(r,point<T>(a,b));
	  return stream;
  }
};
