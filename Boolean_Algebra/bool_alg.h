#pragma once
#include<iostream>
#include<list>
#include<vector>
#include<set>
#include<unordered_set>
#include<map>
#include<algorithm>
using namespace std;

class Point //点类
{
public:
	double x, y;

	Point operator+(const Point rhs)const;
	Point operator-(const Point rhs)const;
	double operator*(const Point rhs)const;
	double cross(const Point rhs)const;
	friend Point operator*(double t, const Point rhs);

	bool operator<(const Point rhs)const;//字典序
	bool operator==(const Point rhs)const;//控制浮点精度的影响

	friend ostream& operator<<(ostream& out, const Point& p);

	Point();
	Point(double x, double y);
	~Point();
};

class Line //线段类
{
public:
	Point P, Q;//线段端点

	static Point E;//扫描线事件点，与操作符<的结果相关联

	bool operator<(const Line rhs)const;//扫描线相关的序结构
	bool operator==(const Line rhs)const;

	friend void findIntersection(list<Line>& lines, map<Point, vector<Line>>& intersections);

	Line();
	Line(Point P, Point Q);
	Line(double Px, double Py, double Qx, double Qy);
	~Line();
};

class Polygon //有向非自交多边形类
{
public:
	struct Vertex //包含定点坐标，邻接指针，相交标签等属性
	{
		Point p;
		Vertex* next=0, * last=0;
		bool isMarked;
	};
	Vertex* head = 0;//指向双向链表起始节点
	bool orientation = false;//多边形的绕向。注：当顶点集不空时，与多边形的绕向一致

	bool interiorTest(Point c);//判断点是否在多边形内部。符号相关
	void append(Point c);
	void refreshSign();

	Polygon();
	Polygon(Point* pg, int n);
	Polygon(const Polygon& obj);
	~Polygon();
};

class Yin //Yin集类
{
public:
	list<Polygon> spadjor;
	bool sign = false;

	Yin inverse();
	Yin meet(Yin& rhs);      
	Yin join(Yin& rhs);

	void intersect(Yin& obj);
	bool interiorTest(Point c);
	bool onTest(Point c);
	void cut(list<list<Point>>& out);
	void getBettiNum(int& b0,int& b1);
	void append(Polygon PL);
	void clearLabel();

	Yin();
	Yin(list<list<Point>>& segment);
	~Yin();
};