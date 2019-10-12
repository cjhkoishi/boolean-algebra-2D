#pragma once
#include<iostream>
#include<fstream>
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
	int ID = 0;

	static Point E;//扫描线事件点，与操作符<的结果相关联

	bool operator<(const Line rhs)const;//扫描线相关的序结构
	bool operator==(const Line rhs)const;//无向意义下的相等比较

	bool isSameLine(const Line rhs)const;
	bool isInside(Point c)const;

	//friend void findIntersection(list<Line>& lines, map<Point, vector<Line>>& intersections);

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
		Vertex* next = 0, * last = 0;
		bool isMarked = false;//是否是分割点
	};
	Vertex* head = 0;//指向双向链表起始节点
	bool orientation = false;//多边形的绕向，被动变量

	bool interiorTest(Point c);//判断点是否在多边形内部。符号相关
	void append(Point c);//在head顶点之前插入一个顶点
	void refreshOri();//计算定向并刷新变量

	void reverse();
	bool split(list<Polygon>& result);

	Polygon();
	Polygon(Point* pg, int n);
	Polygon(string filename);
	Polygon(const Polygon& obj);
	~Polygon();
};

class Yin //Yin集类
{
public:
	list<Polygon> spadjor;
	bool sign = false;
	//布尔运算
	Yin inverse();
	Yin meet(Yin rhs);
	Yin join(Yin rhs);

	void intersect(Yin& obj);//求多边形集合的交点，并插入新交点数据
	bool interiorTest(Point c);//Yin set的内部检测
	bool onTest(Point c, Point d);//有向线段cd是否与多边形重合且同向
	bool onTestInv(Point c, Point d);//有向线段cd是否与多边形重合且逆向
	void cut(list<list<Point>>& out);//利用标记点分割成折线段
	void getBettiNum(int& b0, int& b1);//计算betti数
	void append(Polygon PL);//添加多边形
	void clearLabel();//清空标记点的标记

	void load(string datafiles[], int num);
	void move(Point p);

	Yin();
	Yin(list<list<Point>>& segment);//利用线段集合生成spadjor
	Yin(string datafiles[], int num);
	~Yin();
};