#pragma once
#include<iostream>
#include<list>
#include<vector>
#include<set>
#include<unordered_set>
#include<map>
#include<algorithm>
using namespace std;

class Point //����
{
public:
	double x, y;

	Point operator+(const Point rhs)const;
	Point operator-(const Point rhs)const;
	double operator*(const Point rhs)const;
	double cross(const Point rhs)const;
	friend Point operator*(double t, const Point rhs);

	bool operator<(const Point rhs)const;//�ֵ���
	bool operator==(const Point rhs)const;//���Ƹ��㾫�ȵ�Ӱ��

	friend ostream& operator<<(ostream& out, const Point& p);

	Point();
	Point(double x, double y);
	~Point();
};

class Line //�߶���
{
public:
	Point P, Q;//�߶ζ˵�

	static Point E;//ɨ�����¼��㣬�������<�Ľ�������

	bool operator<(const Line rhs)const;//ɨ������ص���ṹ
	bool operator==(const Line rhs)const;

	friend void findIntersection(list<Line>& lines, map<Point, vector<Line>>& intersections);

	Line();
	Line(Point P, Point Q);
	Line(double Px, double Py, double Qx, double Qy);
	~Line();
};

class Polygon //������Խ��������
{
public:
	struct Vertex //�����������꣬�ڽ�ָ�룬�ཻ��ǩ������
	{
		Point p;
		Vertex* next=0, * last=0;
		bool isMarked;
	};
	Vertex* head = 0;//ָ��˫��������ʼ�ڵ�
	bool orientation = false;//����ε�����ע�������㼯����ʱ�������ε�����һ��

	bool interiorTest(Point c);//�жϵ��Ƿ��ڶ�����ڲ����������
	void append(Point c);
	void refreshSign();

	Polygon();
	Polygon(Point* pg, int n);
	Polygon(const Polygon& obj);
	~Polygon();
};

class Yin //Yin����
{
public:
	list<Polygon> spadjor;
	bool sign = false;

	Yin inverse();
	Yin meet(Yin& rhs);      
	Yin join(Yin& rhs);

	void intersect(Yin& obj);
	bool interiorTest(Point c);
	bool onTest(Point c,Point d);
	void cut(list<list<Point>>& out);
	void getBettiNum(int& b0,int& b1);
	void append(Polygon PL);
	void clearLabel();

	Yin();
	Yin(list<list<Point>>& segment);
	~Yin();
};