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
	int ID = 0;

	static Point E;//ɨ�����¼��㣬�������<�Ľ�������

	bool operator<(const Line rhs)const;//ɨ������ص���ṹ
	bool operator==(const Line rhs)const;//���������µ���ȱȽ�

	bool isSameLine(const Line rhs)const;
	bool isInside(Point c)const;

	//friend void findIntersection(list<Line>& lines, map<Point, vector<Line>>& intersections);

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
		Vertex* next = 0, * last = 0;
		bool isMarked = false;//�Ƿ��Ƿָ��
	};
	Vertex* head = 0;//ָ��˫��������ʼ�ڵ�
	bool orientation = false;//����ε����򣬱�������

	bool interiorTest(Point c);//�жϵ��Ƿ��ڶ�����ڲ����������
	void append(Point c);//��head����֮ǰ����һ������
	void refreshOri();//���㶨��ˢ�±���

	void reverse();
	bool split(list<Polygon>& result);

	Polygon();
	Polygon(Point* pg, int n);
	Polygon(string filename);
	Polygon(const Polygon& obj);
	~Polygon();
};

class Yin //Yin����
{
public:
	list<Polygon> spadjor;
	bool sign = false;
	//��������
	Yin inverse();
	Yin meet(Yin rhs);
	Yin join(Yin rhs);

	void intersect(Yin& obj);//�����μ��ϵĽ��㣬�������½�������
	bool interiorTest(Point c);//Yin set���ڲ����
	bool onTest(Point c, Point d);//�����߶�cd�Ƿ��������غ���ͬ��
	bool onTestInv(Point c, Point d);//�����߶�cd�Ƿ��������غ�������
	void cut(list<list<Point>>& out);//���ñ�ǵ�ָ�����߶�
	void getBettiNum(int& b0, int& b1);//����betti��
	void append(Polygon PL);//��Ӷ����
	void clearLabel();//��ձ�ǵ�ı��

	void load(string datafiles[], int num);
	void move(Point p);

	Yin();
	Yin(list<list<Point>>& segment);//�����߶μ�������spadjor
	Yin(string datafiles[], int num);
	~Yin();
};