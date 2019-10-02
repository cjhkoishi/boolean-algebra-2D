#include "bool_alg.h"

struct PointInfo {
	list<Line> U;
	list<Line> L;
	list<Line> C;
};
void findNextEvent(Line L1, Line L2, Point p, map<Point, PointInfo>& Q);

Point Point::operator+(const Point rhs) const
{
	return Point(x + rhs.x, y + rhs.y);
}

Point Point::operator-(const Point rhs) const
{
	return Point(x - rhs.x, y - rhs.y);
}

double Point::operator*(const Point rhs) const
{
	return x * rhs.x + y * rhs.y;
}

double Point::cross(const Point rhs) const
{
	return x * rhs.y - y * rhs.x;
}

bool Point::operator<(const Point rhs) const
{
	return (y == rhs.y) ? (x > rhs.x) : (y < rhs.y);
}

bool Point::operator==(const Point rhs) const
{
	return (abs(x - rhs.x) < 1e-5) && (abs(y - rhs.y) < 1e-5);
}

Point::Point() :x(0), y(0)
{
}

Point::Point(double x, double y) : x(x), y(y)
{
}

Point::~Point()
{
}

Point operator*(double t, const Point rhs)
{
	return Point(t * rhs.x, t * rhs.y);
}

ostream& operator<<(ostream& out, const Point& p)
{
	return out << "(" << p.x << "," << p.y << ")";
}

void findIntersection(list<Line>& lines, map<Point, vector<Line>>& intersections)
{
	map<Point, PointInfo> Q;
	set<Line> T;
	for (list<Line>::iterator i = lines.begin(); i != lines.end(); i++) {
		Point* upper, * lower;
		if (i->Q < i->P) {
			upper = &i->P;
			lower = &i->Q;
		}
		else {
			upper = &i->Q;
			lower = &i->P;
		}
		Q[*upper].U.push_back(*i);
		Q[*lower].L.push_back(*i);
	}
	while (!Q.empty())
	{
		Point p = Q.rbegin()->first;
		PointInfo pi = Q[p];
		Q.erase(--Q.end());
		//working

		if (pi.C.size() + pi.L.size() + pi.U.size() > 1)
		{
			intersections[p].insert(intersections[p].end(), pi.C.begin(), pi.C.end());
			intersections[p].insert(intersections[p].end(), pi.L.begin(), pi.L.end());
			intersections[p].insert(intersections[p].end(), pi.U.begin(), pi.U.end());
		}
		for (list<Line>::iterator i = pi.L.begin(); i != pi.L.end(); i++)
			T.erase(*i);
		for (list<Line>::iterator i = pi.C.begin(); i != pi.C.end(); i++)
			T.erase(*i);
		Line::E = p;
		for (list<Line>::iterator i = pi.U.begin(); i != pi.U.end(); i++)
			T.insert(*i);
		for (list<Line>::iterator i = pi.C.begin(); i != pi.C.end(); i++)
			T.insert(*i);
		if (pi.U.empty() && pi.C.empty())
		{
			set<Line>::iterator sl = T.begin();
			set<Line>::iterator sr = T.begin();
			T.insert(*pi.L.begin());
			sl = sr = T.find(*pi.L.begin());
			sl--;
			sr++;
			T.erase(*pi.L.begin());

			if (sl != T.end() && sr != T.end())
				findNextEvent(*sl, *sr, p, Q);
		}
		else {
			set<Line> UC;
			UC.insert(pi.U.begin(), pi.U.end());
			UC.insert(pi.C.begin(), pi.C.end());
			set<Line>::iterator i_smin = T.find(*UC.begin());
			set<Line>::iterator i_smax = T.find(*UC.rbegin());
			Line smin = *i_smin;
			Line smax = *i_smax;
			set<Line>::iterator sl = --i_smin;
			set<Line>::iterator sr = ++i_smax;
			if (sl != T.end())
				findNextEvent(*sl, smin, p, Q);
			if (sr != T.end())
				findNextEvent(smax, *sr, p, Q);
		}
		//working

	}
}

void findNextEvent(Line L1, Line L2, Point p, map<Point, PointInfo>& Q) {
	Point v = L1.Q - L1.P;
	Point w = L2.Q - L2.P;
	double cm = v.cross(w);
	double err = 1e-6;
	if (cm != 0)
	{
		double t = (L2.P - L1.P).cross(w) / cm;
		double s = (L2.P - L1.P).cross(v) / cm;
		if (t > -err && t < 1 + err && s > -err && s < 1 + err)
		{
			Point insection = L1.P + t * v;
			if (insection < p || insection == p)
			{
				PointInfo& pi = Q[insection];
				if ((find(pi.L.begin(), pi.L.end(), L1) == pi.L.end()) && (find(pi.U.begin(), pi.U.end(), L1) == pi.U.end()))
					pi.C.push_back(L1);
				if ((find(pi.L.begin(), pi.L.end(), L2) == pi.L.end()) && (find(pi.U.begin(), pi.U.end(), L2) == pi.U.end()))
					pi.C.push_back(L2);
				pi.C.unique();
			}
		}
	}
}

Point Line::E;

bool Line::operator<(const Line rhs) const
{
	bool isParallel1 = P.y == Q.y;
	bool isParallel2 = rhs.P.y == rhs.Q.y;
	double M = isParallel1 ? E.x : P.x + (Q.x - P.x) * (P.y - E.y) / (P.y - Q.y);
	double N = isParallel2 ? E.x : rhs.P.x + (rhs.Q.x - rhs.P.x) * (rhs.P.y - E.y) / (rhs.P.y - rhs.Q.y);
	if (abs(M - N) > 1e-5)
		return M < N;
	else if (isParallel2)
		return !isParallel1;
	else if (isParallel1)
		return false;
	else {
		double k1 = (Q.x - P.x) / (P.y - Q.y);
		double k2 = (rhs.Q.x - rhs.P.x) / (rhs.P.y - rhs.Q.y);
		if (k1 != k2)
			return k1 < k2;
		else
			return P < rhs.P || P == rhs.P && Q < rhs.Q;
	}
}

bool Line::operator==(const Line rhs) const
{
	return P == rhs.P && Q == rhs.Q || P == rhs.Q && Q == rhs.P;
}

Line::Line()
{
}

Line::Line(Point P, Point Q) :P(P), Q(Q)
{
}

Line::Line(double Px, double Py, double Qx, double Qy) : P(Px, Py), Q(Qx, Qy)
{
}

Line::~Line()
{
}

void Polygon::intersect(Polygon& polygon, map<Point, vector<Line>>& intersections)
{
	list<Line>lines;

	Vertex* ii = head, * jj = polygon.head;
	do {
		lines.push_back(Line(ii->p, ii->next->p));
		ii = ii->next;
	} while (ii != head);
	do {
		lines.push_back(Line(jj->p, jj->next->p));
		jj = jj->next;
	} while (jj != polygon.head);

	map<Point, PointInfo> Q;
	set<Line> T;
	for (list<Line>::iterator i = lines.begin(); i != lines.end(); i++) {
		Point* upper, * lower;
		if (i->Q < i->P) {
			upper = &i->P;
			lower = &i->Q;
		}
		else {
			upper = &i->Q;
			lower = &i->P;
		}
		Q[*upper].U.push_back(*i);
		Q[*upper].U.unique();
		Q[*lower].L.push_back(*i);
		Q[*upper].L.unique();
	}
	while (!Q.empty())
	{
		Point p = Q.rbegin()->first;
		PointInfo pi = Q[p];
		Q.erase(--Q.end());
		//working

		if (pi.C.size() > 0 || pi.L.size() + pi.U.size() > 2)
		{
			intersections[p].insert(intersections[p].end(), pi.C.begin(), pi.C.end());
			intersections[p].insert(intersections[p].end(), pi.L.begin(), pi.L.end());
			intersections[p].insert(intersections[p].end(), pi.U.begin(), pi.U.end());
		}
		for (list<Line>::iterator i = pi.L.begin(); i != pi.L.end(); i++)
			T.erase(*i);
		for (list<Line>::iterator i = pi.C.begin(); i != pi.C.end(); i++)
			T.erase(*i);
		Line::E = p;
		for (list<Line>::iterator i = pi.U.begin(); i != pi.U.end(); i++)
			T.insert(*i);
		for (list<Line>::iterator i = pi.C.begin(); i != pi.C.end(); i++)
			T.insert(*i);
		if (pi.U.empty() && pi.C.empty())
		{
			set<Line>::iterator sl = T.begin();
			set<Line>::iterator sr = T.begin();
			T.insert(*pi.L.begin());
			sl = sr = T.find(*pi.L.begin());
			sl--;
			sr++;
			T.erase(*pi.L.begin());

			if (sl != T.end() && sr != T.end())
				findNextEvent(*sl, *sr, p, Q);
		}
		else {
			set<Line> UC;
			UC.insert(pi.U.begin(), pi.U.end());
			UC.insert(pi.C.begin(), pi.C.end());
			set<Line>::iterator i_smin = T.find(*UC.begin());
			set<Line>::iterator i_smax = T.find(*UC.rbegin());
			Line smin = *i_smin;
			Line smax = *i_smax;
			set<Line>::iterator sl = --i_smin;
			set<Line>::iterator sr = ++i_smax;
			if (sl != T.end())
				findNextEvent(*sl, smin, p, Q);
			if (sr != T.end())
				findNextEvent(smax, *sr, p, Q);
		}
		//working
	}


}

bool Polygon::interiorTest(Point c)
{
	Vertex* i = head;
	if (i == 0) {
		return sign;
	}
	int cn = 0;
	do {
		Point v1 = i->p - c;
		Point v2 = i->next->p - c;
		Point w = v2 - v1;
		if (w.y != 0)
		{
			double t = -v1.y / w.y;
			double s = -v1.cross(w) / w.y;
			if (s > 0)
			{
				if (t >= 0 && t < 1 && w.y>0)
					cn++;
				else if (t > 0 && t <= 1 && w.y < 0)
					cn--;
			}
		}
		i = i->next;
	} while (i != head);
	return abs(cn) % 2 == 1;
}

void Polygon::cut(list<list<Vertex*>>& out)
{
}

Polygon::Polygon()
{
}

Polygon::Polygon(Point* pg, int n)
{
	if (n != 0)
	{
		head = new Vertex();
		head->p = pg[0];
		Vertex* j = head;
		for (int i = 1; i < n; i++) {
			j->next = new Vertex();
			j->next->last = j;
			j->next->p = pg[i];
			j = j->next;
		}
		j->next = head;
		head->last = j;
	}
}

Polygon::Polygon(const Polygon& obj)
{
	if (obj.head != 0) {
		head = new Vertex();
		head->p = obj.head->p;
		Vertex* i = obj.head->next;
		Vertex* k = head;
		while (obj.head != i) {
			Vertex* j = new Vertex();
			j->p = i->p;
			j->isMarked = i->isMarked;
			k->next = j;
			j->last = k;
			k = j;
			i = i->next;
		}
		k->next = head;
		head->last = k;
	}
}

Polygon::~Polygon()
{
	Vertex* i = head;
	if (i == 0)
		return;
	do {
		Vertex* j = i->next;
		delete i;
		i = j;
	} while (i != head);
}

Yin Yin::inverse()
{
	return Yin();
}

Yin Yin::meet(Yin& rhs)
{
	list<list<Point>> out1, out2, fin;
	intersect(rhs);
	cut(out1);
	rhs.cut(out2);
	for (list<list<Point>>::iterator i = out1.begin(); i != out1.end(); i++) {
		Point testp = 0.5 * (*i->begin() + *(++i->begin()));
		if (rhs.interiorTest(testp) || rhs.onTest(testp))
			fin.push_back(*i);
	}
	for (list<list<Point>>::iterator i = out2.begin(); i != out2.end(); i++) {
		Point testp = 0.5 * (*i->begin() + *(++i->begin()));
		if (interiorTest(testp) || onTest(testp))
			fin.push_back(*i);
	}
	return Yin(fin);
}

Yin Yin::join(Yin& rhs)
{
	return Yin();
}

void Yin::intersect(Yin& obj)//功能：进行多边形相交算法，获得非连接点交点，并插入原数据集合
{
	struct Comp {
		bool operator() (const Line& lhs, const Line& rhs) const {
			return lhs.P < rhs.P || lhs.P == rhs.P && lhs.Q < rhs.Q;
		}
	};
	map<Line, Polygon::Vertex*, Comp> intersectTable;

	//构造线段链表
	list<Line> lines;
	for (list<Polygon>::iterator i = spadjor.begin(); i != spadjor.end(); i++) {
		Polygon::Vertex* ii = i->head;
		do {
			Line l(ii->p, ii->next->p);
			lines.push_back(l);
			intersectTable[l] = ii;
			ii = ii->next;
		} while (ii != i->head);
	}
	for (list<Polygon>::iterator i = obj.spadjor.begin(); i != obj.spadjor.end(); i++) {
		Polygon::Vertex* ii = i->head;
		do {
			Line l(ii->p, ii->next->p);
			lines.push_back(l);
			intersectTable[l] = ii;
			ii = ii->next;
		} while (ii != i->head);
	}

	//计算交点

	map<Point, PointInfo> intersections;//交点=>所属线段 映射
	map<Point, PointInfo> Q;
	set<Line> T;
	for (list<Line>::iterator i = lines.begin(); i != lines.end(); i++) {
		Point* upper, * lower;
		if (i->Q < i->P) {
			upper = &i->P;
			lower = &i->Q;
		}
		else {
			upper = &i->Q;
			lower = &i->P;
		}
		Q[*upper].U.push_back(*i);
		Q[*upper].U.unique();
		Q[*lower].L.push_back(*i);
		Q[*upper].L.unique();
	}
	while (!Q.empty())
	{
		Point p = Q.rbegin()->first;
		PointInfo pi = Q[p];
		Q.erase(--Q.end());
		//working

		if (pi.C.size() > 0 || pi.L.size() + pi.U.size() > 2)
		{
			intersections[p] = pi;
		}
		for (list<Line>::iterator i = pi.L.begin(); i != pi.L.end(); i++)
			T.erase(*i);
		for (list<Line>::iterator i = pi.C.begin(); i != pi.C.end(); i++)
			T.erase(*i);
		Line::E = p;
		for (list<Line>::iterator i = pi.U.begin(); i != pi.U.end(); i++)
			T.insert(*i);
		for (list<Line>::iterator i = pi.C.begin(); i != pi.C.end(); i++)
			T.insert(*i);
		if (pi.U.empty() && pi.C.empty())
		{
			set<Line>::iterator sl = T.begin();
			set<Line>::iterator sr = T.begin();
			T.insert(*pi.L.begin());
			sl = sr = T.find(*pi.L.begin());
			sl--;
			sr++;
			T.erase(*pi.L.begin());

			if (sl != T.end() && sr != T.end())
				findNextEvent(*sl, *sr, p, Q);
		}
		else {
			set<Line> UC;
			UC.insert(pi.U.begin(), pi.U.end());
			UC.insert(pi.C.begin(), pi.C.end());
			set<Line>::iterator i_smin = T.find(*UC.begin());
			set<Line>::iterator i_smax = T.find(*UC.rbegin());
			Line smin = *i_smin;
			Line smax = *i_smax;
			set<Line>::iterator sl = --i_smin;
			set<Line>::iterator sr = ++i_smax;
			if (sl != T.end())
				findNextEvent(*sl, smin, p, Q);
			if (sr != T.end())
				findNextEvent(smax, *sr, p, Q);
		}
		//working
	}

	//插入&标记
	map<Polygon::Vertex*, list<Point>> buffer;
	for (map<Point, PointInfo>::iterator i = intersections.begin(); i != intersections.end(); i++) {
		Point p = i->first;
		PointInfo pi = i->second;
		for (list<Line>::iterator j = pi.U.begin(); j != pi.U.end(); j++) {
			if (j->P == p)
				intersectTable[*j]->isMarked = true;
			if (j->Q == p)
				intersectTable[*j]->next->isMarked = true;
		}
		for (list<Line>::iterator j = pi.L.begin(); j != pi.L.end(); j++) {
			if (j->P == p)
				intersectTable[*j]->isMarked = true;
			if (j->Q == p)
				intersectTable[*j]->next->isMarked = true;
		}

		for (list<Line>::iterator j = pi.C.begin(); j != pi.C.end(); j++) {
			buffer[intersectTable[*j]].push_back(p);
		}
	}

	for (map<Polygon::Vertex*, list<Point>>::iterator i = buffer.begin(); i != buffer.end(); i++) {
		list<Point>& np = i->second;
		np.sort([&](Point& p1, Point& p2)->bool {Point s1 = p1 - i->first->p, s2 = p2 - i->first->p; return abs(s1.x) + abs(s1.y) < abs(s2.x) + abs(s2.y); });
		Polygon::Vertex* start = i->first;
		Polygon::Vertex* end = i->first->next;
		for (list<Point>::iterator j = np.begin(); j != np.end(); j++) {
			Polygon::Vertex* nv = new Polygon::Vertex();
			nv->p = *j;
			nv->next = end;
			nv->last = start;
			nv->isMarked = true;
			start->next = nv;
			end->last = nv;
			start = nv;
		}
	}
}

bool Yin::interiorTest(Point c)
{
	bool res = false;
	for (list<Polygon>::iterator i = spadjor.begin(); i != spadjor.end(); i++) {
		bool b = i->interiorTest(c);
		res ^= b;
	}
	return res;
}

bool Yin::onTest(Point c)
{
	for (list<Polygon>::iterator i = spadjor.begin(); i != spadjor.end(); i++) {
		Polygon::Vertex* j = i->head;
		do {
			if (abs((j->p - c).cross(j->next->p - c)) < 1e-12)
				return true;
		} while (j != i->head);
	}
	return false;
}

void Yin::cut(list<list<Point>>& out)
{
	out.clear();
	for (list<Polygon>::iterator i = spadjor.begin(); i != spadjor.end(); i++) {
		Polygon::Vertex* start = i->head;
		if (start == 0)
			continue;
		do {
			if (start->isMarked)
				break;
			if (start->next == i->head)
			{
				start->isMarked = true;
				break;
			}
			start = start->next;
		} while (true);

		Polygon::Vertex* it = start->next;
		out.push_front(list<Point>());
		list<list<Point>>::iterator linker = out.begin();
		linker->push_back(start->p);
		do {
			linker->push_back(it->p);
			if (it->isMarked) {
				out.push_front(list<Point>());
				linker = out.begin();
				linker->push_back(it->p);
			}
			it = it->next;
		} while (it != start);
		linker->push_back(it->p);
	}
}

Yin::Yin()
{
}

Yin::Yin(list<list<Point>>& segments)
{
	struct SegmentInfo {
		list<list<Point>*> O;
		list<list<Point>*> I;
	};
	map<Point, SegmentInfo> si;
	for (list<list<Point>>::iterator i = segments.begin(); i != segments.end(); i++) {
		si[*i->begin()].O.push_back(&*i);
		si[*i->rbegin()].I.push_back(&*i);
	}

	while (!si.empty()) {
		list<Point>* beta = *si.begin()->second.O.begin();//当前折线段
		list<Point>::iterator s = beta->begin();//顶点遍历
		Point ls = *s, m = *s;//折线段起点 初始点
		spadjor.push_front(Polygon());
		list<Polygon>::iterator J = spadjor.begin();//被构造的目标多边形
		Polygon::Vertex** v = &J->head;//被修改的向前指针
		do {
			*v = new Polygon::Vertex();
			(*v)->p = *s;
			v = &(*v)->next;
			if (si.find(*s) != si.end())//若是分界点
			{
				si[*s].I.remove(beta);
				si[ls].O.remove(beta);
				if (si[*s].O.size() + si[*s].I.size() == 0)
					si.erase(*s);
				if (si[ls].O.size() + si[ls].I.size() == 0)
					si.erase(ls);
				SegmentInfo currentsi = si[*s];
				double minarg;//working
				for (list<list<Point>*>::iterator j = currentsi.O.begin(); j != currentsi.O.end(); j++) {
					beta = *j;
				}
				s = beta->begin();
				ls = *s;
			}
			s++;
		} while (!(*s == m));
		si[*s].I.remove(beta);
		si[ls].O.remove(beta);
		if (si[*s].O.size() + si[*s].I.size() == 0)
			si.erase(*s);
		if (si[ls].O.size() + si[ls].I.size() == 0)
			si.erase(ls);
		*v = J->head;
		Polygon::Vertex* fix = J->head;
		do {
			fix->next->last = fix;
			fix = fix->next;
		} while (fix != J->head);
	}
}

Yin::~Yin()
{
}
