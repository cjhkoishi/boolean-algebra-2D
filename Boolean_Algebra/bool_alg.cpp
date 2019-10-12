#define PI 3.14159265358
#include "bool_alg.h"


struct PointInfo {
	list<Line> U;
	list<Line> L;
	list<Line> C;
};
void findNextEvent(Line L1, Line L2, Point p, map<Point, PointInfo>& Q);
double mod(double a, double b);

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
	return (abs(x - rhs.x) < 1e-8) && (abs(y - rhs.y) < 1e-8);
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

/*void findIntersection(list<Line>& lines, map<Point, vector<Line>>& intersections)
{
	map<Point, PointInfo> Q;
	set<Line> T;
	for (auto i = lines.begin(); i != lines.end(); i++) {
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
		for (auto i = pi.L.begin(); i != pi.L.end(); i++)
			T.erase(*i);
		for (auto i = pi.C.begin(); i != pi.C.end(); i++)
			T.erase(*i);
		Line::E = p;
		for (auto i = pi.U.begin(); i != pi.U.end(); i++)
			T.insert(*i);
		for (auto i = pi.C.begin(); i != pi.C.end(); i++)
			T.insert(*i);
		if (pi.U.empty() && pi.C.empty())
		{
			auto sl = T.begin();
			auto sr = T.begin();
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
			auto i_smin = T.find(*UC.begin());
			auto i_smax = T.find(*UC.rbegin());
			Line smin = *i_smin;
			Line smax = *i_smax;
			auto sl = --i_smin;
			auto sr = ++i_smax;
			if (sl != T.end())
				findNextEvent(*sl, smin, p, Q);
			if (sr != T.end())
				findNextEvent(smax, *sr, p, Q);
		}
		//working

	}
}*/

void findNextEvent(Line L1, Line L2, Point p, map<Point, PointInfo>& Q) {
	Point v = L1.Q - L1.P;
	Point w = L2.Q - L2.P;
	double cm = v.cross(w);
	if (abs(cm) > 1e-12)
	{
		double t = (L2.P - L1.P).cross(w) / cm;
		Point intersection = L1.P + t * v;
		bool isInsideLine1 = (intersection < L1.P ^ intersection < L1.Q);
		bool isInsideLine2 = (intersection < L2.P ^ intersection < L2.Q);
		bool isEndPoint1 = (intersection == L1.P || intersection == L1.Q);
		bool isEndPoint2 = (intersection == L2.P || intersection == L2.Q);
		if ((isInsideLine1 || isEndPoint1) && (isInsideLine2 || isEndPoint2))
		{
			if (intersection < p || intersection == p)
			{
				PointInfo& pi = Q[intersection];
				bool CF1 = find(pi.C.begin(), pi.C.end(), L1) == pi.C.end();
				bool CF2 = find(pi.C.begin(), pi.C.end(), L2) == pi.C.end();
				if (CF1 && !isEndPoint1)
					pi.C.push_back(L1);
				if (CF2 && !isEndPoint2)
					pi.C.push_back(L2);
			}
		}
	}
}

double mod(double a, double b)
{
	while (!(a < b / 2 && a >= -b / 2)) {
		a += a >= b / 2 ? -b : b;
	}
	return a;
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
	else {
		Point v1 = Q - P;
		Point v2 = rhs.Q - rhs.P;
		v1 = v1 < Point(0, 0) ? v1: Point(0, 0) - v1;
		v2 = v2 < Point(0, 0) ? v2: Point(0, 0) - v2;
		double angleless = v1.cross(v2);
		if (angleless != 0)
			return (E.x > M + 2e-5) ^ (angleless > 0);
		else if (P == rhs.P && Q == rhs.Q)
			return ID < rhs.ID;
		else
			return P < rhs.P || P == rhs.P && Q < rhs.Q;
	}
	/*else if (isParallel2)
		return !isParallel1;
	else if (isParallel1)
		return false;
	else {
		double k1 = (Q.x - P.x) / (P.y - Q.y);
		double k2 = (rhs.Q.x - rhs.P.x) / (rhs.P.y - rhs.Q.y);
		if (k1 != k2)
			return (E.x > M + 2e-5) ^ (k1 < k2);
		else if (P == rhs.P && Q == rhs.Q)
			return ID < rhs.ID;
		else
			return P < rhs.P || P == rhs.P && Q < rhs.Q;
	}*/
}

bool Line::operator==(const Line rhs) const
{
	return P == rhs.P && Q == rhs.Q;
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


bool Polygon::interiorTest(Point c)
{
	Vertex* i = head;
	if (i == 0) {
		return orientation;
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

void Polygon::append(Point c)
{
	if (head == 0) {
		head = new Vertex();
		head->p = c;
		head->next = head;
		head->last = head;
	}
	else {
		Vertex* v = new Vertex();
		v->p = c;
		v->next = head;
		v->last = head->last;
		head->last->next = v;
		head->last = v;
	}
}

void Polygon::refreshOri()
{
	if (head == 0)
		orientation = true;
	else {
		Vertex* i = head;
		Vertex* sidePoint = head;
		do {
			if (i->p < sidePoint->p)
				sidePoint = i;
			i = i->next;
		} while (i != head);
		Point v1 = sidePoint->p - sidePoint->last->p;
		Point v2 = sidePoint->next->p - sidePoint->p;
		orientation = v1.cross(v2) > 0;
	}
}

void Polygon::reverse()
{
	Vertex* i = head;
	if (head == 0)
		return;
	do {
		Vertex* t = i->next;
		i->next = i->last;
		i->last = t;
		i = i->next;
	} while (i != head);
}

bool Polygon::split(list<Polygon>& result)
{
	int num = 0;
	map<Point, list<Point>::iterator> pool;
	list<Point> chain;
	chain.push_back(head->p);
	pool[head->p] = --chain.end();
	auto i = head;
	do {
		bool t = pool.find(i->next->p) != pool.end();
		if (t)
		{
			result.push_front(Polygon());
			auto newPl = result.begin();
			auto splitpos = pool[i->next->p];
			for (auto j = splitpos; j != chain.end(); j++) {
				if (j != splitpos)
					pool.erase(*j);
				newPl->append(*j);

			}
			newPl->refreshOri();
			num++;
			chain.erase(++splitpos, chain.end());
		}
		else
		{
			chain.push_back(i->next->p);
			pool[i->next->p] = --chain.end();
		}
		i = i->next;
	} while (i != head);
	return num > 1;
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

Polygon::Polygon(string filename)
{
	fstream f(filename, ios::in);
	double input;
	vector<double> buffer;
	while (!f.eof()) {
		f >> input;
		buffer.push_back(input);
	}
	f.close();
	int n = buffer.size() / 2;
	for (int i = 0; i < n - 1; i++) {
		//Point p((buffer[i] * 100 - 84.88986354429699-400)*20+400, (buffer[n + i] * 100 - 23.22797462977070-300)*20+300);
		Point p(buffer[i] * 100 + 50, buffer[n + i] * 100 + 50);
		append(p);
	}
}

Polygon::Polygon(const Polygon& obj)
{
	if (obj.head != 0) {
		orientation = obj.orientation;
		head = new Vertex();
		head->p = obj.head->p;
		head->isMarked = obj.head->isMarked;
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
	Yin null, lhs = *this;
	lhs.intersect(null);
	list<list<Point>> out;
	lhs.cut(out);
	for (auto i = out.begin(); i != out.end(); i++) {
		list<Point> buf = *i;
		i->clear();
		while (!buf.empty())
		{
			i->push_back(*buf.rbegin());
			buf.pop_back();
		}
	}
	Yin res(out);
	res.sign = lhs.sign ^ true;
	//clearLabel();
	return res;
}

Yin Yin::meet(Yin rhs)
{
	Yin lhs = *this;
	list<list<Point>> out1, out2, fin;
	lhs.intersect(rhs);
	lhs.cut(out1);
	rhs.cut(out2);
	for (auto i = out1.begin(); i != out1.end(); i++) {
		Point c = *i->begin();
		Point d = *(++i->begin());
		Point testp = 0.5 * (c + d);
		bool isInterior = rhs.interiorTest(testp);
		bool isCoincide = rhs.onTest(c, d);
		if (isInterior || isCoincide)
			fin.push_back(*i);
	}
	for (auto i = out2.begin(); i != out2.end(); i++) {
		Point c = *i->begin();
		Point d = *(++i->begin());
		Point testp = 0.5 * (c + d);
		bool isInterior = lhs.interiorTest(testp);
		bool isCoincide = lhs.onTest(c, d);
		bool isInv = lhs.onTestInv(c, d);
		if (isInterior && !isCoincide && !isInv)
			fin.push_back(*i);
	}

	Yin res(fin);
	res.sign = lhs.sign && rhs.sign;
	//clearLabel();
	//rhs.clearLabel();
	return res;
}

Yin Yin::join(Yin rhs)
{
	Yin rl = inverse();
	Yin rr = rhs.inverse();
	Yin rres = rl.meet(rr);
	return rres.inverse();
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
	int id = 1;
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		Polygon::Vertex* ii = i->head;
		do {
			Line l(ii->p, ii->next->p);
			l.ID = id;
			lines.push_back(l);
			intersectTable[l] = ii;
			ii = ii->next;
		} while (ii != i->head);
		id++;
	}
	for (auto i = obj.spadjor.begin(); i != obj.spadjor.end(); i++) {
		Polygon::Vertex* ii = i->head;
		do {
			Line l(ii->p, ii->next->p);
			l.ID = id;
			lines.push_back(l);
			intersectTable[l] = ii;
			ii = ii->next;
		} while (ii != i->head);
		id++;
	}

	//计算交点

	map<Point, PointInfo> intersections;//交点=>所属线段 映射
	map<Point, PointInfo> Q;
	set<Line> T;
	for (auto i = lines.begin(); i != lines.end(); i++) {
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

		//working




		for (auto i = pi.L.begin(); i != pi.L.end(); i++)
			T.erase(*i);
		for (auto i = pi.C.begin(); i != pi.C.end(); i++)
			T.erase(*i);
		Line::E = p;
		for (auto i = pi.U.begin(); i != pi.U.end(); i++)
			T.insert(*i);
		for (auto i = pi.C.begin(); i != pi.C.end(); i++)
			T.insert(*i);
		if (pi.U.empty() && pi.C.empty())
		{
			auto sl = T.begin();
			auto sr = T.begin();
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
			auto i_smin = T.find(*UC.begin());
			auto i_smax = T.find(*UC.rbegin());
			Line smin = *i_smin;
			Line smax = *i_smax;
			auto sl = --i_smin;
			auto sr = ++i_smax;
			if (sl != T.end())
				findNextEvent(*sl, smin, p, Q);
			if (sr != T.end())
				findNextEvent(smax, *sr, p, Q);
		}
		//working
		if (pi.C.size() > 0 || pi.L.size() + pi.U.size() > 2)
		{
			intersections[p] = Q[p];
		}
		Q.erase(--Q.end());
	}

	//插入&标记
	map<Polygon::Vertex*, list<Point>> buffer;
	for (auto i = intersections.begin(); i != intersections.end(); i++) {
		Point p = i->first;
		PointInfo pi = i->second;
		for (auto j = pi.U.begin(); j != pi.U.end(); j++) {
			if (j->P == p)
				intersectTable[*j]->isMarked = true;
			if (j->Q == p)
				intersectTable[*j]->next->isMarked = true;
		}
		for (auto j = pi.L.begin(); j != pi.L.end(); j++) {
			if (j->P == p)
				intersectTable[*j]->isMarked = true;
			if (j->Q == p)
				intersectTable[*j]->next->isMarked = true;
		}

		for (auto j = pi.C.begin(); j != pi.C.end(); j++) {
			buffer[intersectTable[*j]].push_back(p);
		}
	}

	for (auto i = buffer.begin(); i != buffer.end(); i++) {
		list<Point>& np = i->second;
		np.sort([&](Point& p1, Point& p2)->bool {Point s1 = p1 - i->first->p, s2 = p2 - i->first->p; return abs(s1.x) + abs(s1.y) < abs(s2.x) + abs(s2.y); });
		Polygon::Vertex* start = i->first;
		Polygon::Vertex* end = i->first->next;
		for (auto j = np.begin(); j != np.end(); j++) {
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
	if (spadjor.empty())
		return sign;
	bool res = false;
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		bool b = i->interiorTest(c);
		res ^= b;
	}
	return res ^ sign;
}

bool Yin::onTest(Point c, Point d)
{
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		Polygon::Vertex* j = i->head;
		do {
			if (c == j->p && abs((d - c).cross(j->next->p - c)) < 1e-12 && (d - c) * (j->next->p - c) > 0)
				return true;
			j = j->next;
		} while (j != i->head);
	}
	return false;
}

bool Yin::onTestInv(Point c, Point d)
{
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		Polygon::Vertex* j = i->head;
		do {
			if (c == j->p && abs((d - c).cross(j->next->p - c)) < 1e-12 && (d - c) * (j->next->p - c) < 0)
				return true;
			j = j->next;
		} while (j != i->head);
	}
	return false;
}

void Yin::cut(list<list<Point>>& out)
{
	out.clear();
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
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
		auto linker = out.begin();
		linker->push_back(start->p);
		do {
			linker->push_back(it->p);
			if (it->isMarked) {
				out.push_front(list<Point>());
				linker = out.begin();
				linker->push_back(it->p);
				//?
				if (it->p == Point(380.83679955096261, 152.66266035223680))
					int sss = 1;
				//?
			}
			it = it->next;
		} while (it != start);
		linker->push_back(it->p);
	}
}

void Yin::getBettiNum(int& b0, int& b1)
{
	b0 = sign;
	b1 = 0;
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		if (i->orientation)
			b0++;
		else
			b1++;
	}
}

void Yin::append(Polygon PL)
{
	PL.refreshOri();
	bool ori = PL.orientation;
	if (!interiorTest(PL.head->p) ^ ori)
		PL.reverse();
	PL.refreshOri();
	spadjor.push_back(PL);
}

void Yin::clearLabel()
{
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		Polygon::Vertex* j = i->head;
		do {
			j->isMarked = false;
			j = j->next;
		} while (j != i->head);
	}
}

void Yin::load(string datafiles[], int num)
{
	for (int i = 0; i < num; i++) {
		Polygon PL(datafiles[i]);
		append(PL);
	}
}

void Yin::move(Point p)
{
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		auto j = i->head;
		do {
			j->p = j->p + p;
			j = j->next;
		} while (j != i->head);
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
	for (auto i = segments.begin(); i != segments.end(); i++) {
		si[*i->begin()].O.push_back(&*i);
		si[*i->rbegin()].I.push_back(&*i);
	}

	while (true) {
		auto t = si.begin();
		while (t != si.end()) {
			if (!t->second.O.empty())
				break;
			t++;
		}
		if (t == si.end())
			return;

		list<Point>* beta = *t->second.O.begin();
		Polygon Jor;
		Point start;
		auto current = beta->begin();
		start = *current;
		bool flag = true;

		do {
			Jor.append(*current);
			current++;
			if (current == beta->end()) {
				si[*beta->begin()].O.remove(beta);
				si[*beta->rbegin()].I.remove(beta);
				list<list<Point>*>& sobo = si[*beta->rbegin()].O;
				if (sobo.empty())
				{
					flag = false;
					break;
				}
				//working
				Point Ivec = *beta->rbegin() - *(++beta->rbegin());
				beta = *sobo.begin();
				Point Ovec = *(++beta->begin()) - *beta->begin();
				double Iarg = atan2(Ivec.y, Ivec.x);
				double Oarg = atan2(Ovec.y, Ovec.x);
				double maxArg = mod(Oarg - Iarg, 2 * PI);
				for (auto it = ++sobo.begin(); it != sobo.end(); it++)
				{
					Ovec = *(++(*it)->begin()) - *(*it)->begin();
					Oarg = atan2(Ovec.y, Ovec.x);
					double arg = mod(Oarg - Iarg, 2 * PI);
					if (arg > maxArg) {
						maxArg = arg;
						beta = *it;
					}
				}
				current = ++(beta->begin());
			}

		} while (!(*current == start));
		if (!flag) {

			continue;
		}
		si[*beta->begin()].O.remove(beta);
		si[*beta->rbegin()].I.remove(beta);
		list<Polygon> atom;
		Jor.split(atom);
		spadjor.insert(spadjor.end(), atom.begin(),atom.end());
	}

}

Yin::Yin(string datafiles[], int num)
{
	load(datafiles, num);
}

Yin::~Yin()
{
}
