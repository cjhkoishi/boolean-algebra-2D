#define PI 3.14159265358
#include "bool_alg.h"

double err = 1e-10;

struct PointInfo {
	list<Line> U;
	list<Line> L;
	list<Line> C;
};
void findIntersection(list<Line>& lines, map<Point, PointInfo>& intersections);
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

double Point::norm()
{
	return sqrt(x * x + y * y);
}

bool Point::operator<(const Point rhs) const
{
	return (abs(y - rhs.y) < err ? (x > rhs.x) : (y < rhs.y)) && !(*this == rhs);
}

bool Point::operator==(const Point rhs) const
{
	return (abs(x - rhs.x) < err) && (abs(y - rhs.y) < err);
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

void findIntersection(list<Line>& lines, map<Point, PointInfo>& intersections)
{
	intersections.clear();

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
		PointInfo& pi = Q[p];

		//working
		for (auto i = T.begin(); i != T.end(); i++) {
			if (i->isInside(p) && find(pi.C.begin(), pi.C.end(), *i) == pi.C.end())
				pi.C.push_back(*i);
		}

		if (pi.C.size() + pi.L.size() + pi.U.size() > 1)
		{
			intersections[p] = pi;
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
		Q.erase(--Q.end());
	}

}

void findNextEvent(Line L1, Line L2, Point p, map<Point, PointInfo>& Q) {
	Point v = L1.Q - L1.P;
	Point w = L2.Q - L2.P;
	double cm = v.cross(w);
	if (abs(cm) != 0)
	{
		double t = (L2.P - L1.P).cross(w) / cm;
		Point intersection = L1.P + t * v;
		bool isInsideLine1 = (intersection < L1.P ^ intersection < L1.Q);
		bool isInsideLine2 = (intersection < L2.P ^ intersection < L2.Q);
		bool isEndPoint1 = (intersection == L1.P || intersection == L1.Q);
		bool isEndPoint2 = (intersection == L2.P || intersection == L2.Q);
		if ((isInsideLine1 || isEndPoint1) && (isInsideLine2 || isEndPoint2))
		{
			if (intersection < p)
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
	/*else if ((L1.P - L2.P).cross(v) == 0) {
		bool isIn[4];
		Point endpoint[4] = { L1.P,L1.Q,L2.P,L2.Q };
		isIn[0] = (L1.P < L2.P ^ L1.P < L2.Q) && !(L1.P == L2.P || L1.P == L2.Q);
		isIn[1] = (L1.Q < L2.P ^ L1.Q < L2.Q) && !(L1.Q == L2.P || L1.Q == L2.Q);
		isIn[2] = (L2.P < L1.P ^ L2.P < L1.Q) && !(L2.P == L1.P || L2.P == L1.Q);
		isIn[3] = (L2.Q < L1.P ^ L2.Q < L1.Q) && !(L2.Q == L1.P || L2.Q == L1.Q);
		for (int i = 0; i < 4; i++)
			if (isIn[i] && (endpoint[i] < p || endpoint[i] == p))
			{
				int tak = i / 2;
				Line tkl = tak == 0 ? L2 : L1;
				PointInfo& pi = Q[endpoint[i]];
				bool CF = find(pi.C.begin(), pi.C.end(), tkl) == pi.C.end();
				if (CF)
					pi.C.push_back(tkl);
			}
	}*/
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
	Point v = P - Q;
	Point w = rhs.P - rhs.Q;
	double cm = v.cross(w);
	if (cm != 0) {
		double t = (rhs.P - P).cross(w) / cm;
		Point inte = t * v + P;
		bool vd = v < Point(0, 0);
		bool wd = w < Point(0, 0);
		bool less = (vd ^ wd) ^ (cm > 0);
		if (inte == E || E < inte)
			return less;
		else
			return !less;
	}
	else {
		v = v < Point(0, 0) ? v : Point(0, 0) - v;
		double dist = v.cross(rhs.P - P);
		if (dist != 0)
			return dist > 0;
	}
	//对于重合线段进行简单区分
	if (P == rhs.P && Q == rhs.Q)
		return ID < rhs.ID;
	else
		return P < rhs.P || P == rhs.P && Q < rhs.Q;
}

bool Line::operator==(const Line rhs) const
{
	return P == rhs.P && Q == rhs.Q;
}

bool Line::isInside(Point c)const
{
	Point v1 = P - c;
	Point v2 = Q - c;
	Point w = Q - P;
	double dist = abs(v1.cross(w) / w.norm());
	bool closed = dist < err && v1 * w <= 0 && v2 * w >= 0;
	return closed && !(c == P) && !(c == Q);
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

int Polygon::postion(Point c)
{
	Vertex* i = head;
	if (i == 0) {
		return orientation;
	}
	int cn = 0;
	bool flag = false;
	do {
		Line l(i->p, i->next->p);
		Point v1 = l.P - c;
		Point v2 = l.Q - c;
		Point w = v2 - v1;
		if (v1.y * v2.y > 0) {
			i = i->next;
			continue;
		}

		if (l.isInside(c) || c == l.P || c == l.Q) {
			flag = true;
			break;
		}
		if (w.y != 0)
		{
			bool dir = w.y > 0;
			double s = -v1.cross(w);
			double t = -v1.y;
			bool uptest = dir && s > 0 && t >= 0 && t < w.y;
			bool downtest = !dir && s < 0 && t < 0 && t >= w.y;
			if (uptest || downtest)
				cn++;
		}

		i = i->next;
	} while (i != head);
	return flag ? 2 : abs(cn) % 2;
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

bool Polygon::check()
{
	list<Line> lines;

	Polygon::Vertex* i = head;
	do {
		Line l(i->p, i->next->p);
		lines.push_back(l);
		i = i->next;
	} while (i != head);

	map<Point, PointInfo> intersections;
	findIntersection(lines,intersections);

	for (auto i = intersections.begin(); i != intersections.end(); i++) {
		if (i->second.C.size() > 0 || i->second.U.size() + i->second.L.size() > 2) {
			return false;
		}
	}

	return true;
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
	int n = (int)buffer.size() / 2;
	for (int i = 0; i < n - 1; i++) {
		Point p(buffer[i], buffer[n + i]);
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
		int r = rhs.postion(c, d);
		if (r == 1 || r == 2)
			fin.push_back(*i);
	}
	for (auto i = out2.begin(); i != out2.end(); i++) {
		Point c = *i->begin();
		Point d = *(++i->begin());
		int r = lhs.postion(c, d);
		if (r == 1)
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
	/*map<Point, PointInfo> Q;
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
		PointInfo& pi = Q[p];

		//working
		for (auto i = T.begin(); i != T.end(); i++) {
			if (i->isInside(p) && find(pi.C.begin(), pi.C.end(), *i) == pi.C.end())
				pi.C.push_back(*i);
		}

		if (pi.C.size() > 0 || pi.L.size() + pi.U.size() > 2)
		{
			intersections[p] = pi;
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
		Q.erase(--Q.end());
	}*/
	findIntersection(lines, intersections);
	for (auto i = intersections.begin(); i != intersections.end(); ) {
		if (!(i->second.C.size() > 0 || i->second.U.size() + i->second.L.size() > 2)) {
			auto j = i;
			i++;
			intersections.erase(j);
		}
		else {
			i++;
		}
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


int Yin::postion(Point c)
{
	if (spadjor.empty())
		return sign;
	bool res = false;
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		int b = i->postion(c);
		if (b == 2)
			return 2;
		res ^= (b == 1);
	}
	return (res ^ sign) ? 1 : 0;
}

int Yin::postion(Point c, Point d)
{
	Point center = 0.5 * (c + d);
	int r = postion(center);
	if (r != 2)
		return r;
	else
	{
		for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
			Polygon::Vertex* j = i->head;
			do {
				if (c == j->p && (d - c) * (j->next->p - c) > 0)
					return 2;
				j = j->next;
			} while (j != i->head);
		}
		return 3;
	}
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
	if ((postion(PL.head->p) == 0) ^ ori)
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

void Yin::resetSign()
{
	if (spadjor.empty())
		return;
	auto m = spadjor.begin()->head;
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		auto j = i->head;
		do {
			if (m->p < j->p)
				m = j;
			j = j->next;
		} while (j != i->head);
	}
	Point v1 = m->p - m->last->p;
	Point v2 = m->next->p - m->p;
	sign = v1.cross(v2) < 0;
}

bool Yin::checkPad()//检测是否两两几乎不交
{
	list<Line> lines;
	int id = 1;
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		Polygon::Vertex* ii = i->head;
		do {
			Line l(ii->p, ii->next->p);
			l.ID = id;
			lines.push_back(l);
			ii = ii->next;
		} while (ii != i->head);
		id++;
	}
	//判断是否有自交点或proper交点
	map<Point, PointInfo> intersections;
	map<double, bool> tiktak;
	findIntersection(lines, intersections);

	for (auto i = intersections.begin(); i != intersections.end(); ) {
		if (!(i->second.C.size() > 0 || i->second.U.size() + i->second.L.size() > 2)) {
			auto j = i;
			i++;
			intersections.erase(j);
		}
		else {
			i++;
		}
	}

	for (auto i = intersections.begin(); i != intersections.end(); i++) {
		list<Line> UL;
		UL.insert(UL.begin(), i->second.U.begin(), i->second.U.end());
		UL.insert(UL.begin(), i->second.L.begin(), i->second.L.end());
		for (auto j = UL.begin(); j != UL.end(); j++) {
			bool t = j->P == i->first;
			Point v = j->Q - j->P;
			double arg = t ? atan2(v.y, v.x) : atan2(-v.y, -v.x);
			if (tiktak.find(arg) != tiktak.end())
				return false;
			tiktak[arg] = t;
		}
		for (auto j = i->second.C.begin(); j != i->second.C.end(); j++) {
			Point v = j->Q - j->P;
			double arg1 = atan2(v.y, v.x);
			double arg2 = atan2(-v.y, -v.x);
			if (tiktak.find(arg1) != tiktak.end() || tiktak.find(arg2) != tiktak.end())
				return false;
			tiktak[arg1] = true;
			tiktak[arg2] = false;
		}
		//滴答测试
		bool status = tiktak.rbegin()->second;
		for (auto i = tiktak.begin(); i != tiktak.end(); i++) {
			if (i->second ^ status) {
				status = i->second;
			}
			else {
				return false;
			}
		}
	}

	return true;
}

bool Yin::checkOri()//检测曲线定向是否合法
{
	int n = spadjor.size();
	Yin cop = *this, null;
	cop.intersect(null);
	for (int i = 0; i < n; i++) {
		Polygon obj = *cop.spadjor.begin();
		cop.spadjor.pop_front();
		Polygon::Vertex* vp = obj.head;
		int pos = 0;
		do {
			Point p = 0.5 * (vp->p + vp->next->p);
			pos = cop.postion(p);
			vp = vp->next;
		} while (pos == 2);
		if (!(obj.orientation ^ (pos == 1)))
			return false;
		cop.spadjor.push_back(obj);
	}
	return true;
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

void Yin::OutPut(string filename)
{
	int n = 0;
	stringstream ss;
	fstream f(filename, ios::out);
	f << setprecision(16);
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		ss.str("");
		ss << filename << ++n << ".txt";

		Polygon::Vertex* j = i->head;
		do {
			f << j->p.x << " ";
			j = j->next;
		} while (j != i->head);
		f << j->p.x << endl;
		do {
			f << j->p.y << " ";
			j = j->next;
		} while (j != i->head);
		f << j->p.y << endl;

	}
	f.close();
}

void Yin::InPut(string filename)
{
	spadjor.clear();
	fstream f(filename, ios::in);
	stringstream ss;
	string str;
	while (getline(f, str))
	{
		ss.clear();
		ss.str(str);
		vector<double> xs, ys;
		double number;
		while (true) {
			ss >> number;
			if (ss.fail())
				break;
			xs.push_back(number);
		}
		getline(f, str);
		ss.clear();
		ss.str(str);
		while (true) {
			ss >> number;
			if (ss.fail())
				break;
			ys.push_back(number);
		}
		if (xs.size() != ys.size())
		{
			spadjor.clear();
			cout << "错误的文件格式" << endl;
			return;
		}
		else if (xs.size() == 0)
			return;
		spadjor.push_back(Polygon());
		auto Jor = --spadjor.end();
		for (int i = 0; i < xs.size() - 1; i++) {
			Jor->append(Point(xs[i], ys[i]));
		}
		Jor->refreshOri();
	}
	resetSign();
}

bool Yin::check(string& info)
{
	for (auto i = spadjor.begin(); i != spadjor.end(); i++) {
		if (!i->check()){
			info = "One or more polygons have self intersection(s)!";
			return false;
		}
	}
	if (!checkPad()){
		info = "There are two or more polygons with proper intersections!";
		return false;
	}
	if (!checkOri()) {
		info = "The orientation of polygons is illegal!";
		return false;
	}
	info = "The object is a well-defined Yin Set.";
	return true;
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
		spadjor.insert(spadjor.end(), atom.begin(), atom.end());
	}

}

Yin::Yin(string datafiles[], int num)
{
	load(datafiles, num);
}

Yin::~Yin()
{
}
