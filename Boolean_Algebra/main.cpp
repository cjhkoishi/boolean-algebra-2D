#include"bool_alg.h"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#define N 100

void drawPoint(cv::Mat img, Point P);
void drawLine(cv::Mat img, Line L);
void drawPolygon(cv::Mat img, Polygon PL);
void drawYin(cv::Mat img, Yin& Y, int mode);

void slideBar(int val, void*);


int value1 = 1;

int main()
{
	srand(unsigned(time(0)));
	cv::Mat img = cv::Mat::zeros(600, 800, CV_8UC3);
	Point pg[N], qg[N], rg[N];
	for (int i = 0; i < N; i++) {
		double arg = (i * 2.0) / N * 3.1415926;
		pg[i].x = 200 + (100 + 35 * cos(arg * 7)) * cos(arg);
		pg[i].y = 200 + (100 + 35 * cos(arg * 5)) * sin(arg);
	}
	for (int i = 0; i < N; i++) {
		double arg = -(i * 2.0) / N * 3.1415926;
		qg[i].x = 200 + (50) * cos(arg);
		qg[i].y = 200 + (50) * sin(arg);
	}
	for (int i = 0; i < N; i++) {
		//double arg = (i * 2.0 + 1) / N * 3.1415;
		double arg = (i * 2.0) / N * 3.1415926;
		rg[i].x = 160 + (100 + 90 * cos(arg * 4)) * cos(arg);
		rg[i].y = 160 + (100 + 90 * cos(arg * 4)) * sin(arg);
	}
	/*pg[0] = Point(100,100);
	pg[1] = Point(200, 200);
	pg[2] = Point(100, 300);
	pg[3] = Point(150, 200);

	qg[0] = Point(100, 100);
	qg[1] = Point(50, 200);
	qg[2] = Point(100, 300);
	qg[3] = Point(0, 200);

	rg[0] = Point(120, 90);
	rg[1] = Point(200, 90);
	rg[2] = Point(200, 210);
	rg[3] = Point(120, 210);*/

	Polygon PL1(pg, N);
	Polygon PL2(qg, N);
	Polygon PL3(rg, N);
	Yin Y1, Y2;
	Y1.spadjor.push_back(PL1);
	Y1.spadjor.push_back(PL2);
	Y2.spadjor.push_back(PL3);
	//Y1.intersect(Y2);
	//drawYin(img, Y1);

	//list<list<Point>> out;
	//Y1.cut(out);
	//Yin Y4 = Y1.inverse();
	Yin Y3 = Y1.meet(Y2);
	//Y4.spadjor.pop_back();
	//drawYin(img, Y4,0);
	//drawYin(img, Y2, 0);
	//drawYin(img, Y1, 0);
	//drawYin(img, Y2,0);
	drawYin(img, Y1, 0);
	drawYin(img, Y2, 0);

	cv::imshow("test", img);
	cv::createTrackbar("test1", "test", &value1, 100, slideBar);
	cv::waitKey(0);
	return 0;
}

void drawPoint(cv::Mat img, Point P)
{
	cv::circle(img, cv::Point((int)P.x, (int)P.y), 3, cv::Scalar(0, 0, 255), -1);
}

void drawLine(cv::Mat img, Line L)
{
	cv::line(img, cv::Point((int)L.P.x, (int)L.P.y), cv::Point((int)L.Q.x, (int)L.Q.y), cv::Scalar(255, 255, 255), 1, 16);
}

void drawPolygon(cv::Mat img, Polygon PL)
{
	Polygon::Vertex* i = PL.head;
	do {
		drawLine(img, Line(i->p, i->next->p));
		i = i->next;
	} while (i != PL.head);

}

void drawYin(cv::Mat img, Yin& Y, int mode)
{
	if (mode == 1)
		for (int i = 0; i < img.rows; i++) {
			for (int j = 0; j < img.cols; j++) {
				if (Y.interiorTest(Point(j, i)))
					img.at<cv::Vec3b>(i, j) = cv::Vec3b(128, 128, 128);
			}
		}
	else {
		for (auto i = Y.spadjor.begin(); i != Y.spadjor.end(); i++) {
			Polygon::Vertex* j = i->head;
			do {
				drawLine(img, Line(j->p, j->next->p));
				j = j->next;
			} while (j != i->head);
			do {
				if (j->isMarked)
					drawPoint(img, j->p);
				j = j->next;
			} while (j != i->head);
		}
	}
}

void slideBar(int val, void*)
{
	cv::Mat img = cv::Mat::zeros(600, 800, CV_8UC3);
	Point pg[N], qg[N], rg[N];
	for (int i = 0; i < N; i++) {
		double arg = (i * 2.0) / N * 3.1415926;
		pg[i].x = 200 + (100 + 35 * cos(arg * 7)) * cos(arg);
		pg[i].y = 200 + (100 + 35 * cos(arg * 5)) * sin(arg);
	}
	for (int i = 0; i < N; i++) {
		double arg = -(i * 2.0) / N * 3.1415926;
		qg[i].x = 200 + (50) * cos(arg);
		qg[i].y = 200 + (50) * sin(arg);
	}
	for (int i = 0; i < N; i++) {
		double arg = (i * 2.0) / N * 3.1415926;
		rg[i].x = 200 + (val - 50) * 2 + (100 + 90 * cos(arg * 4)) * cos(arg);
		rg[i].y = 200 + (val - 50) * 2 + (100 + 90 * cos(arg * 4)) * sin(arg);
	}

	Polygon PL1(pg, N);
	Polygon PL2(qg, N);
	Polygon PL3(rg, N);
	Yin Y1, Y2;
	Y1.spadjor.push_back(PL1);
	Y1.spadjor.push_back(PL2);
	Y2.spadjor.push_back(PL3);
	Yin Y4 = Y1.inverse();
	Yin Y3 = Y2.join(Y4);
	//drawYin(img, Y1, 0);
	//drawYin(img, Y2, 0);
	drawYin(img, Y3, 0);
	int b0, b1;
	Y3.getBettiNum(b0,b1);
	cout << b0 << " " << b1 << endl;
	cv::imshow("test", img);
}
