#include"bool_alg.h"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#define N 100

void drawPoint(cv::Mat img, Point P);
void drawLine(cv::Mat img, Line L);
void drawPolygon(cv::Mat img, Polygon PL);
void drawYin(cv::Mat img, Yin& Y, int mode);
Point transform(Point p);

void slideBar(int val, void*);

Yin Y1, Y2;

int value1 = 1;

int main() {
	cv::Mat img = cv::Mat::zeros(600, 800, CV_8UC3);
	Y2.InPut("Input/mickey.txt");
	Y1.InPut("Input/six_hole.txt");
	
	string info;
	Y1.check(info);
	cout << info << endl;
	Y2.check(info);
	cout << info << endl;

	//Y1.move(Point(0.5, 0));

	//Y1.intersect(Y2);
	//drawYin(img, Y2, 0);
	Yin Y3 = Y1.join(Y2);
	Y3.OutPut("res.txt");
	drawYin(img, Y3, 0);


	//Y4.check(info);
	//cout << info << endl;
	/*cout << Y1.checkOri() << endl;
	cout << Y2.checkOri() << endl;
	cout << Y3.checkOri() << endl;
	cout << Y4.checkPad() << endl;*/

	cv::imshow("test", img);
	cv::createTrackbar("test1", "test", &value1, 100, slideBar);
	cv::waitKey(0);
	return 0;
}

void drawPoint(cv::Mat img, Point P)
{
	P = transform(P);
	cv::circle(img, cv::Point((int)P.x, (int)P.y), 3, cv::Scalar(0, 0, 255), -1);
}

void drawLine(cv::Mat img, Line L)
{
	L.P = transform(L.P);
	L.Q = transform(L.Q);
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
				if (Y.postion(Point(j * 0.01, i * 0.01)) == 1)
					img.at<cv::Vec3b>(i, j) += cv::Vec3b(64, 64, 64);
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

Point transform(Point p)
{
	return Point((p.x) * 100, (p.y) * 100);
}

void slideBar(int val, void*)
{
	cv::Mat img = cv::Mat::zeros(600, 800, CV_8UC3);
	Y1.move(Point(0.01, 0));
	Yin Y3 = Y1.join(Y2);
	drawYin(img, Y3, 0);
	cv::imshow("test", img);
}
