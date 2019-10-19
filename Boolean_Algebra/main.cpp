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

int main3() {
	cv::Mat img = cv::Mat::zeros(600, 800, CV_8UC3);
	Polygon PL1, PL2, PL3;
	PL1.append(Point(0, 100));
	PL1.append(Point(100, 100));
	PL1.append(Point(100, 200));
	PL1.append(Point(0, 200));

	PL2.append(Point(300, 0));
	PL2.append(Point(500, 200));
	PL2.append(Point(300, 400));
	PL2.append(Point(100, 200));

	PL3.append(Point(300, 100));
	PL3.append(Point(200, 200));
	PL3.append(Point(300, 300));
	PL3.append(Point(400, 200));
	Y1.append(PL1);
	Y2.append(PL2);
	Y2.append(PL3);

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			int r = Y2.postion(Point(j, i));
			switch (r)
			{
			case 0:
				img.at<cv::Vec3b>(i, j) += cv::Vec3b(0, 0, 0);
				break;
			case 1:
				img.at<cv::Vec3b>(i, j) += cv::Vec3b(64, 64, 64);
				break;
			case 2:
				img.at<cv::Vec3b>(i, j) += cv::Vec3b(0, 0, 255);
				break;
			default:
				break;
			}
		}
	}

	cv::imshow("test", img);
	cv::waitKey(0);

	return 0;
}

int main1() {

	cv::Mat img = cv::Mat::zeros(600, 800, CV_8UC3);

	Polygon PL1, PL2, PL3;
	PL1.append(Point(0, 10));
	PL1.append(Point(100, 10));
	PL1.append(Point(100, 210));
	PL1.append(Point(0, 210));

	PL2.append(Point(100, 10));
	PL2.append(Point(400, 10));
	PL2.append(Point(400, 310));
	PL2.append(Point(100, 310));

	PL3.append(Point(200, 110));
	PL3.append(Point(200, 230));
	PL3.append(Point(300, 210));
	PL3.append(Point(300, 110));
	Y1.append(PL1);
	Y2.append(PL2);
	Y2.append(PL3);

	Y1.move(Point(196, 0));
	Yin Y3 = Y1.join(Y2);
	//Y1.intersect(Y2);
	drawYin(img, Y3, 0);

	cv::imshow("test", img);
	cv::createTrackbar("test1", "test", &value1, 100, slideBar);
	cv::waitKey(0);
	return 0;
}

int main0() {

	cv::Mat img = cv::Mat::zeros(600, 800, CV_8UC3);

	Polygon PL1, PL2, PL3;
	PL1.append(Point(0.5, 0));
	PL1.append(Point(0.5, 0.5));
	PL1.append(Point(1, 0.5));
	PL1.append(Point(1, 1));
	PL1.append(Point(0, 1));
	PL1.append(Point(0, 0));


	PL2.append(Point(0.5, 0));
	PL2.append(Point(1, 0));
	PL2.append(Point(1, 0.5));
	PL2.append(Point(0.5, 0.5));

	Y1.append(PL1);
	Y2.append(PL2);

	Yin Y3 = Y1.join(Y2);
	//Y1.intersect(Y2);
	drawYin(img, Y3, 0);

	cv::imshow("test", img);
	cv::createTrackbar("test1", "test", &value1, 100, slideBar);
	cv::waitKey(0);
	return 0;
}

int main5()
{
	cv::Mat img = cv::Mat::zeros(600, 800, CV_8UC3);
	stringstream ss;
	string s1[6], s2[10];
	for (int i = 0; i < 6; i++) {
		ss.str("");
		ss << "mickeydata/Data_" << i + 1 << ".txt";
		s1[i] = ss.str();
	}
	for (int i = 0; i < 10; i++) {
		ss.str("");
		ss << "pandadata/Data_" << i + 1 << ".txt";
		s2[i] = ss.str();
	}
	Y1.load(s1, 6);
	Y2.load(s2, 10);
	Yin Y3 = Y1.join(Y2);
	Yin Y4;
	Y4.InPut("panda_inv.txt");
	//Y1.intersect(Y2);
	//drawYin(img, Y1, 1);
	//drawYin(img, Y2, 1);

	//drawYin(img, Y3, 0);
	//Y1.OutPut("mickey.txt");
	//Y2.OutPut("panda.txt");
	drawYin(img, Y4, 1);

	cv::imshow("test", img);
	cv::createTrackbar("test1", "test", &value1, 100, slideBar);
	cv::waitKey(0);
	return 0;
}

int main() {
	cv::Mat img = cv::Mat::zeros(600, 800, CV_8UC3);
	Y1.InPut("Input/two_inv_square.txt");
	Y2.InPut("Input/imp.txt");
	Y1.move(Point(1.0000001,0));
	Yin Y3 = Y1.join(Y2);
	Y3.OutPut("res.txt");
	drawYin(img,Y3,0);

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
				if (Y.postion(Point(j*0.01, i*0.01)) == 1)
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
