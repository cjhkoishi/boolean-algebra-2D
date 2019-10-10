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

Yin Y1,Y2;
Polygon::Vertex* te;

int value1 = 1;

int main()
{
	srand(unsigned(time(0)));
	cv::Mat img = cv::Mat::zeros(600, 800, CV_8UC3);

	/*rg[0] = Point(120, 90);
	rg[1] = Point(200, 90);
	rg[2] = Point(200, 210);
	rg[3] = Point(120, 210);*/



	Polygon PL1("mickeydata/Data001.txt");
	Polygon PL2("mickeydata/Data002.txt");
	Polygon PL3("mickeydata/Data003.txt");
	Polygon PL4("mickeydata/Data004.txt");
	Polygon PL5("mickeydata/Data005.txt");
	Polygon PL6("mickeydata/Data006.txt");

	Polygon QL1("pandadata/Data001.txt");
	Polygon QL2("pandadata/Data002.txt");
	Polygon QL3("pandadata/Data003.txt");
	Polygon QL4("pandadata/Data004.txt");
	Polygon QL5("pandadata/Data005.txt");
	Polygon QL6("pandadata/Data006.txt");
	Polygon QL7("pandadata/Data007.txt");
	Polygon QL8("pandadata/Data008.txt");
	Polygon QL9("pandadata/Data009.txt");
	Polygon QL10("pandadata/Data010.txt");
	
	Y1.append(PL1);
	Y1.append(PL2);
	Y1.append(PL3);
	Y1.append(PL4);
	Y1.append(PL5);
	Y1.append(PL6);

	Y2.append(QL1);
	Y2.append(QL2);
	Y2.append(QL3);
	Y2.append(QL4);
	Y2.append(QL5);
	Y2.append(QL6);
	Y2.append(QL7);
	Y2.append(QL8);
	Y2.append(QL9);
	Y2.append(QL10);

	te = (Y2.spadjor.begin())->head;
	//drawYin(img, Y1);

	//list<list<Point>> out;
	//Y1.cut(out);
	//Yin Y4 = Y1;
	Yin Y3 = Y1.join(Y2);
	drawYin(img, Y1, 1);
	drawYin(img, Y2, 1);
	drawYin(img, Y3, 0);
	//drawPoint(img,Point(380.83679955096261, 152.66266035223680));
	//drawPoint(img, Point(308.51901102099629 , 154.40530934868514));
	//drawPoint(img, Point(385.24060265364523 , 145.56119825716991));
	
	//Y4.spadjor.pop_back();
	//drawYin(img, Y4,0);
	//drawYin(img, Y2, 0);
	//drawYin(img, Y1, 0);
	//drawYin(img, Y2,0);
	
	//drawYin(img, Y3, 0);

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

void slideBar(int val, void*)
{
	cv::Mat img = cv::Mat::zeros(600, 800, CV_8UC3);
	te = te->next;

	drawYin(img, Y2, 0);
	drawPoint(img,te->p);
	cv::imshow("test", img);
}
