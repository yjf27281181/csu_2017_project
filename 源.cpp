#include "gdal_priv.h"
#include<iostream>
#include "AnyTest.h"
using namespace std;
int main()
{
	AnyTest anyTest;
	anyTest.test();
	anyTest.readLRImage();
	anyTest.readRPCCOEFFCIENTAndRPCImAffine();

	SATPoint2D test2D;
	test2D.line = 10;
	test2D.sample = 20;
	anyTest.calCounterpartLines(test2D,10);
	return 0;
}

