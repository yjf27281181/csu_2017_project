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

	anyTest.createEpipolarImage();
	return 0;
}

