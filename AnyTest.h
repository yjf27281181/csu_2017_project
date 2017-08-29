#pragma once
#include <vector>

#include "gdal_priv.h"
#include "cpl_conv.h" //for CPLMalloc()
#include "RPCProcessing.h"

class AnyTest
{
public:
	AnyTest();
	~AnyTest();
	void test();
	void readLRImage();
	void readRPCCOEFFCIENTAndRPCImAffine();
	void calCounterpartLines(SATPoint2D point, int layer);
private:
	float * leftImage;
	float * rightImage;

	int rImageSizeWidth;

	RPCCOEFFCIENT lCoefficient;
	RPCCOEFFCIENT rCoefficient;

	RPCImAffine lAffine;
	RPCImAffine rAffine;

	void calCoefficients(vector<SATPoint2D> points, double& k, double& b);


};

