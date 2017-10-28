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
	void createEpipolarImage();
private:
	float * leftImage;
	float * rightImage;


	RPCCOEFFCIENT lCoefficient;
	RPCCOEFFCIENT rCoefficient;

	RPCImAffine lAffine;
	RPCImAffine rAffine;

	int lineCount = 0;
	int sampleCount = 0;

	double lk = 0, lb = 0;
	double rk = 0, rb = 0;

	float* rightEpipolarImage;
	float* leftEpipolarImage;
	void calCoefficients(vector<SATPoint2D> points, double& k, double& b);


};

