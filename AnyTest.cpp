#include "AnyTest.h"
#include "gdal_priv.h"
#include "cpl_conv.h" //for CPLMalloc()
#include "RPCProcessing.h"
//#include "cpl_conv.h" //for CPLMalloc()

typedef unsigned char BYTE;

using namespace std;



AnyTest::AnyTest()
{
}


AnyTest::~AnyTest()
{
}

/**
*测试函数，可以随便改
*/
void AnyTest::test()
{
	RPCProcessing processing;
	string path = "E:\\data\\testRPC.txt";
	RPCCOEFFCIENT rpcCoeff;
	processing.readRPCfile(path, rpcCoeff);

	path = "E:\\data\\testAffine.txt";
	RPCImAffine rpcImAffine;
	processing.readaffinepara(path, rpcImAffine);


	SATPoint2D test2Dpoint;
	test2Dpoint.line = 10;
	test2Dpoint.sample = 20;

	SATPoint3D test3Dpoint;
	processing.RPCImg2Obj(rpcCoeff, 10, rpcImAffine, test2Dpoint, test3Dpoint);

	RPCImAffine inverseAffine;
	processing.GetInverseAffPara(rpcImAffine, inverseAffine);

	SATPoint2D point = processing.RPCObj2Img(rpcCoeff, test3Dpoint, inverseAffine);

}
/**
*读取左右影像灰度值，并保存在float数组中
*/
void AnyTest::readLRImage()
{
	GDALDataset *demFileData;
	GDALDriver *poDriver1;
	GDALAllRegister();

	char* demName;
	demName = "E:\\data\\1.tif";
	demFileData = (GDALDataset *)GDALOpen(demName, GA_ReadOnly);
	int nImgSizeX1 = demFileData->GetRasterXSize();
	int nImgSizeY1 = demFileData->GetRasterYSize();

	leftImage = new float[nImgSizeX1*nImgSizeY1];
	demFileData->RasterIO(GF_Read, 0, 0, nImgSizeX1, nImgSizeY1, leftImage, nImgSizeX1, nImgSizeY1, GDT_Float32, 1, 0, 0, 0, 0);

	demName = "E:\\data\\2.tif";
	demFileData = (GDALDataset *)GDALOpen(demName, GA_ReadOnly);
	nImgSizeX1 = demFileData->GetRasterXSize();
	nImgSizeY1 = demFileData->GetRasterYSize();

	rightImage = new float[nImgSizeX1*nImgSizeY1];
	demFileData->RasterIO(GF_Read, 0, 0, nImgSizeX1, nImgSizeY1, rightImage, nImgSizeX1, nImgSizeY1, GDT_Float32, 1, 0, 0, 0, 0);
}
