#include "AnyTest.h"
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
	//SATPoint2D test2Dpoint;
	//test2Dpoint.line = 10;
	//test2Dpoint.sample = 20;

	//SATPoint3D test3Dpoint;
	//processing.RPCImg2Obj(rpcCoeff, 10, rpcImAffine, test2Dpoint, test3Dpoint);

	//RPCImAffine inverseAffine;
	//processing.GetInverseAffPara(rpcImAffine, inverseAffine);

	//SATPoint2D point = processing.RPCObj2Img(rpcCoeff, test3Dpoint, inverseAffine);

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
	demFileData = (GDALDataset *)GDALOpen(demName, GA_ReadOnly);		//将tif文件存入进去
	int nImgSizeX1 = demFileData->GetRasterXSize();										//图片的长宽
	int nImgSizeY1 = demFileData->GetRasterYSize();

	leftImage = new float[nImgSizeX1*nImgSizeY1];											//用一个float类型的数组存储图片的灰度值
	demFileData->RasterIO(GF_Read, 0, 0, nImgSizeX1, nImgSizeY1, leftImage, nImgSizeX1, nImgSizeY1, GDT_Float32, 1, 0, 0, 0, 0);		//将图片数据存入数组

	//同上，读取右影像
	demName = "E:\\data\\2.tif";
	demFileData = (GDALDataset *)GDALOpen(demName, GA_ReadOnly);
	nImgSizeX1 = demFileData->GetRasterXSize();
	nImgSizeY1 = demFileData->GetRasterYSize();

	rightImage = new float[nImgSizeX1*nImgSizeY1];
	demFileData->RasterIO(GF_Read, 0, 0, nImgSizeX1, nImgSizeY1, rightImage, nImgSizeX1, nImgSizeY1, GDT_Float32, 1, 0, 0, 0, 0);
	rImageSizeWidth = nImgSizeX1;
}

//依据一系列points，拟合一个y=kx+b形式的直线
void AnyTest::calCoefficients(vector<SATPoint2D> points, double & k, double & b)
{
	double SigmaX(0);
	double SigmaY(0);
	double SigmaXY(0);
	double SigmaX2(0);
	int n = points.size();
	for (int i = 0; i < n; i++)
	{
		SigmaX += points[i].line;
		SigmaY += points[i].sample;
		SigmaXY += points[i].line * points[i].sample;
		SigmaX2 += points[i].line * points[i].line;
	}

	k = (n * SigmaXY - SigmaX * SigmaY) / (n * SigmaX2 - SigmaX * SigmaX);
	b = SigmaY / n - k * SigmaX / n;
}

/**
*读取相应的参数
**/
void AnyTest::readRPCCOEFFCIENTAndRPCImAffine()
{
	RPCProcessing processing;
	string path = "E:\\data\\testRPC.txt";
	processing.readRPCfile(path, lCoefficient); //读取左影像参数
	path = "E:\\data\\testAffine.txt";
	processing.readaffinepara(path, lAffine);

	path = "E:\\data\\testRPC.txt";
	processing.readRPCfile(path, rCoefficient); //读取左影像参数
	path = "E:\\data\\testAffine.txt";
	processing.readaffinepara(path, rAffine);
}

/**
*计算某个point在右影像上的直线
*point:需要的点
*layer:网格的层数
**/
void AnyTest::calCounterpartLines(SATPoint2D point, int layer)
{
	RPCProcessing processing;
	vector<SATPoint2D> rPoints;

	double maxHeight = 50;
	double minHeight = 0;
	for (int i = 0; i < layer; i++) {
		double tmpHeight = i* (maxHeight - minHeight) / (layer*1.0) + minHeight;	//根据层数选取相应的高度H
		SATPoint3D tmp3Dpoint;
		processing.RPCImg2Obj(lCoefficient, tmpHeight, lAffine, point, tmp3Dpoint);		//根据点的坐标和高度反算成实际的3维坐标
		
		RPCImAffine inverseRAffine;	//三维算二维是此参数需要取逆
		processing.GetInverseAffPara(rAffine, inverseRAffine);
		SATPoint2D tmpRPoint = processing.RPCObj2Img(rCoefficient, tmp3Dpoint, inverseRAffine);	//三维点换算到右影像的二维点
		rPoints.push_back(tmpRPoint);	//存入数组，以便计算直线
	}

	//计算右边对应点拟合的直线
	double rk = 0, rb = 0;
	calCoefficients(rPoints,rk,rb);
	
	SATPoint2D tmpRPoint2;//随机找一个右影像线上的点，这里是横坐标取全图的一半，纵坐标用y=kx+b计算
	tmpRPoint2.line = rImageSizeWidth / 2;
	tmpRPoint2.sample = rk*tmpRPoint2.line + rb;
	SATPoint3D tmp3Dpoint;
	processing.RPCImg2Obj(rCoefficient, 0, rAffine, tmpRPoint2, tmp3Dpoint);
	RPCImAffine inverseLAffine;
	processing.GetInverseAffPara(lAffine, inverseLAffine);
	SATPoint2D lPoint2 = processing.RPCObj2Img(lCoefficient, tmp3Dpoint, inverseLAffine); //左影像上的点

	vector<SATPoint2D> lPoints;
	lPoints.push_back(point);
	lPoints.push_back(lPoint2);
	double lk = 0, lb = 0;
	//根据两个点拟合成直线
	calCoefficients(lPoints, lk, lb);
}
