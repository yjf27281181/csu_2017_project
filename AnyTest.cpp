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
*���Ժ�������������
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
*��ȡ����Ӱ��Ҷ�ֵ����������float������
*/
void AnyTest::readLRImage()
{
	GDALDataset *demFileData;
	GDALDriver *poDriver1;
	GDALAllRegister();

	char* demName;
	demName = "E:\\data\\1.tif";
	demFileData = (GDALDataset *)GDALOpen(demName, GA_ReadOnly);		//��tif�ļ������ȥ
	int nImgSizeX1 = demFileData->GetRasterXSize();										//ͼƬ�ĳ���
	int nImgSizeY1 = demFileData->GetRasterYSize();

	leftImage = new float[nImgSizeX1*nImgSizeY1];											//��һ��float���͵�����洢ͼƬ�ĻҶ�ֵ
	demFileData->RasterIO(GF_Read, 0, 0, nImgSizeX1, nImgSizeY1, leftImage, nImgSizeX1, nImgSizeY1, GDT_Float32, 1, 0, 0, 0, 0);		//��ͼƬ���ݴ�������

	//ͬ�ϣ���ȡ��Ӱ��
	demName = "E:\\data\\2.tif";
	demFileData = (GDALDataset *)GDALOpen(demName, GA_ReadOnly);
	nImgSizeX1 = demFileData->GetRasterXSize();
	nImgSizeY1 = demFileData->GetRasterYSize();

	rightImage = new float[nImgSizeX1*nImgSizeY1];
	demFileData->RasterIO(GF_Read, 0, 0, nImgSizeX1, nImgSizeY1, rightImage, nImgSizeX1, nImgSizeY1, GDT_Float32, 1, 0, 0, 0, 0);
	rImageSizeWidth = nImgSizeX1;
}

//����һϵ��points�����һ��y=kx+b��ʽ��ֱ��
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
*��ȡ��Ӧ�Ĳ���
**/
void AnyTest::readRPCCOEFFCIENTAndRPCImAffine()
{
	RPCProcessing processing;
	string path = "E:\\data\\testRPC.txt";
	processing.readRPCfile(path, lCoefficient); //��ȡ��Ӱ�����
	path = "E:\\data\\testAffine.txt";
	processing.readaffinepara(path, lAffine);

	path = "E:\\data\\testRPC.txt";
	processing.readRPCfile(path, rCoefficient); //��ȡ��Ӱ�����
	path = "E:\\data\\testAffine.txt";
	processing.readaffinepara(path, rAffine);
}

/**
*����ĳ��point����Ӱ���ϵ�ֱ��
*point:��Ҫ�ĵ�
*layer:����Ĳ���
**/
void AnyTest::calCounterpartLines(SATPoint2D point, int layer)
{
	RPCProcessing processing;
	vector<SATPoint2D> rPoints;

	double maxHeight = 50;
	double minHeight = 0;
	for (int i = 0; i < layer; i++) {
		double tmpHeight = i* (maxHeight - minHeight) / (layer*1.0) + minHeight;	//���ݲ���ѡȡ��Ӧ�ĸ߶�H
		SATPoint3D tmp3Dpoint;
		processing.RPCImg2Obj(lCoefficient, tmpHeight, lAffine, point, tmp3Dpoint);		//���ݵ������͸߶ȷ����ʵ�ʵ�3ά����
		
		RPCImAffine inverseRAffine;	//��ά���ά�Ǵ˲�����Ҫȡ��
		processing.GetInverseAffPara(rAffine, inverseRAffine);
		SATPoint2D tmpRPoint = processing.RPCObj2Img(rCoefficient, tmp3Dpoint, inverseRAffine);	//��ά�㻻�㵽��Ӱ��Ķ�ά��
		rPoints.push_back(tmpRPoint);	//�������飬�Ա����ֱ��
	}

	//�����ұ߶�Ӧ����ϵ�ֱ��
	double rk = 0, rb = 0;
	calCoefficients(rPoints,rk,rb);
	
	SATPoint2D tmpRPoint2;//�����һ����Ӱ�����ϵĵ㣬�����Ǻ�����ȡȫͼ��һ�룬��������y=kx+b����
	tmpRPoint2.line = rImageSizeWidth / 2;
	tmpRPoint2.sample = rk*tmpRPoint2.line + rb;
	SATPoint3D tmp3Dpoint;
	processing.RPCImg2Obj(rCoefficient, 0, rAffine, tmpRPoint2, tmp3Dpoint);
	RPCImAffine inverseLAffine;
	processing.GetInverseAffPara(lAffine, inverseLAffine);
	SATPoint2D lPoint2 = processing.RPCObj2Img(lCoefficient, tmp3Dpoint, inverseLAffine); //��Ӱ���ϵĵ�

	vector<SATPoint2D> lPoints;
	lPoints.push_back(point);
	lPoints.push_back(lPoint2);
	double lk = 0, lb = 0;
	//������������ϳ�ֱ��
	calCoefficients(lPoints, lk, lb);
}
