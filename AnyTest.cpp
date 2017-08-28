#include "AnyTest.h"
#include "gdal_priv.h"
#include "cpl_conv.h" //for CPLMalloc()
//#include "cpl_conv.h" //for CPLMalloc()

typedef unsigned char BYTE;

using namespace std;



AnyTest::AnyTest()
{
}


AnyTest::~AnyTest()
{
}

void AnyTest::readImage()
{
	GDALDataset *demFileData;
	GDALDriver *poDriver1;
	GDALAllRegister();

	char* demName;
	// 		char* demFirstName = "ADS_DEM_";
	demName = "E:\\data\\2.tif";

	demFileData = (GDALDataset *)GDALOpen(demName, GA_ReadOnly);
	const int nImgSizeX1 =500;
	const int nImgSizeY1 = 500;

	float *pDem = new float[nImgSizeX1*nImgSizeY1];
	demFileData->RasterIO(GF_Read, 0, 0, nImgSizeX1, nImgSizeY1, pDem, nImgSizeX1, nImgSizeY1, GDT_Float32, 1, 0, 0, 0, 0);

	int bandNum = demFileData->GetRasterCount();    //²¨¶ÎÊý
	int depth = GDALGetDataTypeSize(demFileData->GetRasterBand(1)->GetRasterDataType()) / 8;

	GDALDriver *pDriver = GetGDALDriverManager()->GetDriverByName("GTIFF"); //Í¼ÏñÇý¶¯
	char** ppszOptions = NULL;
	ppszOptions = CSLSetNameValue(ppszOptions, "BIGTIFF", "IF_NEEDED"); //ÅäÖÃÍ¼ÏñÐÅÏ¢
	const char* dstPath = "E:\\data\\2_part.tif";
	GDALDataset* dst = pDriver->Create(dstPath, nImgSizeX1, nImgSizeY1, bandNum, GDT_Byte, ppszOptions);
	
	//ÉêÇëbuf

	dst->RasterIO(GF_Write, 0, 0, nImgSizeX1, nImgSizeY1, pDem, nImgSizeX1, nImgSizeY1,
		GDT_Byte, bandNum, nullptr, bandNum*depth, nImgSizeX1*bandNum*depth, depth);

	delete[] pDem;

}

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
