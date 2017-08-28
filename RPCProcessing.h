#pragma once
#include <string>


//class rpc
//author zys 
//data 2010
using namespace std;


typedef struct RPCCOEFFCIENT
{
	double LINE_OFF;
	double SAMP_OFF;
	double LAT_OFF;
	double LONG_OFF;
	double HEIGHT_OFF;
	double LINE_SCALE;
	double SAMP_SCALE;
	double LAT_SCALE;
	double LONG_SCALE;
	double HEIGHT_SCALE;
	double LINE_NUM_COEFF[20];
	double LINE_DEN_COEFF[20];
	double SAMP_NUM_COEFF[20];
	double SAMP_DEN_COEFF[20];	
}RPCcoeffcient;


//利用RPC 投影，还要使用以下仿射变换变换到影像空间才能获取对应像点的坐标
struct RPCImAffine
{
	//L=lineb0+lineb1*L+lineb2*S
	//S=samplea0+samplea1*S+samplea2*L
	double samplea0; //shift parameters
	double samplea1;
	double samplea2;

	double lineb0; //shift parameters
	double lineb1;
	double lineb2;
};
//


struct SATPoint2D
{
	// ------>sample
	// |
	// |
	// |
	// line

	double sample,line;
};

struct SATPoint3D
{
	//
	double L,P,H;
};

class RPCProcessing
{

public:
	RPCProcessing(void);
	~RPCProcessing(void);
	

	//设置左右影像的 RPC 参数
	void SetRPC(RPCcoeffcient &lRPCcoef,RPCcoeffcient &rRPCcoef);
	void SetImgAffine(RPCImAffine &lRPCAff,RPCImAffine &rRPCAff);

	//前方交会
	void RPCInterSection(SATPoint2D *lpt,SATPoint2D *rpt,int npt,SATPoint3D *ptObj);
	//3D 投影到2D
	SATPoint2D RPCObj2Img(RPCcoeffcient &RPCcoef,SATPoint3D &ObjPt,RPCImAffine &affpara);

	//给定地物点高程确定地面点坐标
	void RPCImg2Obj(RPCcoeffcient &RPCcoef,double H,RPCImAffine &affpara,SATPoint2D pimgpt,SATPoint3D &ObjPt);
	//输入一个仿射变换参数 求出逆仿射变换参数
	bool GetInverseAffPara(RPCImAffine &srcAffPara,RPCImAffine &dstAffPara);

	bool readRPCfile(string &rpcfile,RPCcoeffcient &rpc);
	bool readaffinepara(string &afffile,RPCImAffine &affpara);

	bool CalDLT(SATPoint2D *p2dt,SATPoint3D *p3dt,int npt,double DLTpara[12]);
private:
	int inv(double *m1,int n);
	void matrixmulti(double *r1,double *r2,double *r,int m,int n,int p);

	//偏导数
	double getpartialderivativeofL1(double *rpc, double L, double P, double H);
	double getpartialderivativeofP1(double *rpc, double L, double P, double H);
	double getpartialderivativeofH1(double *rpc, double L, double P, double H);
	//sum
	double getaccumulation1(double *rpc, double L, double P, double H);
	double getaccumulation1(double rpc[20], SATPoint3D &ObjPt);
	//像方到物方
	double GetPartialDerivativeofP(double Numrpc[20],double Denrpc[20],SATPoint3D &objpt, double SL);
	double GetPartialDerivativeofL(double Numrpc[20],double Denrpc[20],SATPoint3D &objpt, double SL);



	RPCcoeffcient m_lRPC,m_rRPC;
	RPCImAffine m_lAffine,m_rAffine;
};

