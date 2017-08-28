#include "afx.h"
#include "RPCProcessing.h"
#include "malloc.h"
#include "math.h"

RPCProcessing::RPCProcessing(void)
{
}


RPCProcessing::~RPCProcessing(void)
{
}

//前方交会
void RPCProcessing::RPCInterSection(SATPoint2D *lpt,SATPoint2D *rpt,int npt,SATPoint3D *ptObj)
{

	double A[4][3],l[4][1],d[3][1],AT[3][4],ATA[3][3],ATl[3][1];
	d[0][0]=1;
	d[1][0]=1;
	d[2][0]=50;
	double Pn=0,Ln=0,Hn=0;
	double NumL=0,DenL=0,NumS=0,DenS=0;	
	double dNumLdLn=0.0,dNumLdPn=0.0,dNumLdHn=0.0;
	double dDenLdLn=0.0,dDenLdPn=0.0,dDenLdHn=0.0;
	double dNumSdLn=0.0,dNumSdPn=0.0,dNumSdHn=0.0;
	double dDenSdLn=0.0,dDenSdPn=0.0,dDenSdHn=0.0;  //上面四行用来记录正解法偏导数

	int numofiterative=0;
	double Ltmp,Stmp;	

	for(int i=0; i<npt ; i++)
	{
		ptObj[i].L=(m_lRPC.LONG_OFF+m_rRPC.LONG_OFF)/2;
		ptObj[i].P=(m_lRPC.LAT_OFF+m_rRPC.LAT_OFF)/2;
		ptObj[i].H=(m_lRPC.HEIGHT_OFF+m_rRPC.HEIGHT_OFF)/2;	
	}

	for(int i=0; i<npt; i++)
	{
		// while(numofiterative<30&&(fabs(d[0][0])>0.00001||fabs(d[1][0])>0.00001||fabs(d[2][0])>0.001))
		while(numofiterative<15)
		{	 
			//……………………处理左片………………	

			Pn=(ptObj[i].P-m_lRPC.LAT_OFF)/m_lRPC.LAT_SCALE;
			Ln=(ptObj[i].L-m_lRPC.LONG_OFF)/m_lRPC.LONG_SCALE;
			Hn=(ptObj[i].H-m_lRPC.HEIGHT_OFF)/m_lRPC.HEIGHT_SCALE; 

			NumL=getaccumulation1(m_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			DenL=getaccumulation1(m_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			NumS=getaccumulation1(m_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			DenS=getaccumulation1(m_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);

			dNumLdLn=getpartialderivativeofL1(m_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdPn=getpartialderivativeofP1(m_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdHn=getpartialderivativeofH1(m_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);

			dDenLdLn=getpartialderivativeofL1(m_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdPn=getpartialderivativeofP1(m_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdHn=getpartialderivativeofH1(m_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);

			dNumSdLn=getpartialderivativeofL1(m_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdPn=getpartialderivativeofP1(m_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdHn=getpartialderivativeofH1(m_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);

			dDenSdLn=getpartialderivativeofL1(m_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdPn=getpartialderivativeofP1(m_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdHn=getpartialderivativeofH1(m_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);

			A[0][0]=(m_lRPC.LINE_SCALE/m_lRPC.LONG_SCALE)*((dNumLdLn*DenL-dDenLdLn*NumL)/(DenL*DenL));
			A[0][1]=(m_lRPC.LINE_SCALE/m_lRPC.LAT_SCALE)*((dNumLdPn*DenL-dDenLdPn*NumL)/(DenL*DenL));
			A[0][2]=(m_lRPC.LINE_SCALE/m_lRPC.HEIGHT_SCALE)*((dNumLdHn*DenL-dDenLdHn*NumL)/(DenL*DenL));
			A[1][0]=(m_lRPC.SAMP_SCALE/m_lRPC.LONG_SCALE)*((dNumSdLn*DenS-dDenSdLn*NumS)/(DenS*DenS));
			A[1][1]=(m_lRPC.SAMP_SCALE/m_lRPC.LAT_SCALE)*((dNumSdPn*DenS-dDenSdPn*NumS)/(DenS*DenS));
			A[1][2]=(m_lRPC.SAMP_SCALE/m_lRPC.HEIGHT_SCALE)*((dNumSdHn*DenS-dDenSdHn*NumS)/(DenS*DenS));  

			Ltmp = ((NumL/DenL)*m_lRPC.LINE_SCALE+m_lRPC.LINE_OFF);
			Stmp = ((NumS/DenS)*m_lRPC.SAMP_SCALE+m_lRPC.SAMP_OFF);

			double L=m_lAffine.lineb0+m_lAffine.lineb1*Stmp+m_lAffine.lineb2*Ltmp;
			double S=m_lAffine.samplea0+m_lAffine.samplea1*Stmp+m_lAffine.samplea2*Ltmp;

			l[0][0]=lpt[i].line- L;                           
			l[1][0]=lpt[i].sample- S;

			//……………………处理右片………………………………

			Pn=(ptObj[i].P-m_rRPC.LAT_OFF)/m_rRPC.LAT_SCALE;
			Ln=(ptObj[i].L-m_rRPC.LONG_OFF)/m_rRPC.LONG_SCALE;		          
			Hn=(ptObj[i].H-m_rRPC.HEIGHT_OFF)/m_rRPC.HEIGHT_SCALE;  


			NumL=getaccumulation1(m_rRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			DenL=getaccumulation1(m_rRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			NumS=getaccumulation1(m_rRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			DenS=getaccumulation1(m_rRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);


			dNumLdLn=getpartialderivativeofL1(m_rRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdPn=getpartialderivativeofP1(m_rRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdHn=getpartialderivativeofH1(m_rRPC.LINE_NUM_COEFF, Ln, Pn, Hn);

			dDenLdLn=getpartialderivativeofL1(m_rRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdPn=getpartialderivativeofP1(m_rRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdHn=getpartialderivativeofH1(m_rRPC.LINE_DEN_COEFF, Ln, Pn, Hn);

			dNumSdLn=getpartialderivativeofL1(m_rRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdPn=getpartialderivativeofP1(m_rRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdHn=getpartialderivativeofH1(m_rRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);

			dDenSdLn=getpartialderivativeofL1(m_rRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdPn=getpartialderivativeofP1(m_rRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdHn=getpartialderivativeofH1(m_rRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);

			A[2][0]=(m_rRPC.LINE_SCALE/m_rRPC.LONG_SCALE)*((dNumLdLn*DenL-dDenLdLn*NumL)/(DenL*DenL));
			A[2][1]=(m_rRPC.LINE_SCALE/m_rRPC.LAT_SCALE)*((dNumLdPn*DenL-dDenLdPn*NumL)/(DenL*DenL));
			A[2][2]=(m_rRPC.LINE_SCALE/m_rRPC.HEIGHT_SCALE)*((dNumLdHn*DenL-dDenLdHn*NumL)/(DenL*DenL));
			A[3][0]=(m_rRPC.SAMP_SCALE/m_rRPC.LONG_SCALE)*((dNumSdLn*DenS-dDenSdLn*NumS)/(DenS*DenS));
			A[3][1]=(m_rRPC.SAMP_SCALE/m_rRPC.LAT_SCALE)*((dNumSdPn*DenS-dDenSdPn*NumS)/(DenS*DenS));
			A[3][2]=(m_rRPC.SAMP_SCALE/m_rRPC.HEIGHT_SCALE)*((dNumSdHn*DenS-dDenSdHn*NumS)/(DenS*DenS));


			Ltmp = ((NumL/DenL)*m_rRPC.LINE_SCALE+m_rRPC.LINE_OFF);
			Stmp = ((NumS/DenS)*m_rRPC.SAMP_SCALE+m_rRPC.SAMP_OFF);

			L=m_rAffine.lineb0+m_rAffine.lineb1*Stmp+m_rAffine.lineb2*Ltmp;
			S=m_rAffine.samplea0+m_rAffine.samplea1*Stmp+m_rAffine.samplea2*Ltmp;

			l[2][0]=rpt[i].line- L;                           
			l[3][0]=rpt[i].sample- S;


			for(int p=0; p<4; p++)
			{
				for(int q=0;q<3;q++)
				{
					AT[q][p]=A[p][q];					   
				}
			}

			matrixmulti(*AT,*A,*ATA,3,4,3);
			matrixmulti(*AT,*l,*ATl,3,4,1);


			if(inv(*ATA,3))
			{
				matrixmulti(*ATA,*ATl,*d,3,3,1);

			}

			else break;

			ptObj[i].L+=d[0][0];

			ptObj[i].P+=d[1][0];

			ptObj[i].H+=d[2][0];

			numofiterative++;

		}//end of while

		//lpoint<<i<<ptl.r[i]<<ptl.c[i]<<endl;
		numofiterative=0;
		d[0][0]=1;
		d[1][0]=1;
		d[2][0]=50;

	}//end of for


}
//3D 投影到2D

SATPoint2D RPCProcessing::RPCObj2Img(RPCcoeffcient &RPCcoef,SATPoint3D &ObjPt,RPCImAffine &affpara)
{
	SATPoint2D ptImg;
	double P = (ObjPt.P-RPCcoef.LAT_OFF)/RPCcoef.LAT_SCALE;
	double L = (ObjPt.L-RPCcoef.LONG_OFF)/RPCcoef.LONG_SCALE;
	double H = (ObjPt.H-RPCcoef.HEIGHT_OFF)/RPCcoef.HEIGHT_SCALE;
	
	double M[20];
	M[0]=1,M[1]=L,M[2]=P,M[3]=H,
	M[4]=L*P,M[5]=L*H,M[6]=P*H,M[7]=L*L,M[8]=P*P,M[9]=H*H,
	M[10]=P*L*H,M[11]=L*L*L,M[12]=L*P*P,M[13]=L*H*H,M[14]=L*L*P,
	M[15]=P*P*P,M[16]=P*H*H,M[17]=L*L*H,M[18]=P*P*H,M[19]=H*H*H;

	double NumL=0,DenL=0,NumS=0,DenS=0;
	for(int i=0; i<20; i++)
	{
		NumL += RPCcoef.LINE_NUM_COEFF[i]*M[i];
		DenL += RPCcoef.LINE_DEN_COEFF[i]*M[i];
		NumS += RPCcoef.SAMP_NUM_COEFF[i]*M[i];
		DenS += RPCcoef.SAMP_DEN_COEFF[i]*M[i];
	}

	double linetmp = (NumL/DenL)*RPCcoef.LINE_SCALE+RPCcoef.LINE_OFF;
	double sampetmp = (NumS/DenS)*RPCcoef.SAMP_SCALE+RPCcoef.SAMP_OFF;	

	ptImg.sample = affpara.samplea0+affpara.samplea1*sampetmp+affpara.samplea2*linetmp;
	ptImg.line = affpara.lineb0+affpara.lineb1*sampetmp+affpara.lineb2*linetmp;
	return ptImg;
}
//求逆矩阵
int RPCProcessing::inv(double *m1,int n)
{ 
	int *is,*js;
	int i,j,k,l,u,v;
	double temp,max_v;
	is=(int *)malloc(n*sizeof(int));
	js=(int *)malloc(n*sizeof(int));
	if(is==NULL||js==NULL)
	{
		printf("out of memory!\n");
		return(0);
	}
	for(k=0;k<n;k++)
	{
		max_v=0.0;
		for(i=k;i<n;i++)
			for(j=k;j<n;j++)
			{
				temp=fabs(m1[i*n+j]);
				if(temp>max_v)
				{
					max_v=temp; is[k]=i; js[k]=j;
				}
			}
			if(max_v==0.0)
			{
				free(is); free(js);
				printf("invers is not availble!\n");
				return(0);
			}
			if(is[k]!=k)
				for(j=0;j<n;j++)
				{
					u=k*n+j; v=is[k]*n+j;
					temp=m1[u]; m1[u]=m1[v]; m1[v]=temp;
				}
				if(js[k]!=k)
					for(i=0;i<n;i++)
					{
						u=i*n+k; v=i*n+js[k];
						temp=m1[u]; m1[u]=m1[v]; m1[v]=temp;
					}
					l=k*n+k;
					m1[l]=1.0/m1[l];
					for(j=0;j<n;j++)
						if(j!=k)
						{
							u=k*n+j;
							m1[u]*=m1[l];
						}
						for(i=0;i<n;i++)
							if(i!=k)
								for(j=0;j<n;j++)
									if(j!=k)
									{
										u=i*n+j;
										m1[u]-=m1[i*n+k]*m1[k*n+j];
									}
									for(i=0;i<n;i++)
										if(i!=k)
										{
											u=i*n+k;
											m1[u]*=-m1[l];
										}
	}
	for(k=n-1;k>=0;k--)
	{
		if(js[k]!=k)
			for(j=0;j<n;j++)
			{
				u=k*n+j; v=js[k]*n+j;
				temp=m1[u]; m1[u]=m1[v]; m1[v]=temp;
			}
			if(is[k]!=k)
				for(i=0;i<n;i++)
				{
					u=i*n+k; v=i*n+is[k];
					temp=m1[u]; m1[u]=m1[v]; m1[v]=temp;
				}
	}
	free(is); 
	free(js);
	return(1);
}

//矩阵相乘
void RPCProcessing::matrixmulti(double *r1,double *r2,double *r,int m,int n,int p)
{
	for(int i=0;i<m*p;i++)
		*(r+i)=0;
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<p;j++)
		{
			for(int k=0;k<n;k++)
			{
				*(r+i*p+j)+=*(r1+i*n+k)*(*(r2+k*p+j));
			}

		}
	}
}


//////////////////////////////////////////////////////////////////////////

double RPCProcessing::getpartialderivativeofL1(double *rpc, double L, double P, double H)
{
	return(rpc[1]+rpc[4]*P+rpc[5]*H+2*rpc[7]*L+rpc[10]*P*H+3*rpc[11]*L*L
		+rpc[12]*P*P+rpc[13]*H*H+2*rpc[14]*P*L+2*rpc[17]*L*H);
}
double RPCProcessing::getpartialderivativeofP1(double *rpc, double L, double P, double H)
{
	return(rpc[2]+rpc[4]*L+rpc[6]*H+2*rpc[8]*P+rpc[10]*L*H+2*rpc[12]*L*P
		+rpc[14]*L*L+3*rpc[15]*P*P+rpc[16]*H*H+2*rpc[18]*P*H);
}
double RPCProcessing::getpartialderivativeofH1(double *rpc, double L, double P, double H)
{
	return(rpc[3]+rpc[5]*L+rpc[6]*P+2*rpc[9]*H+rpc[10]*P*L+2*rpc[13]*L*H
		+2*rpc[16]*P*H+rpc[17]*L*L+rpc[18]*P*P+3*rpc[19]*H*H);
}
double RPCProcessing::getaccumulation1(double *rpc, double L, double P, double H)
{
	double M[20];
	double S=0;
	M[0]=1,M[1]=L,M[2]=P,M[3]=H,
		M[4]=L*P,M[5]=L*H,M[6]=P*H,M[7]=L*L,M[8]=P*P,M[9]=H*H,
		M[10]=P*L*H,M[11]=L*L*L,M[12]=L*P*P,M[13]=L*H*H,M[14]=L*L*P,
		M[15]=P*P*P,M[16]=P*H*H,M[17]=L*L*H,M[18]=P*P*H,M[19]=H*H*H;

	for(int j=0; j<20; j++)
	{
		S+=rpc[j]*M[j];
	}
	return S;

}
double RPCProcessing::getaccumulation1(double rpc[20], SATPoint3D &ObjPt)
{
	double M[20];
	double S=0;
	double L = ObjPt.L;
	double P = ObjPt.P;
	double H = ObjPt.H;
	M[0]=1,M[1]=L,M[2]=P,M[3]=H,
		M[4]=L*P,M[5]=L*H,M[6]=P*H,M[7]=L*L,M[8]=P*P,M[9]=H*H,
		M[10]=P*L*H,M[11]=L*L*L,M[12]=L*P*P,M[13]=L*H*H,M[14]=L*L*P,
		M[15]=P*P*P,M[16]=P*H*H,M[17]=L*L*H,M[18]=P*P*H,M[19]=H*H*H;

	for(int j=0; j<20; j++)
	{
		S+=rpc[j]*M[j];
	}
	return S;

}
void RPCProcessing::SetRPC(RPCcoeffcient &lRPCcoef,RPCcoeffcient &rRPCcoef)
{
	m_lRPC = lRPCcoef;
	m_rRPC = rRPCcoef;
	
}
void RPCProcessing::SetImgAffine(RPCImAffine &lRPCAff,RPCImAffine &rRPCAff)
{
	m_lAffine=lRPCAff;
	m_rAffine=rRPCAff;
}
//获取仿射变换的逆矩阵
// 	L0=lineb0+lineb1*L+lineb2*S      //右边地面坐标的投影
//  S0=samplea0+samplea1*S+samplea2*L //原始影像上的坐标
//  L = dstlineb0+dstlineb1*L0+dstlineb2*S0
//  S = dstsamplea0+dstsamplea1*L0+dstsamplea2*s0
bool RPCProcessing::GetInverseAffPara(RPCImAffine &srcAffPara,RPCImAffine &dstAffPara)
{
	double b0=srcAffPara.lineb0,b1=srcAffPara.lineb1,b2=srcAffPara.lineb2;
	double a0=srcAffPara.samplea0,a1=srcAffPara.samplea1,a2=srcAffPara.samplea2;

	double fenmu=a1*b2-a2*b1;

	if (fenmu == 0)
	{
		return false;
	}
	//
	//
	double inverA[6];
	//Line 相关
	inverA[0] = (a2*b0 - a0*b2)/fenmu;		
	inverA[1] = b2/fenmu;
	inverA[2] = a2/fenmu;


	//SAMPLE 相关
	inverA[3] = (a0*b1-a1*b0)/fenmu;
	inverA[4] = -b1/fenmu;
	inverA[5] = a1/fenmu;

	dstAffPara.samplea0 = inverA[0];
	dstAffPara.samplea1 = inverA[1];
	dstAffPara.samplea2 = inverA[2];

	dstAffPara.lineb0 = inverA[3];
	dstAffPara.lineb1 = inverA[4];
	dstAffPara.lineb2 = inverA[5];

	return true;

}

//zys
//data 2012-12-12
// Fs = NumS(P,L,H)-Sample*DenS(P,L,H) = 0;
// Fl = NumL(P,L,H)-Line*DenL(P,L,H)=0;
//求与P相关的偏导数
double RPCProcessing::GetPartialDerivativeofP(double Numrpc[20],double Denrpc[20],SATPoint3D &objpt, double SL)
{
	double P=objpt.P;
	double L=objpt.L;
	double H=objpt.H;
	double sum = Numrpc[2]+Numrpc[4]*L+Numrpc[6]*H+2*Numrpc[8]*P+
		         Numrpc[10]*L*H+2*Numrpc[12]*L*P+Numrpc[14]*L*L+3*Numrpc[15]*P*P+Numrpc[16]*H*H+
				 2*Numrpc[18]*P*H
				 -SL*(Denrpc[2]+Denrpc[4]*L+Denrpc[6]*H+2*Denrpc[8]*P+
				 Denrpc[10]*L*H+2*Denrpc[12]*L*P+Denrpc[14]*L*L+3*Denrpc[15]*P*P+Denrpc[16]*H*H+
				 2*Denrpc[18]*P*H);

	return sum;
}
//求与L相关的偏导数
double RPCProcessing::GetPartialDerivativeofL(double Numrpc[20],double Denrpc[20],SATPoint3D &objpt, double SL)
{
	double P=objpt.P;
	double L=objpt.L;
	double H=objpt.H;
	double sum = Numrpc[1]+Numrpc[4]*P+Numrpc[5]*H+2*Numrpc[7]*L+
		Numrpc[10]*P*H+3*Numrpc[11]*L*L+Numrpc[12]*P*P+Numrpc[13]*H*H+2*Numrpc[14]*L*P+
		2*Numrpc[17]*L*H
		-SL*( Denrpc[1]+Denrpc[4]*P+Denrpc[5]*H+2*Denrpc[7]*L+
		Denrpc[10]*P*H+3*Denrpc[11]*L*L+Denrpc[12]*P*P+Denrpc[13]*H*H+2*Denrpc[14]*L*P+
		2*Denrpc[17]*L*H);

	return sum;
}

//给定一个高程，以及对应的影像RPC参数，计算影像（sample，line）对应的 地面坐标

void RPCProcessing::RPCImg2Obj(RPCcoeffcient &RPCcoef,double H,RPCImAffine &affpara,SATPoint2D pimgpt,SATPoint3D &ObjPt)
{
	//坐标归一化

	//给定的仿射变换矩阵 的逆变换

	double sample = affpara.samplea0+affpara.samplea1*pimgpt.sample+affpara.samplea2*pimgpt.line;
	double line = affpara.lineb0+affpara.lineb1*pimgpt.sample+affpara.lineb2*pimgpt.line;

	sample = (sample-RPCcoef.SAMP_OFF)/RPCcoef.SAMP_SCALE;
	line = (line-RPCcoef.LINE_OFF)/RPCcoef.LINE_SCALE;

	ObjPt.L = 0;
	ObjPt.P = 0;
	
	ObjPt.H = (H-RPCcoef.HEIGHT_OFF)/RPCcoef.HEIGHT_SCALE;

	//

	//v=bx-lcoef
	double B[2][2],lcoef[2][1],x[2][1],BT[2][2],BTB[2][2],BTl[2][1];
	int noitrative = 0;

	//x[2][1], det P det L
	while (noitrative<15)
	{
		B[0][0] = GetPartialDerivativeofL(RPCcoef.SAMP_NUM_COEFF,RPCcoef.SAMP_DEN_COEFF,ObjPt,sample);
		B[0][1] = GetPartialDerivativeofP(RPCcoef.SAMP_NUM_COEFF,RPCcoef.SAMP_DEN_COEFF,ObjPt,sample);
		B[1][0] = GetPartialDerivativeofL(RPCcoef.LINE_NUM_COEFF,RPCcoef.LINE_DEN_COEFF,ObjPt,line);
		B[1][1] = GetPartialDerivativeofP(RPCcoef.LINE_NUM_COEFF,RPCcoef.LINE_DEN_COEFF,ObjPt,line);

		lcoef[0][0] = -(getaccumulation1(RPCcoef.SAMP_NUM_COEFF,ObjPt)-sample*getaccumulation1(RPCcoef.SAMP_DEN_COEFF,ObjPt));
		lcoef[1][0] = -(getaccumulation1(RPCcoef.LINE_NUM_COEFF,ObjPt)-line*getaccumulation1(RPCcoef.LINE_DEN_COEFF,ObjPt));

		for(int p=0; p<2; p++)
		{
			for(int q=0;q<2;q++)
			{
				BT[q][p]=B[p][q];					   
			}
		}

		matrixmulti(*BT,*B,*BTB,2,2,2);
		matrixmulti(*BT,*lcoef,*BTl,2,2,1);


		if(inv(*BTB,2))
		{
			matrixmulti(*BTB,*BTl,*x,2,2,1);

		}
		else break;

		ObjPt.L+=x[0][0];
		ObjPt.P+=x[1][0];
		/*if(fabs(x[0][0])<1.0e-10 && fabs(x[1][0])<1.0e-10)
			break;*/
		noitrative++;
	}

	ObjPt.H = H;
	ObjPt.L = ObjPt.L*RPCcoef.LONG_SCALE+RPCcoef.LONG_OFF;
	ObjPt.P = ObjPt.P*RPCcoef.LAT_SCALE+RPCcoef.LAT_OFF;
}


//给定一个高程，以及对应的影像RPC参数，计算影像（sample，line）对应的 地面坐标

bool RPCProcessing::CalDLT(SATPoint2D *p2dt,SATPoint3D *p3dt,int npt,double DLTpara[12])
{
	if (npt < 6)
	{
		return false;
	}

	double A[2][11],AT[11][2],l[2][1],x[11][1],ATA[11][11],ATl[11][1],SATA[11][11],SATl[11][1];

	for (int m1=0; m1<11; m1++ )
	{
		SATl[m1][0]= 0;
		for (int n1=0; n1<11;n1++)
		{
			SATA[m1][n1] = 0;
		}
	}

	A[0][4]=0,A[0][5]=0,A[0][6]=0,A[0][7]=0;
	A[1][0]=0,A[1][1]=0,A[1][2]=0,A[1][3]=0;	

	for (int i=0; i<npt; i++)
	{
		A[0][0] = p3dt[i].L,A[0][1]=p3dt[i].P,A[0][2]=p3dt[i].H,A[0][3]=1;
		A[0][8] = -p2dt[i].sample*p3dt[i].L,A[0][9]=-p2dt[i].sample*p3dt[i].P,A[0][10]=-p2dt[i].sample*p3dt[i].H;
		
	
		A[1][4] = p3dt[i].L,A[1][5]=p3dt[i].P,A[1][6]=p3dt[i].H,A[1][7]=1;
		A[1][8] = -p2dt[i].line*p3dt[i].L,A[1][9]=-p2dt[i].line*p3dt[i].P,A[1][10]=-p2dt[i].line*p3dt[i].H;
	
		l[0][0] = p2dt[i].sample;
		l[1][0] = p2dt[i].line;
	

		for (int i1=0; i1<2; i1++)
		{
			for (int j1=0; j1<11; j1++)
			{			
				AT[j1][i1] = A[i1][j1];
			}
		}
        
		matrixmulti(*AT,*l,*ATl,11,2,1);
        matrixmulti(*AT,*A,*ATA,11,2,11);
		for (int m2=0; m2<11; m2++ )
		{
			SATl[m2][0]+=ATl[m2][0];
			for (int n2=0; n2<11;n2++)
			{
				SATA[m2][n2] += ATA[m2][n2];
			}
		}
		
		
	}//end of for

	//********????*********** 
	    
    if(inv(*SATA,11))
	{	
		matrixmulti(*SATA,*SATl,*x,11,11,1);
	}
	else
	{
		return false;
	}
	for (int i=0; i<11; i++)
	{
		DLTpara[i] = x[i][0];
	}
	DLTpara[11] = 1;

	return true;	
}

bool RPCProcessing::readRPCfile(string &rpcfile,RPCcoeffcient &rpc)
{
	FILE *input = fopen(rpcfile.c_str(),"r");
	if (input == NULL)
	{
		return false;
	}
	char temp[40];
	fscanf(input, "%s %lf %s", &temp,&rpc.LINE_OFF,&temp); //fscanf(input, "%s", &temp);
	fscanf(input, "%s %lf %s", &temp, &rpc.SAMP_OFF, &temp);
	fscanf(input, "%s %lf %s", &temp, &rpc.LAT_OFF, &temp);
	fscanf(input, "%s %lf %s", &temp, &rpc.LONG_OFF, &temp);
	fscanf(input, "%s %lf %s", &temp, &rpc.HEIGHT_OFF, &temp);
	fscanf(input, "%s %lf %s", &temp, &rpc.LINE_SCALE, &temp);
	fscanf(input, "%s %lf %s", &temp, &rpc.SAMP_SCALE, &temp);
	fscanf(input, "%s %lf %s", &temp, &rpc.LAT_SCALE, &temp);
	fscanf(input, "%s %lf %s", &temp, &rpc.LONG_SCALE, &temp);
	fscanf(input, "%s %lf %s", &temp, &rpc.HEIGHT_SCALE, &temp);


	// fscanf(input,"%lf",&rpc.SAMP_OFF);
	//fscanf(input,"%lf",&rpc.LAT_OFF);
	//fscanf(input,"%lf",&rpc.LONG_OFF);
	//fscanf(input,"%lf",&rpc.HEIGHT_OFF);
	//fscanf(input,"%lf",&rpc.LINE_SCALE);
	//fscanf(input,"%lf",&rpc.SAMP_SCALE);
	//fscanf(input,"%lf",&rpc.LAT_SCALE);
	//fscanf(input,"%lf",&rpc.LONG_SCALE);
	//fscanf(input,"%lf",&rpc.HEIGHT_SCALE);

	for(int i=0;i<20;i++)
	{
		fscanf(input,"%s %le",&temp,&rpc.LINE_NUM_COEFF[i]);
	}
	for(int i=0;i<20;i++)
	{
		fscanf(input, "%s %le", &temp,&rpc.LINE_DEN_COEFF[i]);
	}
	for(int i=0;i<20;i++)
	{
		fscanf(input, "%s %le", &temp,&rpc.SAMP_NUM_COEFF[i]);
	}
	for(int i=0;i<20;i++)
	{
		fscanf(input, "%s %le", &temp,&rpc.SAMP_DEN_COEFF[i]);
	}
	fclose(input);

	return true;
}

bool RPCProcessing::readaffinepara(string &afffile,RPCImAffine &affpara)
{

	FILE *ffaaffine=fopen(afffile.c_str(),"r");
	char temp[40];
	if (ffaaffine == NULL)
	{
		return false;
	}
	fscanf(ffaaffine,"%s",temp);
	fscanf(ffaaffine, "%s %lf %lf %lf", &temp, &affpara.samplea0, &affpara.samplea1, &affpara.samplea2);
	fscanf(ffaaffine, "%s %lf %lf %lf", &temp, &affpara.lineb0, &affpara.lineb1, &affpara.lineb2);
	fclose(ffaaffine);

	return true;
}