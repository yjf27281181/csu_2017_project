#pragma once
class AnyTest
{
public:
	AnyTest();
	~AnyTest();
	void readImage();
	void readLRImage();
private:
	float * leftImage;
	float * rightImage;
};

