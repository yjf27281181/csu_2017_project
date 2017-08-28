#pragma once
class AnyTest
{
public:
	AnyTest();
	~AnyTest();
	void test();
	void readLRImage();
private:
	float * leftImage;
	float * rightImage;
};

