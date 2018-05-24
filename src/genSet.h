#pragma once
#ifndef __genSet__
#define __getSet__

#include <string>
#include <iostream>

class genSet
{
public:

	genSet() : online(false),exportObjs(false), RESOURCE_DIR(""), OUTPUT_DIR("") {};
	virtual ~genSet() {};

	bool online;
	bool exportObjs;
	bool exportTimings;
	std::string RESOURCE_DIR;
	std::string OUTPUT_DIR;

	void printGenSet() {
		std::cout << "General settings" << std::endl;
		std::cout << "	online: " << printGenBool(online) << std::endl;
		std::cout << "	exportObjs: " << printGenBool(exportObjs) << std::endl;
		std::cout << "	exportTimings: " << printGenBool(exportTimings) << std::endl;
		std::cout << "	RESOURCE_DIR: " << RESOURCE_DIR << std::endl;
		std::cout << "	OUTPUT_DIR: " << OUTPUT_DIR << std::endl;
	}

private:

	std::string printGenBool(bool b) {
		if (b) return "True";
		return "False";
	}

};

#endif