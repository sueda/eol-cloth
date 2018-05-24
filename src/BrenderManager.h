/*
* @author: Gustavo Lopez 10-21-17
*
* @version: 1.0
*/

#pragma once

#include <vector>
#include <memory>
#include <string>
#include <ostream>

class Brenderable;

class BrenderManager
{
private:
	static bool instanceFlag;
	static BrenderManager *manager;
	int frame;
	std::string EXPORT_DIR;
	std::vector<std::shared_ptr<Brenderable> > brenderables;
	BrenderManager()
	{
		//private constructor
		EXPORT_DIR = ".";
		frame = 0;
	}
public:
	static BrenderManager* getInstance();
	void setExportDir(std::string export_dir);
	int getFrame() const;
	void exportBrender(double time = 0.0);
	void add(std::shared_ptr<Brenderable> brenderable);
	~BrenderManager()
	{
		instanceFlag = false;
	}
};