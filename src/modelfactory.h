/*
 *  modelfactory.h
 *  Model DLL manager
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/MODEL
 *
 */ 
 
#ifndef modelfactory_h
#define modelfactory_h

#include <string>
#include <map>
#include <vector>

#include "model.h"

class iWQModelFactoryEntry
{
public:
	iWQModelFactoryMethod create;
	iWQModelDestructorMethod destroy;
};

//--------------------------------------------------------------------------------

typedef std::map<std::string, iWQModelFactoryEntry> iWQModelFactoryEntryMap;

//--------------------------------------------------------------------------------

class iWQModelFactory
{
private:
	int mPluginInterfaceMajorVersion;
	int mPluginInterfaceMinorVersion;
	iWQModelFactoryEntryMap mModelFactoryMethods;
	std::vector<void *> mLibraries;
public:
	iWQModelFactory(std::string pluginpath);
	~iWQModelFactory();
	iWQModel * newModelOfType(std::string type);
	void deleteModel(iWQModel * model);
};

#endif
