/*
 *  modelfactory.cpp
 *  Model DLL manager
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/MODEL
 *
 */ 
 
#include "modelfactory.h"
#include <dirent.h>
#include <stdio.h>

//--------------------------------------------------------------------------------

#pragma mark Cross-platform dynamic loading 

#ifdef _WIN32
	#include <windows.h>
	#define PATH_DELIMITER "\\"
	#define PLUGIN_EXTENSION ".dll"
#else
	#include <dlfcn.h>
	#define PATH_DELIMITER "/"
	#define PLUGIN_EXTENSION ".so"
#endif

void * loadPlugin(std::string path)
{
#ifdef _WIN32
	HINSTANCE hdll = LoadLibrary(path.c_str());
	if(hdll <= (HINSTANCE)HINSTANCE_ERROR){
		return NULL;
	}
	return (void *)hdll;
#else
	return dlopen(path.c_str(), RTLD_LAZY);
#endif
}

void * loadMethod(void * library, std::string method_name)
{
#ifdef _WIN32
	HINSTANCE hdll = (HINSTANCE)library;
	if(hdll){
		return (void *)GetProcAddress(hdll, method_name.c_str());
	}
	return NULL;
#else
	return dlsym(library, method_name.c_str());
#endif
}

void unloadPlugin(void * library)
{
#ifdef _WIN32
	HINSTANCE hdll = (HINSTANCE)library;
	FreeLibrary(hdll);
#else
	dlclose(library);
#endif
}

//--------------------------------------------------------------------------------

#pragma mark ModelFactory

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

iWQModelFactory::iWQModelFactory(std::string pluginpath)
{
	//set preferred interface version
	mPluginInterfaceMajorVersion = IWQ_PLUGIN_INTERFACE_VERSION_MAJOR;
	mPluginInterfaceMinorVersion = IWQ_PLUGIN_INTERFACE_VERSION_MINOR;
	
	//load plugins
	//printf("Loading models...\n");
	
	DIR * d;
	struct dirent * dir;
	d = opendir(pluginpath.c_str());
	if(d){
		while ((dir = readdir(d)) != NULL){
			//determine if possible plugin
			std::string act_fname=dir->d_name;
			
			if(act_fname.find(PLUGIN_EXTENSION)!=std::string::npos){
				std::string fullpluginpath=pluginpath + PATH_DELIMITER + act_fname;
				//printf("\t%s\t\t", act_fname.c_str());
				
				/* load library */
				void * act_lib=loadPlugin(fullpluginpath);
				if(act_lib){
					//check interface version
					int imajversion=0;
					int iminversion=0;
					iWQPluginVersionMethod vfunc=(iWQPluginVersionMethod)loadMethod(act_lib, QUOTEME( IWQ_PLUGIN_MAJOR_VERSION_METHOD ));
					if(vfunc){
						imajversion=vfunc();
					}
					vfunc=(iWQPluginVersionMethod)loadMethod(act_lib, QUOTEME( IWQ_PLUGIN_MINOR_VERSION_METHOD ));
					if(vfunc){
						iminversion=vfunc();
					}
					if(imajversion==mPluginInterfaceMajorVersion){
						//should be OK, get model methods
						iWQModelFactoryMethod ffunc=(iWQModelFactoryMethod)loadMethod(act_lib, QUOTEME( IWQ_MODEL_CREATOR_METHOD ));
						iWQModelDestructorMethod dfunc=(iWQModelDestructorMethod)loadMethod(act_lib, QUOTEME( IWQ_MODEL_DESTROY_METHOD ));
						iWQModelIdentifierMethod ifunc=(iWQModelIdentifierMethod)loadMethod(act_lib, QUOTEME( IWQ_MODEL_IDENTIFIER_METHOD ));
						
						if(ffunc && dfunc && ifunc){
							iWQModelFactoryEntry e;
							e.create=ffunc;
							e.destroy=dfunc;
							std::string act_id=ifunc();
							if(mModelFactoryMethods.find(act_id)!=mModelFactoryMethods.end()){
                                printf("\t%s\t\t", act_fname.c_str());
								printf("Fail: Name \"%s\" is already loaded\n",act_id.c_str());
							}
							else{
								mModelFactoryMethods[act_id]=e;
								//printf("OK: (%s)\n",act_id.c_str());
							}
						}
						else{
                            printf("\t%s\t\t", act_fname.c_str());
							printf("Fail: Invalid plugin.\n");
						}							
					}
					else{
                        printf("\t%s\t\t", act_fname.c_str());
						printf("Fail: Invalid interface version (%d.%d instead of %d.%d)\n",imajversion, iminversion, mPluginInterfaceMajorVersion, mPluginInterfaceMinorVersion);
					}
					//store it anyway
					mLibraries.push_back(act_lib);
				}
				else{
					printf("Fail for %s\n",fullpluginpath.c_str());
				}
			}
		}
		closedir(d);
	}
	//printf("Ready.\n");
}

//--------------------------------------------------------------------------------

iWQModelFactory::~iWQModelFactory()
{
	//unload plugins
	for(int i=0; i<mLibraries.size(); i++){
		unloadPlugin(mLibraries[i]);
	}
}

//--------------------------------------------------------------------------------

iWQModel * iWQModelFactory::newModelOfType(std::string type)
{
	//get proper entry
	if(mModelFactoryMethods.find(type)==mModelFactoryMethods.end()){
		printf("[Error]: Cannot create model \"%s\" becuse it is of unknown type.\n", type.c_str());
		return NULL;
	}
	iWQModelFactoryEntry e = mModelFactoryMethods[type];
	if(e.create){
		return e.create();
	}
	else{
		return NULL;
	}
}

//--------------------------------------------------------------------------------

void iWQModelFactory::deleteModel(iWQModel * model)
{
	if(!model){
		return;
	}
	
	//ask the model about its type
	std::string type=model->modelType();
	
	//look for it in the respository
	if(mModelFactoryMethods.find(type)==mModelFactoryMethods.end()){
		printf("[Error]: Cannot destroy model \"%s\" becuse it is of unknown type.\n", type.c_str());
		return;
	}
	iWQModelFactoryEntry e = mModelFactoryMethods[type];
	if(e.destroy){
		e.destroy(model);
	}
}

//--------------------------------------------------------------------------------
