#include "spec_model.h"
#include "model.h"
#include <string>

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

#ifdef _WIN32
	#define EXPORT extern "C" __declspec(dllexport)
#else
	#define EXPORT extern "C"
#endif

//################################################################################

// Plugin interface definition 1.0 

#define IWQ_PLUGIN_INTERFACE_VERSION_MAJOR 1
#define IWQ_PLUGIN_INTERFACE_VERSION_MINOR 0

#define IWQ_MODEL_CREATOR_METHOD iWQModelCreate
#define IWQ_MODEL_DESTROY_METHOD iWQModelDestroy
#define IWQ_MODEL_IDENTIFIER_METHOD iWQModelIdentifier
#define IWQ_PLUGIN_MAJOR_VERSION_METHOD iWQPluginMajorVersion
#define IWQ_PLUGIN_MINOR_VERSION_METHOD iWQPluginMinorVersion

//################################################################################

#ifdef IWQ_MODEL_NAME

EXPORT iWQModel * IWQ_MODEL_CREATOR_METHOD ()
{
	iWQModel * model = new IWQ_MODEL_NAME; 
	return model; 
}

EXPORT void IWQ_MODEL_DESTROY_METHOD (iWQModel * obj)
{
	delete obj;
}

EXPORT std::string IWQ_MODEL_IDENTIFIER_METHOD ()
{
	return QUOTEME( IWQ_MODEL_NAME );
}

EXPORT int IWQ_PLUGIN_MAJOR_VERSION_METHOD ()
{
	return IWQ_PLUGIN_INTERFACE_VERSION_MAJOR;
}

EXPORT int IWQ_PLUGIN_MINOR_VERSION_METHOD ()
{
	return IWQ_PLUGIN_INTERFACE_VERSION_MINOR;
}

#endif
