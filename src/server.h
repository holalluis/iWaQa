/*
 *  server.h
 *  Network functionality for keepalive running (works with Client)
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/SERVER
 *
 */ 

#ifndef server_h
#define server_h

typedef std::string (*iWQProcessCallback)(std::string command);

class iWQServer
{
public: 
	static int run(int port, iWQProcessCallback func);
};

#endif
