/*
 *  server.cpp
 *  Network functionality for keepalive running (works with Client)
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/SERVER
 *
 */ 

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>

#include "server.h"

#ifdef _WIN32

#include <windows.h>
#include <winsock2.h>

#else

#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#endif

#define MAXMSG	512
#define MAXPENDING 5

#ifndef _WIN32
	#define SOCKET int 
	#define INVALID_SOCKET (-1)
#endif

//-----------------------------------------------------------------------------------------------

void cleanup_sockets()
{
#ifdef _WIN32 
	WSACleanup();
#endif
}

//-----------------------------------------------------------------------------------------------

bool serr(int rc)
{
#ifdef _WIN32
	return rc==SOCKET_ERROR; 	//win
#else
	return rc<0; 				//osx
#endif
}

//-----------------------------------------------------------------------------------------------

bool sinval(SOCKET sock)
{
#ifdef _WIN32
	return sock==INVALID_SOCKET;//win
#else
	return sock<0; 			//osx
#endif
}

//-----------------------------------------------------------------------------------------------

void init_sockets()
{
#ifdef _WIN32
	WSADATA wsa;
	int e=WSAStartup(MAKEWORD(2,0),&wsa);
	if(e!=0){
		printf("[Error]: Failed to init winsock (error code=%d)\n",e);
		exit(EXIT_FAILURE);
	}
#endif
}

//-----------------------------------------------------------------------------------------------

int getErrorCode()
{
#ifdef _WIN32
	return WSAGetLastError();
#else
	return errno;
#endif
}

//-----------------------------------------------------------------------------------------------

// Replicate some Unix methods on Windows
#ifdef _WIN32
void close(SOCKET filedes)
{
	closesocket(filedes);
}

int write(SOCKET filedes, char * str, int length)
{
	return send(filedes, str, length, 0);
}

int read(SOCKET filedes, char * buffer, int size)
{
	return recv(filedes, buffer, size, 0);
}
#endif

//-----------------------------------------------------------------------------------------------

int read_from_client (SOCKET filedes)
{
	char echoBuffer[MAXMSG];        // Buffer for echo string
	int recvMsgSize;                // Size of received message

	// Receive message from client
	if ((recvMsgSize = read(filedes, echoBuffer, MAXMSG)) < 0){
		printf("[Error]: read() failed\n");
		return 1;
	}
	// Send received string and receive again until end of transmission 
	while (recvMsgSize > 0)      // zero indicates end of transmission
	{
		// Echo message back to client
		if (send(filedes, echoBuffer, recvMsgSize, 0) != recvMsgSize){
			printf("[Error]: send() failed\n");
			return 1;
		}
		// See if there is more data to receive
		if ((recvMsgSize = read(filedes, echoBuffer, MAXMSG)) < 0){
			printf("[Error]: read() failed\n");
			return 1;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------

SOCKET make_socket (uint16_t port)
{
	SOCKET sock;
	struct sockaddr_in name;
	
	init_sockets();
	
	// Create the socket.
	sock = socket (PF_INET, SOCK_STREAM, 0);
	if(sinval(sock)){
		printf("[Error]: Failed to create socket (%d)\n",getErrorCode());
		return -1;
	}
	
	// Give the socket a name.
	memset(&name,0,sizeof(sockaddr_in));
	name.sin_family = AF_INET;
	name.sin_port = htons (port);
	name.sin_addr.s_addr = htonl (INADDR_ANY);
	int rc=bind (sock, (struct sockaddr *) &name, sizeof (name));
	if(serr(rc)){
		printf("[Error]: Failed to bind to port %d (%d)\n",port, getErrorCode());
		return -1;
	}
	
	return sock;
}

//-----------------------------------------------------------------------------------------------

std::string makeCommandFancy(std::string command)
{
	std::string fancycmd=command;
	fancycmd=fancycmd.substr(1,fancycmd.size()-2);
	std::string::size_type pos=0;
	while(1){
		pos=fancycmd.find("|");
		if(pos==std::string::npos){
			break;
		}
		fancycmd.replace(pos,1," ");
	}
	return fancycmd;
}

//-----------------------------------------------------------------------------------------------

int iWQServer::run(int port, iWQProcessCallback func)
{
	// Code adopted from the WinSock tutorial at http://www.c-worker.ch/tuts/wstut_op.php
	SOCKET sock;
	//SOCKET connectedSocket;
	//int status;
	int rc;
	char buf[MAXMSG];
	fd_set active_fd_set, read_fd_set;
	int i;
	//struct sockaddr_in clientname;
	//int size;
	SOCKET clients [MAXPENDING];
	for(i=0;i<MAXPENDING;i++){
		clients[i]=INVALID_SOCKET;
	}
	
	// Create the socket and set it up to accept connections.
	sock = make_socket (port);
	if(sinval(sock)){
		printf("[Error]: make_socket() failed\n");
		return 1;
	}
	rc=listen (sock, MAXPENDING);
	if(serr(rc)){
		printf("[Error]: listen failed (%d)\n",getErrorCode());
		return 1;
	}
	else{
		printf("Listening on port #%d....\n",port); 
	}

	// Initialize the set of active sockets.
	FD_ZERO (&active_fd_set);
	FD_SET (sock, &active_fd_set);

	while (1){
		FD_ZERO(&read_fd_set); // Inhalt leeren
		FD_SET(sock,&read_fd_set); // Den Socket der verbindungen annimmt hinzuf¸gen
   
		// alle g¸ltigen client sockets hinzuf¸gen (nur die die nicht INVALID_SOCKET sind)
		for(i=0;i<MAXPENDING;i++){
			if(!sinval(clients[i])){
				FD_SET(clients[i],&read_fd_set);
			}
		}
		
		rc=select(FD_SETSIZE,&read_fd_set,NULL,NULL,NULL); // nicht vergessen den ersten parameter bei anderen betriebssystem anzugeben
		if(serr(rc)){
			printf("[Error]: select failed (%d)\n",getErrorCode());
			return 1;
		}
   
		// acceptSocket is im fd_set? => verbindung annehmen (sofern es platz hat)
		if(FD_ISSET(sock,&read_fd_set)) {
			// einen freien platz f¸r den neuen client suchen, und die verbingung annehmen
			for(i=0;i<MAXPENDING;i++){
				if(sinval(clients[i])){
					clients[i]=accept(sock,NULL,NULL);
					printf("*** Client %d connected ***\n",i);
					break;
				}
			}
		}
		
		// pr¸fen wlecher client sockets im fd_set sind
		for(i=0;i<MAXPENDING;i++){
			if(sinval(clients[i])){
				continue; // ung¸ltiger socket, d.h. kein verbunder client an dieser position im array
			}
			if(FD_ISSET(clients[i],&read_fd_set)){
				rc=recv(clients[i],buf,MAXMSG,0);
				// pr¸fen ob die verbindung geschlossen wurde oder ein fehler auftrat
				if(rc==0 || serr(rc)){
					printf("*** Client %d disconnected ***\n",i);
					close(clients[i]); // socket schliessen         
					clients[i]=INVALID_SOCKET; // seinen platz wieder freigeben
				}
				else{
					buf[rc]='\0';
					// daten ausgeben und eine antwort senden
					std::string fancycmd=makeCommandFancy(buf);
					printf("Client %d: %s\n",i,fancycmd.c_str());
					// antwort senden
					std::string result;
					if(func){
						result=func(buf);
						fancycmd=makeCommandFancy(result);
						printf("Me: %s\n",fancycmd.c_str());
					}
					send(clients[i],result.c_str(),result.size(),0);
				}
			}
		}
	}
	cleanup_sockets();
}

//-----------------------------------------------------------------------------------------------
