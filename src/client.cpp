/*
 *  client.cpp
 *  Client application for talking to a server listening on port 5555
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/CLIENT
 *
 */

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>

#ifdef _WIN32

#include <windows.h>
#include <winsock.h>

#else

#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>

#endif

#define PORT		5555
#define SERVERHOST 	"127.0.0.1"
#define MAXMSG 1024
#ifndef _WIN32
	#define SOCKET int 
	#define INVALID_SOCKET (-1)
#endif

void cleanup_sockets()
{
#ifdef _WIN32 
	WSACleanup();
#endif
}

bool serr(int rc)
{
#ifdef _WIN32
	return rc==SOCKET_ERROR; 	//win
#else
	return rc<0; 				//osx
#endif
}

bool sinval(SOCKET sock)
{
#ifdef _WIN32
	return sock==INVALID_SOCKET;//win
#else
	return sock<0; 			//osx
#endif
}

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

int getErrorCode()
{
#ifdef _WIN32
	return WSAGetLastError();
#else
	return errno;
#endif
}

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

void init_sockaddr (struct sockaddr_in *name, const char *hostname, uint16_t port)
{
	//struct hostent *hostinfo;
	memset(name, 0, sizeof(sockaddr_in));
	name->sin_family = AF_INET;
	name->sin_port = htons (port);
	name->sin_addr.s_addr = inet_addr(hostname); 
}


void write_to_server (SOCKET filedes, char * message)
{
	int nbytes;

	nbytes = write (filedes, message, strlen (message) + 1);
	if (nbytes < 0)
	{
		perror ("write");
		exit (EXIT_FAILURE);
	}
}

void read_from_server(SOCKET filedes)
{
	int nbytes;
	char buf [MAXMSG+1];
	char buf2 [MAXMSG];
	
	nbytes=read(filedes, buf, MAXMSG);
	buf[nbytes]='\0';
	
	if(nbytes>=2 && buf[0]=='@' && buf[nbytes-1]=='\n'){
		for(int i=1; i<nbytes-1; i++){
			buf2[i-1]=buf[i];
		}
		buf2[nbytes-1]='\0';
		printf("Server: %s\n",buf2);
	}
	else{
		printf("Server sent corrupt response (\"%s\")\n",buf);
	}
}

int main (int argc, char * const argv[])
{
	char * message=NULL;
	if(argc>=2){
		std::string command="@";
		for(int i=1; i<argc; i++){
			if(i>1){
				command+="|";
			}
			command+=std::string(argv[i]);
		}
		command+="\n";
		message=new char [command.size()+1];
		strcpy(message,command.c_str());
	}
	else{
		printf("Supply a command as an argument!\n");
		return 1;
	}
	
	int sock;
	struct sockaddr_in servername;
	init_sockets();
	
	/* Create the socket.   */
	sock = socket (PF_INET, SOCK_STREAM, 0);
	if (sinval(sock)){
		perror ("socket (client)");
		return 1;
	}

	/* Connect to the server.   */
	init_sockaddr (&servername, SERVERHOST, PORT);
	int rc=connect (sock, (struct sockaddr *) &servername,	sizeof (servername));
	if(serr(rc)){
		printf("Error: connect failed (%d)\n",getErrorCode());
		return 1;
	}
	else{
		printf("Connected to %s\n",SERVERHOST);
	}

	/* Send data to the server.   */
	write_to_server (sock, message);
	
	read_from_server(sock);
	
	delete [] message;
	
	close (sock);
	return 0;;
}
