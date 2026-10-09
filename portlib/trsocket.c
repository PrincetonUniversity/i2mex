/*
 transp socket server and client
 */
#include <unistd.h>
#include <stdio.h>
#include <sys/socket.h>
#include <stdlib.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <string.h>
#include <errno.h>
#include <endian.h>
#include <bits/byteswap.h>
#include "f77name.h"
int PORTRANGE_LO = 8000;
int PORTRANGE_HI = 8100;

extern int errno ;

int port_not_in_range(int port)
{
  return !(PORTRANGE_LO <= port && port <= PORTRANGE_HI);
}

/*
 setup socket server
 */
int F77NAME(setup_trserver)(int* port)
{
    int server_fd, new_socket;
    struct sockaddr_in address;
    int opt = 1;
    int addrlen = sizeof(address);

    if(port_not_in_range(*port))
    {
      perror("port not in range");
      return -1;
    }

    /* Creating socket file descriptor */
    if ((server_fd = socket(AF_INET, SOCK_STREAM, 0)) == 0)
    {
        perror("socket failed");
        return -1;
    }

    /* Forcefully attaching socket to the port 8080 */
    if (setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT,
                                                  &opt, sizeof(opt)))
    {
        perror("setsockopt");
        return -1;
    }
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons( *port );

    /* Forcefully attaching socket to the port 8080 */
    if (bind(server_fd, (struct sockaddr *)&address,
                                 sizeof(address))<0)
    {
        perror("bind failed");
        return -1;
    }

    if (listen(server_fd, 3) < 0)
    {
        perror("listen");
        return -1;
    }

    if ((new_socket = accept(server_fd, (struct sockaddr *)&address,
                       (socklen_t*)&addrlen))<0)
    {
        perror("accept");
        return -1;
    }

    return new_socket;
}

/*
 * setup socket client
 */
int F77NAME(setup_trclient)(const char* ipaddr, int* socket_id, int* port)
{
    struct sockaddr_in address;
    int sock = 0, valread;
    struct sockaddr_in serv_addr;

    char *hello = "Hello from client";
    char buffer[1024] = {0};

    if(port_not_in_range(*port))
    {
      perror("port not in range");
      return -1;
    }

    if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
      printf("\n Socket creation error = %d\n", errno);
      return -1;
    }

    memset(&serv_addr, '0', sizeof(serv_addr));

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(*port);

    /* Convert IPv4 and IPv6 addresses from text to binary form */
    if(inet_pton(AF_INET, ipaddr, &serv_addr.sin_addr)<=0)
    {
       printf("\nInvalid address/ Address not supported \n");
       return -1;
    }

    if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
    {
       printf("\nConnection Failed \n");
       return -1;
    }

    printf(" sock = %d\n", sock);

  /*  send(sock , hello , strlen(hello) , 0 ); */
  /*  printf("Hello message sent\n");  */

  /*  valread = read( sock , buffer, 1024); */
  /*  printf("%s\n",buffer ); */

    (*socket_id) = sock;
    return 0;

}

int F77NAME(trsocket_close)(int* socket_fd)
{
 return close((*socket_fd));
}

int F77NAME(trsocket_string_to_double)(char* buffer, double* num)
{
  sscanf(buffer, "%lf", num);
  return 0;
}

int F77NAME(trsocket_double_to_string)(double* num, char*buffer)
{
 gcvt((*num), 32, buffer);
 strcat(buffer,"\n");
 return 0;
}

int F77NAME(trsocket_string_to_float)(char* buffer, float* num)
{
  sscanf(buffer, "%g", num);
  return 0;
}


int F77NAME(trsocket_float_to_string)(float* num, char*buffer)
{
   gcvt((*num), 16, buffer);
   strcat(buffer,"\n");
   printf(" %s\n", buffer);
   return 0;
}

int F77NAME(trsocket_int_to_string)(int* num, char* buffer)
{
   snprintf(buffer, sizeof(buffer), "%d", (*num));
   strcat(buffer,"\n");
   return 0;
}

int F77NAME(trsocket_string_to_int)(char* buffer, int* num)
{
 sscanf(buffer, "%d", num);
 return 0;
}

int F77NAME(trsocket_send_str)(char* buffer, int* socket_fd)
{
   if ((send((*socket_fd), buffer, strlen(buffer), 0))== -1) {
       fprintf(stderr, "Failure Sending Message\n");
       int ic = close((*socket_fd));
       return -1;
   }
   else {
    /*   printf("Message being sent: %s\n",buffer); */
       return 0;
  }
}

int F77NAME(trsocket_receive_str)(char* buffer, int* socket_fd)
{
  int rc;
  int i;
  char c[1];
  i = 0;
  do {
     rc = read((*socket_fd), c, 1);
     /* printf(" message received: %d\n", c[0]); */
     if (rc < 0) {
        if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
           /* use select() or epoll() to wait for the socket to be writable again */
        }
        else if (errno != EINTR) {
           return -1;
        }
     } else if((int) (c[0]) != 10) {
         buffer[i] = c[0];
         i++;
     }

  } while ((int)(c[0]) != 10);

  return 0;
}



/* send integer */
int F77NAME(trsocket_send_int)(int* num, int* fd)
{
    int32_t conv = htonl(*num);
    char *data = (char*)&conv;
    int left = sizeof(conv);
    int rc;
    do {
        rc = write((*fd), data, left);
        printf("message sending: %d\n", (*data));
        if (rc < 0) {
            if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
                /* use select() or epoll() to wait for the socket to be writable again */
            }
            else if (errno != EINTR) {
                return -1;
            }
        }
        else {
            data += rc;
            left -= rc;
        }
    }
    while (left > 0);
    return 0;
}

/* receive integer */
int F77NAME(trsocket_receive_int)(int *num, int* fd)
{
    int32_t ret;
    char buf[1];
    int left = sizeof(ret);
    int nsize = sizeof(ret);
    int rc;
    char *data = malloc(nsize*sizeof(char));

    do {
        rc = read((*fd), buf, 1);
     /* printf("message received: %d, left= %d\n", (*buffer), left); */
        if (rc <= 0) {
          if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
           /* use select() or epoll() to wait for the socket to be readable again */
          } else if (errno != EINTR) {
              return -1;
          }
        }
        else if((int) buf[0] != 10) {
            data[nsize-left] = buf[0];
            left -= 1;
        }
        else {
           left -= -1;
        }
    } while ((int)buf[0] != 10);

    sscanf(data, "%d", num);
  /* printf(" received int is: %d\n", (*num)); */
    return 0;
}


/* send integer */
int F77NAME(trsocket_send_float)(float* num, int* fd)
{
    char *data;
    sprintf(data, "%f", num);
    int left = sizeof(data);
    int rc;
    do {
        rc = write((*fd), data, left);
        if (rc < 0) {
            if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
             /* use select() or epoll() to wait for the socket to be writable again */
            }
            else if (errno != EINTR) {
                return -1;
            }
        }
        else {
            data += rc;
            left -= rc;
        }
    }
    while (left > 0);
    return 0;
}

/* receive float */
int F77NAME(trsocket_receive_float)(float *num, int* fd)
{
    float ret;
    char data[80]; /* = (char*)&ret; */
    int left;      /* = sizeof(ret); */
    int rc;
    left = 80;
    do {
        /* rc = read((*fd), data, left); */
        rc = recv((*fd), data, left, 0);
        if (rc <= 0) { /* instead of ret */
            if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
                /* use select() or epoll() to wait for the socket to be readable again */
            }
            else if (errno != EINTR) {
                return -1;
            }
        }
        else {
            left = 0;
        }
      printf("left = %d\n", left);
    }
    while (left > 0);
    sscanf(data, "%f", num);
    return 0;
}


/* send double float */
int F77NAME(trsocket_send_double)(double* num, int* fd)
{
    char *data;
    sprintf(data, "%f", num);
    int left = sizeof(data);
    int rc;
    do {
        rc = write((*fd), data, left);
        if (rc < 0) {
            if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
             /* use select() or epoll() to wait for the socket to be writable again */
            }
            else if (errno != EINTR) {
                return -1;
            }
        }
        else {
            data += rc;
            left -= rc;
        }
    }
    while (left > 0);
    return 0;
}

/* receive double float */
int F77NAME(trsocket_receive_double)(double *num, int* fd)
{
    double ret;
    char *data = (char*)&ret;
    int left = sizeof(ret);
    int rc;
    do {
        rc = read((*fd), data, left);
        if (rc <= 0) { /* instead of ret */
            if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
             /* use select() or epoll() to wait for the socket to be readable again */
            }
            else if (errno != EINTR) {
                return -1;
            }
        }
        else {
            data += rc;
            left -= rc;
        }
    }
    while (left > 0);
    sscanf(data, "%f", num);
    return 0;
}


