#include <string.h>
#include <stdlib.h>
int CSYSTEM(const char *line, int *size)
{
/* could not get it to work on CRAY without allocating string */
  char *string;
  int sts;
  string = (char *)malloc(*size+1);
  strcpy(string,line);
  sts = system(string);
  free((void *)string);
  return sts;
}
