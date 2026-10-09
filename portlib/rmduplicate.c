/*
 get the number of non-dpulicated number in an array
 */
#include "f77name.h"

int F77NAME(rmduplicate)(arr,len)
int *arr;
int len;
{
  int prev = 0;
  int curr = 1;
  int last = len - 1;
  while (curr <= last) {
    for (prev = 0; prev < curr && arr[curr] != arr[prev]; ++prev);
    if (prev == curr) {
      ++curr;
    } else {
      arr[curr] = arr[last];
      --last;
    }
  }
  return curr;
}

