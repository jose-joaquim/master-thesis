#include "basic.h"

int main(int argc, char **argv) {
  freopen(argv[1], "r", stdin);
  read_data();

  T = 1;

  misi dummy;
  Solution ans = vrbsp(dummy);
  return 0;
}
