#include "basic.h"

int main(int argc, char **argv) {
  load_instance(argv[1]);

  instance.T = 1;

  misi dummy;
  Solution ans = vrbsp(dummy);
  return 0;
}
