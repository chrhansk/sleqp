#include "lp/sleqp_lpi_soplex.h"

int main(int argc, char *argv[])
{
  SleqpLPi* lp_interface;
  sleqp_lpi_soplex_create_interface(&lp_interface);
  return 0;
}
