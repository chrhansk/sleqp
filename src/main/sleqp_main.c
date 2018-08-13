#include "sleqp_log.h"

int main(int argc, char *argv[])
{
  sleqp_log_debug("Hello %s", "debug");
  sleqp_log_info("Hello %s", "info");
  sleqp_log_warn("Hello %s", "warn");
  sleqp_log_error("Hello %s", "error");
  return 0;
}
