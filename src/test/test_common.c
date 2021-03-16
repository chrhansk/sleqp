#include "test_common.h"

#include <stdlib.h>
#include <string.h>

int run_suite(int argc, char *argv[], Suite* suite)
{
  char outfile[256];

  char* ret = strrchr(argv[0], '/');

  snprintf(outfile, 256, "%s.xml", ret ? (ret + 1) : argv[0]);

  int num_fails;

  SRunner* srunner;

  srunner = srunner_create(suite);

  if(getenv("TEST_EXPORT_XML"))
  {
    srunner_set_xml(srunner, outfile);
  }

  srunner_set_fork_status(srunner, CK_NOFORK);
  srunner_run_all(srunner, CK_NORMAL);

  num_fails = srunner_ntests_failed(srunner);

  srunner_free(srunner);

  return (num_fails > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
