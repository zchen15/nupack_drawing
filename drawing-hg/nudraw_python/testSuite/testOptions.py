# nudraw.py specific options
# Conrad Steenberg <conrad.steenberg@caltech.edu>
# Nov 20, 2008

import Options
# ------------------------------------------------------------------------------
# Parse command line options
def_option_vars={\
      #~ ""       : "",
      "runall"          : Options.IntList(),
      "suites"          : Options.StringList(),
      "list"            : Options.IntList(),
      "dryrun"          : Options.IntList(),
}

option_vars_doc={ \
      "runall"          : "Run all tests - WARNING MAY TAKE VERY LONG",
      "suites"          : "Comma separated list of suites to run",
      "list"            : "List available test suites",
      "dryrun"          : "Process all options but do not run tests"
}
