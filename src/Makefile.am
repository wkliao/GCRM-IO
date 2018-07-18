#
#  Copyright (C) 2012, Northwestern University.
#  See COPYRIGHT notice COPYING in top-level directory.
#
# @configure_input@

noinst_PROGRAMS = gcrm_io

gcrm_io_SOURCES = gcrm_param_init.c \
                  read_param.c  \
                  init_term.c \
                  ZGrd_output.c \
                  ZGrd_params.c \
                  ZGrd_params_vertical.c \
                  ZGrd_params_tracer.c \
                  ZGrd_params_time.c \
                  ZGrd_vars_diagnostic.c \
                  ZGrd_vars_prognostic.c \
                  ZGrd_clusdet.c \
                  ZGrd_register.c  \
                  phys_vars_diag_out.c  \
                  grid_subdomain.c \
                  grid_params.c \
                  grid_connectivity.c \
                  grid_metrics.c \
                  grid_utilities.c \
                  gio_API_statistics.c \
                  wrap_data.c \
                  dimensions.c \
                  gio_driver.c \
                  gio_sched.c \
                  gio_register.c  \
                  gio_utility.c \
                  average.c \
                  io_netcdfx.c \
                  nonblocking_io.c \
                  copy.c \
                  util.c \
                  main.c

# Pack these header files in distribution
EXTRA_DIST = ZGrd_clusdet.h \
             ZGrd_output.h \
             ZGrd_params.h \
             ZGrd_params_time.h \
             ZGrd_params_tracer.h \
             ZGrd_params_vertical.h \
             ZGrd_vars_diagnostic.h \
             ZGrd_vars_prognostic.h \
             gcrm.h \
             gio.h \
             gio_data.h \
             gio_defaults.h \
             gio_dim.h \
             gio_file.h \
             gio_grid.h \
             gio_prop.h \
             grid_metrics.h \
             grid_params.h \
             grid_subdomain.h \
             grid_utilities.h \
             phys_vars_diag_out.h \
             physical_params.h \
             util.h \
             wrap_data.h

dist-hook:
	$(SED_I) -e "s|RELEASE_DATE|`date '+%e %b %Y'`|g" $(distdir)/main.c

