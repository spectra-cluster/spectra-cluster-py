#################
id_transferer_cli
#################

.. automodule:: spectra_cluster.ui.id_transferer_cli

split_moff_file
---------------

This small script is necessary to use the moFF output files created by the
*id_transfer_cli* tool together with the moFF toolchain.
moFF expects one result file per input file. The
*id_transferer_cli* tool on the other hand creates one result file per clustering file.

The
**splot_moff_file** tool can now be used to create one result file per input file based
on a moFF formatted (use the option *--moff_compatible*) result file from the
*id_transferer_cli* tool.

.. automodule:: spectra_cluster.ui.split_moff_file