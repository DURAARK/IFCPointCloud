C:\Python27\python.exe point_association.py input_data\Haus30.ifc input_data\Haus30\%d.pcd.gz input_data\Haus30_alignment_matrix.json --use-display
C:\Python27\python.exe write_spf_file.py input_data\Haus30.ifc 0 with_unassociated discrete_param raster_grid
C:\Python27\python.exe hdf5_serialization.py generated_files\spf\Haus30-no_grid-with_unassociated-raster_grid-discrete_param.ifc compressed single_precision
C:\Python27\python.exe client.py generated_files\spf\Haus30-no_grid-with_unassociated-raster_grid-discrete_param-compressed-single_precision.ifc
