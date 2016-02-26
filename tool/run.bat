@echo off

REM Associate points to building element surfaces
C:\Python27\python.exe point_association.py input_data\Haus30.ifc input_data\Haus30\%%d.pcd.gz input_data\Haus30_alignment_matrix.json --use-display

REM Apply plane segmentation to the remaining unassociated points
md intermediate_files\Haus30\associated_points\segmented
bin\segmentation\segmentation intermediate_files\Haus30\associated_points\ 0 10000 0.1

REM Write some SPF files
C:\Python27\python.exe write_spf_file.py input_data\Haus30.ifc --all_points_as_unassociated
C:\Python27\python.exe write_spf_file.py input_data\Haus30.ifc 0     unassociated_cartesian continuous_param raster_grid
C:\Python27\python.exe write_spf_file.py input_data\Haus30.ifc 0     unassociated_cartesian discrete_param   raster_grid
C:\Python27\python.exe write_spf_file.py input_data\Haus30.ifc 0     unassociated_segmented discrete_param   raster_grid
C:\Python27\python.exe write_spf_file.py input_data\Haus30.ifc 0.020 unassociated_segmented discrete_param   raster_grid
C:\Python27\python.exe write_spf_file.py input_data\Haus30.ifc 0.100 unassociated_segmented discrete_param   raster_grid

REM Convert them into HDF files. The python script iterates over all ifc files and emits compressed and uncompressed hdf5 files
md generated_files\hdf\
C:\Python27\python.exe convert_to_hdf.py

REM Also write some files to compare against a baseline
copy input_data\Haus30.ifc generated_files\spf\
C:\Python27\python.exe write_spf_file.py input_data\Haus30.ifc --only_include_unassociated
C:\Python27\python.exe write_spf_file.py input_data\Haus30.ifc --only_include_unassociated_segmented
C:\Python27\python.exe convert_to_hdf.py
del generated_files\spf\Haus30.ifc

REM The python HDF5 serializer has been superseded by the IfcConvert executable
REM C:\Python27_x64\python.exe hdf5_serialization.py generated_files\spf\Haus30-no_grid-with_unassociated-raster_grid-discrete_param.ifc compressed single_precision

REM The viewer client has not been updated because of limitation in the h5py api in terms of variable length aggregates
REM C:\Python27_x64\python.exe client.py generated_files\spf\Haus30-no_grid-with_unassociated-raster_grid-discrete_param-compressed-single_precision.ifc
