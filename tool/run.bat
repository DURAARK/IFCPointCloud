@echo off

C:\Python27\python.exe point_association.py input_data\Haus30.ifc input_data\Haus30\%%d.pcd.gz input_data\Haus30_alignment_matrix.json --use-display

REM Apply plane segmentation to the remaining unassociated points
md intermediate_files\Haus30\associated_points\segmented
bin\segmentation\segmentation intermediate_files\Haus30\associated_points\ 0 10000 0.1

REM Write some SPF files
C:\Python27_x64\python.exe write_spf_file.py input_data\Haus30.ifc 0 unassociated_cartesian discrete_param raster_grid
C:\Python27_x64\python.exe write_spf_file.py input_data\Haus30.ifc 0 unassociated_segmented discrete_param raster_grid

md generated_files\hdf

REM Convert them into HDF files
bin\ifcconvert\IfcConvert --hdf5-compress --hdf5-fix-global-id --hdf5-fix-cartesian-point    ^
    generated_files\spf\Haus30-no_grid-unassociated_cartesian-raster_grid-discrete_param.ifc ^
    generated_files\hdf\Haus30-no_grid-unassociated_cartesian-raster_grid-discrete_param.hdf ^
    --hdf5-ref-attributes IfcCartesianPointList3D.CoordList IfcDiscreteParameterValueList.Values IfcContinuousParameterValueList.Values

bin\ifcconvert\IfcConvert --hdf5-compress --hdf5-fix-global-id --hdf5-fix-cartesian-point    ^
    generated_files\spf\Haus30-no_grid-unassociated_segmented-raster_grid-discrete_param.ifc ^
    generated_files\hdf\Haus30-no_grid-unassociated_segmented-raster_grid-discrete_param.hdf ^
    --hdf5-ref-attributes IfcCartesianPointList3D.CoordList IfcDiscreteParameterValueList.Values IfcContinuousParameterValueList.Values

REM The python HDF5 serializer has been superceded by the IfcConvert executable
REM C:\Python27_x64\python.exe hdf5_serialization.py generated_files\spf\Haus30-no_grid-with_unassociated-raster_grid-discrete_param.ifc compressed single_precision

REM The viewer client has not been updated because of limitation in the h5py api in terms of variable length aggregates
REM C:\Python27_x64\python.exe client.py generated_files\spf\Haus30-no_grid-with_unassociated-raster_grid-discrete_param-compressed-single_precision.ifc
