# IFCPointCloud prototype implementation

A prototype tool for point cloud and IFC association and serialization written in Python. The input for this tools contains of three parts:

* An IFC file
* A set of point cloud scans as gzipped PCD-ascii files
* An alignment matrix that transforms points from the point cloud scan into the coordinate system of the IFC model

One of such datasets is provided in this repository.

The prototype tool itself consists of four executables. An example how to invoke these is provided below.

1. Point cloud association (point_association.py)
    
    Input:
    * [**IFC file**] SPF-serialized IFC instance model that describes the building element
    data.
    * [**Point cloud scans**] A filename pattern that points to a set of gzip-compressed
    ASCII-encoded .pcd files.
    * [**Alignment matrix**] An alignment matrix that transforms points from the point cloud scan into the coordinate system of the IFC model in JSON format.

    Output:
    * [**Cached BRep representations**] For every IfcProduct descendant in the
    IFC file, an Open Cascade BRep file is generated that represents the
    "Body" representation of the product in global model coordinates.
    * [**Associated points**] A binary sequence of the Cartesian coordinates (as transformed
    by the alignment matrix) for every point in the .pcd input files.
    If the point is close to an IFC building element surface, the Cartesian
    coordinates are followed by a) a reference to the product b) an index that
    points to a face in the Cached BRep representations c) the parametric
    coordinates and deviation of the point as projected onto the surface. The
    points are subdivided over several passes so that a contiguous subset of
    the points yields a randomly subsampled subset.

2. IFC SPF generation (write_spf_file.py)
    
    Input:
    * [**IFC file**] SPF-serialized IFC instance model that describes the building element
    data.
    * [**Cached BRep representations**] A reference to the geometry cache generated
    earlier.
    * [**Associated points**] The sequence of associated points as generated above.

    Output:
    * [**IFC file**] An SPF-serialized IFC instance model that, in addition to the original
    data, contains point cloud data and association relations conforming to
    the schema extension distributed as part of this repository.

3. HDF5 serialization (hdf5_serialization.py)

    Input:
    * [**IFC file**] SPF-serialized IFC instance model that describes the building element
    data.
    * [**IFC schema**] The EXPRESS (ISO 10303-11) rendition of the complete IFC
    schema.

    Output:
    * [**IFC file**] An HDF5-serialized version of the input file, generated according to
    ISO 10303-26.

4. IFC and point cloud viewer (viewer.py)

    Input:
    * [**IFC file**] An HDF-serialized IFC file with or without point cloud
    data. Configuration of the viewer, i.e what to visualize and in what level
    of detail happens in-code.

Due to some of the binary dependencies, currently the tool only works on Windows.

Example on how to invoke the executables:

    point_association.py input_data\Haus30.ifc input_data\Haus30\%d.pcd.gz input_data\Haus30_alignment_matrix.json --use-display
    write_spf_file.py input_data\Haus30.ifc 0 with_unassociated discrete_param raster_grid
    hdf5_serialization.py generated_files\spf\Haus30-no_grid-with_unassociated-raster_grid-discrete_param.ifc compressed single_precision
    viewer.py generated_files\spf\Haus30-no_grid-with_unassociated-raster_grid-discrete_param-compressed-single_precision.ifc
