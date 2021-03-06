# IFCPointCloud

In this repository an extension to the Industry Foundation Classes schema is introduced with the capabilities to efficiently describe and store point cloud data. The repository contains schemas, prototype tools and documentation.

The prototype tool is capable of associating external point clouds and describing them parametrically according to IFC building surfaces. As such, a unified format of IFC building data is obtained that represents point cloud data more efficiently than prevalent point cloud formats. The decomposition and association between the two data sources enriches the semantics of both. The point clouds are labelled with the building elements they describe and can be navigated according to the spatial subdivision structure of IFC. IFC building elements obtain an additional detailed representation that describes surface characteristics and the exact as-built physical form.

A binary serialization is introduced that writes IFC files (with or without point clouds) in an efficient hierarchical structure for efficient partial retrieval. This binary serialization format is based on ISO-13030-26 and uses the HDF5 file format.
