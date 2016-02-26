/* segmentation.cpp: a component of the IFC - PointCloud schema
 * extension tooling to compress unassociated point clouds.
 *
 * Copyright (C) 2016 Eindhoven University of Technology
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD license. See the LICENSE file for details.
 */

#pragma warning(disable: 4521 4996)

#include <fstream>
#include <algorithm>

#include <boost/lexical_cast.hpp>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/segmentation/sac_segmentation.h>

typedef pcl::PointXYZ point_t;

static int points_read = 0;
static int points_written = 0;

template <typename T>
void read(char*& buffer, T& t) {
	t = *((T*) buffer);
	buffer += sizeof(T);
}

template <typename T>
void write(char*& buffer, T& t) {
	*((T*) buffer) = t;
	buffer += sizeof(T);
}

template <int N, typename T>
void writen(char*& buffer, T& t) {
	for (int i = 0; i < N; ++i) {
		write(buffer, t[i]);
	}
}

std::string make_input_filename(const std::string& dir, int n) {
	return dir + "points-subset-" + boost::lexical_cast<std::string>(n) + ".bin";
}

std::string make_output_filename(const std::string& dir, const std::string& n) {
	return dir + "segmented\\" + n + ".bin";
}

std::string make_output_filename(const std::string& dir, int n) {
	return make_output_filename(dir, "plane-" + boost::lexical_cast<std::string>(n));
}

void plane_to_axis_placement(const std::vector<float>& coef, Eigen::Vector3f& pos, Eigen::Vector3f& x, Eigen::Vector3f& y, Eigen::Vector3f& z) {
	float abc[3];
	for (int i = 0; i < 3; ++i) {
		abc[i] = fabs(coef[i]);
	}
	
	// Plane construction code below taken from:
	// https://github.com/tpaviot/oce/blob/master/src/gp/gp_Pln.cxx#L57

	if (abc[1] <= abc[0] && abc[1] <= abc[2]) {
		if (abc[0] > abc[2]) {
			pos = Eigen::Vector3f(-coef[3]/coef[0],  0.,  0.);
			z = Eigen::Vector3f(coef[0],coef[1],coef[2]);
			x = Eigen::Vector3f(-coef[2],0., coef[0]);
		} else {
			pos = Eigen::Vector3f(  0.,  0.,-coef[3]/coef[2]);
			z = Eigen::Vector3f(coef[0],coef[1],coef[2]);
			x = Eigen::Vector3f( coef[2],0.,-coef[0]);
		}
	}
	else if (abc[0] <= abc[1] && abc[0] <= abc[2]) {
		if (abc[1] > abc[2]) {
			pos = Eigen::Vector3f(  0.,-coef[3]/coef[1],  0.);
			z = Eigen::Vector3f(coef[0],coef[1],coef[2]);
			x = Eigen::Vector3f(0.,-coef[2], coef[1]);
		} else {
			pos = Eigen::Vector3f(  0.,  0.,-coef[3]/coef[2]);
			z = Eigen::Vector3f(coef[0],coef[1],coef[2]);
			x = Eigen::Vector3f(0., coef[2],-coef[1]);
		}
	}
	else {
		if (abc[0] > abc[1]) {
			pos = Eigen::Vector3f(-coef[3]/coef[0],  0.,  0.);
			z = Eigen::Vector3f(coef[0],coef[1],coef[2]);
			x = Eigen::Vector3f(-coef[1], coef[0], 0.);
		} else {
			pos = Eigen::Vector3f(  0.,-coef[3]/coef[1],  0.);
			z = Eigen::Vector3f(coef[0],coef[1],coef[2]);
			x = Eigen::Vector3f( coef[1],-coef[0], 0.);
		}
	}

	z.normalize();
	x.normalize();

	y = x.cross(z).normalized();
}

void output(const std::string& dir, int n, const std::vector<float>& coef, pcl::PointCloud<point_t>& cloud) {
	Eigen::Vector3f pos, x, y, z;
	plane_to_axis_placement(coef, pos, x, y, z);	

	/*
	// Asserts plane coefficients remained the same. These can be negated however.

	Eigen::Hyperplane<float, 3> pln = Eigen::Hyperplane<float, 3>::Through(pos+x, pos, pos+y);
	float* coef2 = pln.coeffs().data();
	for (int i = 0; i < 4; ++i) {
		if (fabs(coef[i] - coef2[i]) > 1.e-5) {
			throw std::exception();
		}
	}
	*/
	
	const size_t data_size = cloud.size() * 3 * sizeof(float) + sizeof(float) * 12;
	char* ptr = new char[data_size];
	const char* buffer = ptr;

	float* vector_data[4] = {x.data(), y.data(), z.data(), pos.data()};
	
	// Write axis 3 placement to the file
	for (int i = 0; i < 4; ++i) {
		writen<3>(ptr, vector_data[i]);
	}

	/*
	// Construct a 4x4 matrix for multiplication
	Eigen::Matrix4f m = Eigen::Matrix4f::Identity();
	float* f = m.data();	
	for (int i = 0; i < 4; ++i) {
		memcpy(&f[4*i], vector_data[i], sizeof(float) * 3);
	}
	Eigen::Matrix4f mi = m.inverse();
	*/

	for (auto it = cloud.begin(); it != cloud.end(); ++it) {
		Eigen::Vector3f p(it->data);
		Eigen::Vector3f op = p - pos;
		float uvw[3] = { op.dot(x), op.dot(y), op.dot(z) };
		writen<3>(ptr, uvw);

		points_written += 1;

		/*
		Eigen::Vector4f pp;
		memcpy(pp.data(), p.data(), sizeof(float) * 3);
		pp.data()[3] = 1.;
		Eigen::Vector4f ppp = mi * pp;
		writen<3>(ptr, ppp.data());		
		*/
	}

	std::string fn = make_output_filename(dir, n);
	std::ofstream fs(fn.c_str(), std::ios_base::binary);
	fs.write(buffer, data_size);

	delete[] buffer;
}

bool load_data(pcl::PointCloud<point_t>& cloud, const std::string& filename) {
	std::cout << filename << std::endl;

	std::ifstream fs(filename.c_str(), std::ios_base::binary);
	
	if (!fs.is_open()) {
		return false;
	}

	fs.seekg(0, std::ios::end); 
	size_t s = fs.tellg();
	fs.seekg(0, std::ios::beg);
	char* buffer = new char[s];
	fs.read(buffer, s);
	
	double x,y,z;
	uint64_t has_assoc;

	char* start = buffer;
	while (buffer < start + s) {
		read(buffer, x);
		read(buffer, y);
		read(buffer, z);
		read(buffer, has_assoc);
		
		if (has_assoc) {
			// Advance to skip association data
			buffer += sizeof(uint64_t) * 2 + sizeof(double) * 3;
		} else {
			// Only interested in un-associated points
			cloud.push_back(point_t((float)x, (float)y, (float)z));
		}
	}

	delete[] start;

	return true;
}

static char* methods[] = {"RANSAC", "LMEDS", "MSAC", "RRANSAC", "RMSAC", "MLESAC", "PROSAC"};
static int num_methods = (sizeof(methods) / sizeof(methods[0]));

int main (int argc, char** argv) {
	bool valid_command_line = argc == 5;
	
	std::string dir_name;
	int method, iterations;
	double distance;
	
	if (valid_command_line) {
		try {
			dir_name = argv[1];
			method = boost::lexical_cast<int>(argv[2]);
			iterations = boost::lexical_cast<int>(argv[3]);
			distance = boost::lexical_cast<double>(argv[4]);
		} catch (...) {
			valid_command_line = false;
		}
	}

	if (valid_command_line) {
		valid_command_line = method >= 0 && method < num_methods;
	}
	
	if (!valid_command_line) {
		std::cerr << "Usage:" << std::endl << argv[0] << " <dir> <method> <iterations> <distance>" << std::endl;
		for (int i = 0; i < num_methods; ++i) {
			std::cerr << "    method #" << i << ": " << methods[i] << std::endl;
		}
		return 1;
	}

	if (dir_name.size() && *dir_name.rbegin() != '\\') {
		dir_name += "\\";
	}

	pcl::PointCloud<point_t>::Ptr cloud, plane, plane_cluster, remaining;

	cloud.reset(new pcl::PointCloud<point_t>);
	plane.reset(new pcl::PointCloud<point_t>);
	remaining.reset(new pcl::PointCloud<point_t>);
	plane_cluster.reset(new pcl::PointCloud<point_t>);

	for(int i = 0; i < 100; ++i) {
		const size_t old_size = cloud->size();
		if (load_data(*cloud, make_input_filename(dir_name, i))) {
			std::cerr << "Loaded " << (cloud->size() - old_size) << " points" << std::endl;
		}
	}
	
	pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients());
	pcl::PointIndices::Ptr inliers(new pcl::PointIndices());
	pcl::SACSegmentation<pcl::PointXYZ> seg;
	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_PLANE);

	std::cerr << "Using " << iterations << " iterations of " << methods[method] << " with a distance of " << distance << std::endl;

	seg.setMethodType(method);
	seg.setMaxIterations(iterations);
	seg.setDistanceThreshold(distance);
	
	pcl::ExtractIndices<pcl::PointXYZ> extract;
	pcl::ExtractIndices<pcl::PointXYZ> extract2;

	/*
	// Debug: Write segmented points in single contiguous file in original coordinates
	const std::string fn2 = make_output_filename(dir_name, "written-associated-points");
	std::ofstream fs2(fn2.c_str(), std::ios_base::binary);
	*/

	size_t i = 0, nr_points = cloud->points.size();
	
	while (cloud->points.size() > (nr_points / 4))
	{
		// Segment the largest planar component from the remaining cloud
		seg.setInputCloud (cloud);
		seg.segment (*inliers, *coefficients);
		if (inliers->indices.size() == 0)
		{
			std::cerr << "Could not estimate a planar model for the given dataset." << std::endl;
			break;
		}

		pcl::IndicesPtr all_indices(new std::vector<int>);
		extract.setIndices (inliers);

		// Extract the inliers
		extract.setInputCloud (cloud);		
		extract.setNegative (false);
		extract.filter (*plane);

		size_t ref_index;

		std::vector<int> labeled(plane->size(), 0);
		pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
		kdtree.setInputCloud (plane);
		std::map<int, int> cluster_sizes;
		int max_cluster_size = 0;
		std::vector< pcl::IndicesPtr > clusters;

		for (auto it = plane->begin(); it != plane->end(); ++it) {
			auto i = std::distance(plane->begin(), it);
			if (!labeled[i]) {
				std::vector<int>* indices = new std::vector<int>;
				std::vector<int>* indices2 = new std::vector<int>;
				std::vector<float> distances;
				size_t indices_offset = 0;
				indices->push_back(i);
				while (true) {
					std::vector<int> indices_copy;
					for (auto jt = indices->begin() + indices_offset;  jt != indices->end(); ++jt) {
						std::vector<int> tmp; // radiusSearch() clears vector (why??)
						kdtree.radiusSearch(*jt, distance * 2, tmp, distances);
						std::copy(tmp.begin(), tmp.end(), std::back_inserter(indices_copy));
					}
					std::set<int> indices_set(indices->begin(), indices->end());
					std::set<int> new_indices(indices_copy.begin(), indices_copy.end());
					indices_offset = indices->size();
					std::set_difference(new_indices.begin(), new_indices.end(), indices_set.begin(), indices_set.end(), std::back_inserter(*indices));
					if (indices_offset == indices->size()) {
						break;
					}
				}
				for (auto jt = indices->begin();  jt != indices->end(); ++jt) {
					auto j = *jt;
					if (!labeled[j]) {
						labeled[j] = i + 1;
						indices2->push_back(j);
					}
				}
				delete indices;
				cluster_sizes[i + 1] = indices2->size();
				max_cluster_size = (std::max)(max_cluster_size, (int) indices2->size());
				clusters.push_back(pcl::IndicesPtr(indices2));				
			}
		}

		std::cerr << "Found plane representing: " << plane->size() << " data points." << std::endl;
		std::cerr << "Clusters: ";
		std::vector<int>* indices_within_clusters = new std::vector<int>;

		for (auto it = clusters.begin(); it != clusters.end(); ++it) {
			if ((**it).size() > max_cluster_size / 4) {
				std::cerr << (**it).size() << " ";
				extract2.setInputCloud (plane);
				extract2.setIndices (*it);
				extract2.setNegative (false);
				extract2.filter (*plane_cluster);

				output(dir_name, i++, coefficients->values, *plane_cluster);

				/*
				// Debug: Write segmented points in single contiguous file in original coordinates

				const size_t data_size = plane_cluster->size() * sizeof(float) * 3;
				char* ptr = new char[data_size];
				const char* buffer = ptr;

				for (auto pit = plane_cluster->begin(); pit != plane_cluster->end(); ++pit) {
					writen<3>(ptr, pit->data);
				}
				
				fs2.write(buffer, data_size);
				delete[] buffer;
				*/

				indices_within_clusters->insert(indices_within_clusters->end(), (**it).begin(), (**it).end());
			}
		}

		std::cerr << std::endl;

		size_t original_inliers_size;
		const std::vector<int>& original_inliers = inliers->indices;
		original_inliers_size = original_inliers.size();
		
		std::cerr << "Inliers shrunk from " << original_inliers_size << " to " << (indices_within_clusters->size()) << std::endl;

		// Swap point cloud with the outliers
		extract.setIndices(pcl::IndicesPtr(indices_within_clusters));
		extract.setNegative(true);
		extract.filter(*remaining);

		std::cerr << (100 * remaining->size() / nr_points) << "% processed." << std::endl;
		
		cloud.swap(remaining);
	}
	
	const size_t data_size = cloud->size() * sizeof(float) * 3;
	char* ptr = new char[data_size];
	const char* buffer = ptr;

	for (auto it = cloud->begin(); it != cloud->end(); ++it) {
		writen<3>(ptr, it->data);
	}

	const std::string fn = make_output_filename(dir_name, "remaining-unassociated-points");
	std::ofstream fs(fn.c_str(), std::ios_base::binary);
	fs.write(buffer, data_size);

	delete[] buffer;

	std::cout << "Written " << points_written << " / " << nr_points;
	std::cin.get();
}