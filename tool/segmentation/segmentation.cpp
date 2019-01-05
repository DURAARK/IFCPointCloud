/* segmentation.cpp: a component of the IFC - PointCloud schema
 * extension tooling to compress unassociated point clouds.
 *
 * Copyright (C) 2016 Eindhoven University of Technology
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD license. See the LICENSE file for details.
 */

#define BOOST_OPTIONAL_USE_OLD_DEFINITION_OF_NONE 1

#include <Eigen/LU>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <boost/container/flat_set.hpp>

#include <ifcgeom/IfcGeom.h>
#include <ifcgeom/IfcGeomIterator.h>

#include <TopExp_Explorer.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <BRepTools_ShapeSet.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Plane.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <ProjLib.hxx>

#include <pcl/point_types.h>
#include <pcl/common/common.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/segmentation/sac_segmentation.h>

#include <iostream>
#include <fstream>
#include <algorithm>

typedef pcl::PointXYZ point_t;


float Ceil(float x) {
	return std::ceil(x);
}

float Floorq(float x) {
	return std::floor(x);
}


class voxel_mapping {
public:

	class voxel {
	public:
		Eigen::Vector3f centroid(const pcl::PointCloud<point_t>::Ptr& cloud) {
			Eigen::Vector3f v = Eigen::Vector3f::Zero();
			for (auto it = indices.begin(); it != indices.end(); ++it) {
				v += Eigen::Vector3f((*cloud)[*it].data);
			}
			return v / indices.size();
		}
		std::vector<int> indices;
		bool visited;
		Eigen::Array3i index;
		size_t n;
	};

	template <int N>
	class sector {
	public:
		sector(voxel_mapping& mapping, const Eigen::Array3i& index)
			: mapping_(mapping)
			, index_(index)
		{
			voxels_.resize(N*N*N);
		}

		/*
		voxel& get(size_t n) {
			if (voxels_[n] == nullptr) {
				mapping_.voxels_flat_.push_back(voxels_[n] = new voxel);
			}
			return *voxels_[n];
		}
		*/

		voxel& get(const Eigen::Array3i& idx) {
			size_t n = idx[0] + idx[1] * N + idx[2] * N * N;
			if (voxels_[n] == nullptr) {
				size_t m = mapping_.voxels_flat_.size();
				mapping_.voxels_flat_.push_back(voxels_[n] = new voxel);
				voxels_[n]->index = index_ * N + idx;
				voxels_[n]->n = m;
			}
			return *voxels_[n];
		}

		voxel* get_if_exists(const Eigen::Array3i& idx) {
			size_t n = idx[0] + idx[1] * N + idx[2] * N * N;
			return voxels_[n];
		}

		voxel_mapping& mapping_;
		// vector for zero initialization
		std::vector<voxel*> voxels_;
		const Eigen::Array3i& index_;
	};

	voxel_mapping(const pcl::PointCloud<point_t>::Ptr& cloud, float size)
		: cloud_(cloud)
	{
		Eigen::Array3f pmin = Eigen::Array3f::Constant(std::numeric_limits<float>::infinity());
		Eigen::Array3f pmax = Eigen::Array3f::Constant(-std::numeric_limits<float>::infinity());
		for (auto it = cloud_->begin(); it != cloud_->end(); ++it) {
			auto p = it->getArray3fMap();
			for (int i = 0; i < 3; ++i) {
				if (p[i] < pmin[i]) {
					pmin[i] = p[i];
				}
				if (p[i] > pmax[i]) {
					pmax[i] = p[i];
				}
 			}
		}
		// To make it into an inclusive boundary value
		pmax += Eigen::Array3f::Constant(1e-3);
		origin_ = pmin;
		size_ = Eigen::Array3f::Constant(size);
		extent_ = ((pmax - pmin) / (size_ * 8)).unaryExpr(&Ceil).cast<int>();
		sectors_.resize(extent_[0] * extent_[1] * extent_[2]);
		for (size_t i = 0; i < cloud_->size(); ++i) {
			voxel& v = get((*cloud_)[i]);
			v.indices.push_back(i);
			// get_if_exists(v.index);
		}
	}

	voxel& get(const point_t& p) {
		Eigen::Array3i idx = ((Eigen::Array3f(p.data) - origin_) / size_).unaryExpr(&Floorq).cast<int>();
		return get(idx);
	}

	voxel& get(size_t n) {
		return *voxels_flat_[n];
	}
	
	voxel& get(const Eigen::Array3i& idx) {
		Eigen::Array3i cluster = idx / 8;
		size_t n = cluster[0] + cluster[1] * extent_[0] + cluster[2] * extent_[0] * extent_[1];
		if (sectors_[n] == nullptr) {
			sectors_[n] = new sector<8>(*this, cluster);
		}
		return sectors_[n]->get(idx - cluster * 8);
	}

	voxel* get_if_exists(const Eigen::Array3i& idx) {
		Eigen::Array3i cluster = idx / 8;
		if (idx[0] < 0 || cluster[0] >= extent_[0] ||
			idx[1] < 0 || cluster[1] >= extent_[1] ||
			idx[2] < 0 || cluster[2] >= extent_[2])
		{
			return nullptr;
		}
		size_t n = cluster[0] + cluster[1] * extent_[0] + cluster[2] * extent_[0] * extent_[1];
		if (sectors_[n] == nullptr) {
			return nullptr;
		}
		return sectors_[n]->get_if_exists(idx - cluster * 8);
	}

	pcl::PointCloud<point_t>::Ptr applyFilter() {
		pcl::PointCloud<point_t>::Ptr cloud(new pcl::PointCloud<point_t>);
		cloud->resize(voxels_flat_.size());
		for (auto it = voxels_flat_.begin(); it != voxels_flat_.end(); ++it) {
			Eigen::Vector3f p = (*it)->centroid(cloud_);
			memcpy((*cloud)[std::distance(voxels_flat_.begin(), it)].data, p.data(), sizeof(float) * 3);
		}
		return cloud;
	}

	const pcl::PointCloud<point_t>::Ptr& cloud() const { return cloud_; }

	std::vector<std::vector<int>> components() {
		for (auto it = voxels_flat_.begin(); it != voxels_flat_.end(); ++it) {
			(**it).visited = false;
		}
		std::vector<std::vector<int>> return_value;
		for (auto it = voxels_flat_.begin(); it != voxels_flat_.end(); ++it) {
			if (!(**it).visited) {
				std::vector<int> current;
				std::queue<voxel*> queue;
				queue.push(*it);
				(**it).visited = true;
				current.push_back((**it).n);
				while (!queue.empty()) {
					voxel* c = queue.front();
					queue.pop();
					for (int i = -1; i < 2; ++i) {
						for (int j = -1; j < 2; ++j) {
							for (int k = -1; k < 2; ++k) {
								if (i == 0 && j == 0 && k == 0) {
									continue;
								}
								Eigen::Array3i d(i, j, k);
								voxel* v = get_if_exists(c->index + d);
								if (v == nullptr || v->visited) {
									continue;
								}
								v->visited = true;
								queue.push(v);
								current.push_back(v->n);
							}
						}
					}
				}
				return_value.push_back(current);
			}
		}
		return return_value;
	}

	const pcl::PointCloud<point_t>::Ptr cloud_;
	Eigen::Array3f origin_;
	Eigen::Array3f size_;
	Eigen::Array3i extent_;
	std::vector<sector<8>*> sectors_;
	std::vector<voxel*> voxels_flat_;
};

class voxel_subset {

public:

	voxel_subset(voxel_mapping& mapping, const std::vector<int>& indices)
		: mapping_(mapping)
		, indices_(indices)
	{}

	class subset_iterator {
	public:
		subset_iterator(voxel_subset& voxels, const std::vector<int>& indices, bool begin)
			: voxels_(voxels)
			, indices_(indices)
		{
			if (begin) {
				it1_ = indices.begin();
				it2_ = voxels_.mapping_.get(*it1_).indices.begin();
			} else {
				it1_ = indices.end();
			}
		}

		subset_iterator& operator++() {
			it2_++;
			if (it2_ == voxels_.mapping_.get(*it1_).indices.end()) {
				it1_++;
				if (it1_ != indices_.end()) {
					it2_ = voxels_.mapping_.get(*it1_).indices.begin();
				}
			}
			return *this;
		}

		/*const point_t& operator*() const {
			return (*voxels_.mapping_.cloud())[*it2_];
		}*/

		const int& operator*() const {
			return *it2_;
		}

		bool operator==(const subset_iterator& other) const {
			if (it1_ == indices_.end()) {
				return other.it1_ == indices_.end();
			} else {
				return it1_ == other.it1_ && it2_ == other.it2_;
			}
		}

		bool operator!=(const subset_iterator& other) const {
			return !(*this == other);
		}

		voxel_subset& voxels_;
		const std::vector<int>& indices_;
		std::vector<int>::const_iterator it1_;
		std::vector<int>::const_iterator it2_;
	};

	typedef subset_iterator const_iterator;

	const_iterator begin() {
		return const_iterator(*this, indices_, true);
	}

	const_iterator end() {
		return const_iterator(*this, indices_, false);
	}

	size_t size() {
		size_t s = 0;
		for (auto it = indices_.begin(); it != indices_.end(); ++it) {
			s += mapping_.get(*it).indices.size();
		}
		return s;
	}

	voxel_mapping& mapping_;
	const std::vector<int>& indices_;
};

typedef std::pair<pcl::PointCloud<point_t>::Ptr, pcl::PointCloud<point_t>::Ptr> pointcloud_pair_t;

/* std::vector<int> */
template <typename T>
pointcloud_pair_t split_pointcloud(const pcl::PointCloud<point_t>::Ptr& cloud, const T& indices, bool only_in = false) {

	std::pair<pcl::PointCloud<point_t>::Ptr, pcl::PointCloud<point_t>::Ptr> ret;
	ret = std::make_pair(
		pcl::PointCloud<point_t>::Ptr(new pcl::PointCloud<point_t>),
		pcl::PointCloud<point_t>::Ptr(new pcl::PointCloud<point_t>));

	ret.first->reserve(indices.size());
	ret.second->reserve(cloud->size() - indices.size());
	std::vector<int> sorted_indices(indices.begin(), indices.end());

	std::sort(sorted_indices.begin(), sorted_indices.end());

	auto cloud_it = cloud->begin();
	auto idx_it = sorted_indices.begin();

	for (; cloud_it < cloud->end(); ++cloud_it) {
		if (idx_it != sorted_indices.end() && std::distance(cloud->begin(), cloud_it) == *idx_it) {
			ret.first->push_back(*cloud_it);
			idx_it++;
		} else if (!only_in) {
			ret.second->push_back(*cloud_it);
		}
	}

	// Make sure all indices have been exhausted
	if (idx_it != sorted_indices.end()) {
		std::cerr << std::distance(idx_it, sorted_indices.end());
		throw std::exception();
	}

	if (!only_in) {
		if (ret.first->size() + ret.second->size() != cloud->size()) {
			throw std::exception();
		}
	}

	return ret;
}

class lazy_pointcloud {
public:
	typedef boost::container::flat_set<int> indices_type;

	template <typename V>
	class point_iterator {
	public:
		point_iterator(const lazy_pointcloud& cloud, const indices_type::const_iterator& it)
			: cloud_(cloud)
			, it1_(it)
		{
			indexed_ = cloud_.indexed();
			if (!indexed_) {
				throw std::exception("Indexed iteration over an unindexed cloud");
			}
		}

		point_iterator(const lazy_pointcloud& cloud, const pcl::PointCloud<point_t>::const_iterator& it)
			: cloud_(cloud)
			, it2_(it)
		{
			indexed_ = cloud_.indexed();
		}

		point_iterator& operator++() {
			if (indexed_) {
				++ it1_;
			} else {
				++ it2_;
			}
			return *this;
		}

		V& operator*() const {
			if (indexed_) {
				return (*cloud_.cloud())[*it1_];
			} else {
				return *it2_;
			}
		}

		size_t original_index() const {
			if (indexed_) {
				// if (cloud_.indices().empty()) {
				if (it1_ == cloud_.indices().end()) {
					return cloud_.cloud()->size();
				} else {
					return *it1_;
				}
			} else {
				return it2_ - cloud_.cloud()->begin();
			}
		}

		bool operator==(const point_iterator& other) const {
			return original_index() == other.original_index();
		}

		bool operator!=(const point_iterator& other) const {
			return !(*this == other);
		}

	private:
		const lazy_pointcloud& cloud_;
		indices_type::const_iterator it1_;
		pcl::PointCloud<point_t>::const_iterator it2_;
		bool indexed_;
	};

	typedef point_iterator<const point_t> const_iterator;

	lazy_pointcloud(const pcl::PointCloud<point_t>::Ptr& cloud)
		: cloud_(cloud)
		, tree_(nullptr)
		, indexed_(false)
	{
		if (!cloud->empty()) {
			std::cerr << "Building tree" << std::endl;
			tree_.reset(new pcl::search::KdTree<point_t>);
			tree_->setInputCloud(cloud);
		}
	}

	lazy_pointcloud(const pcl::PointCloud<point_t>::Ptr& cloud, const pcl::search::KdTree<point_t>::Ptr& tree)
		: cloud_(cloud)
		, tree_(tree)
		, indexed_(false)
	{}

	lazy_pointcloud(const pcl::PointCloud<point_t>::Ptr& cloud, const indices_type& indices)
		: cloud_(cloud)
		, tree_(nullptr)
		, indices_(indices)
		, indexed_(true)
	{
		if (!cloud->empty()) {
			std::cerr << "Building tree" << std::endl;
			tree_.reset(new pcl::search::KdTree<point_t>);
			tree_->setInputCloud(cloud);
		}
	}

	lazy_pointcloud(const pcl::PointCloud<point_t>::Ptr& cloud, const pcl::search::KdTree<point_t>::Ptr& tree, const indices_type& indices)
		: cloud_(cloud)
		, tree_(tree)
		, indices_(indices)
		, indexed_(true)
	{}

	const_iterator begin() const {
		if (indexed_) {
			return const_iterator(*this, indices().begin());
		} else {
			return const_iterator(*this, cloud()->begin());
		}
	}

	const_iterator end() const {
		if (indexed_) {
			return const_iterator(*this, indices().end());
		} else {
			return const_iterator(*this, cloud()->end());
		}
	}

	// friend lazy_pointcloud operator+(const lazy_pointcloud& a, const lazy_pointcloud& b);
	// friend lazy_pointcloud operator-(const lazy_pointcloud& a, const lazy_pointcloud& b);

	pcl::PointCloud<point_t>::Ptr extract() const {
		if (!indexed_) {
			return cloud();
		} else {
			auto x = split_pointcloud(cloud_, indices_, true);
			return x.first;
		}
	}

	/*
	std::pair<lazy_pointcloud, lazy_pointcloud> split() const {
		if (!indexed_) {
			pcl::PointCloud<point_t>::Ptr cloud2(new pcl::PointCloud<point_t>);
			return std::make_pair(lazy_pointcloud(cloud(), tree_), lazy_pointcloud(cloud2));
		} else {
			auto x = split_pointcloud(cloud_, indices_);
			return std::make_pair(lazy_pointcloud(x.first), lazy_pointcloud(x.second));
		}
	}
	*/

	template <typename T>
	std::pair<lazy_pointcloud, lazy_pointcloud> split(const T& idxs) const {
		std::vector<int> sorted(idxs->begin(), idxs->end());
		std::sort(sorted.begin(), sorted.end());
		indices_type indices;
		indices.insert(boost::container::ordered_unique_range, sorted.begin(), sorted.end());

		lazy_pointcloud other(cloud(), tree_, indices);
		lazy_pointcloud x = intersect(other);
		lazy_pointcloud y = subtract(x);
		return std::make_pair(x, y);
	}

	const pcl::PointCloud<point_t>::Ptr& cloud() const {
		return cloud_;
	}

	const indices_type& indices() const {
		return indices_;
	}

	bool indexed() const {
		return indexed_;
	}

	lazy_pointcloud radius_search(const point_t& center, float radius) const {
		
		std::vector<int> indices;
		std::vector<float> distances;
		tree_->radiusSearch(center, radius, indices, distances);
		std::sort(indices.begin(), indices.end());

		/*{
			float mi = std::numeric_limits<float>::infinity();
			float ma = -std::numeric_limits<float>::infinity();

			for (auto it = indices.begin(); it != indices.end(); ++it) {
				const float* xyz1 = (*cloud_)[*it].data;
				const float* xyz2 = center.data;

				float accum = 0.f;
				for (int i = 0; i < 3; ++i) {
					const float d = xyz1[i] - xyz2[i];
					accum += d * d;
				}

				accum = sqrt(accum);
				if (accum < mi) {
					mi = accum;
				}
				if (accum > ma) {
					ma = accum;
				}
			}
			std::cerr << "Range [" << mi << "," << ma << ")";
		}*/

		indices_type intersection;
		if (indexed()) {
			std::vector<int> v_intersection((std::max)(indices.size(), indices_.size()));
			auto it = std::set_intersection(indices.begin(), indices.end(), indices_.begin(), indices_.end(), v_intersection.begin());
			v_intersection.resize(it - v_intersection.begin());
			intersection.insert(boost::container::ordered_unique_range, v_intersection.begin(), v_intersection.end());
		} else {			
			intersection.insert(boost::container::ordered_unique_range, indices.begin(), indices.end());
		}
		
		lazy_pointcloud result(cloud(), tree_, intersection);
		/*{
			float mi = std::numeric_limits<float>::infinity();
			float ma = -std::numeric_limits<float>::infinity();
			
			for (auto it = result.begin(); it != result.end(); ++it) {
				const float* xyz1 = (*it).data;
				const float* xyz2 = center.data;

				float accum = 0.f;
				for (int i = 0; i < 3; ++i) {
					const float d = xyz1[i] - xyz2[i];
					accum += d * d;
				}

				accum = sqrt(accum);
				if (accum < mi) {
					mi = accum;
				}
				if (accum > ma) {
					ma = accum;
				}
			}
			std::cerr << " [" << mi << "," << ma << ")" << std::endl;
		}*/

		return result;
	}

	size_t size() const {
		return indexed_ ? indices_.size() : cloud_->size();
	}

	lazy_pointcloud intersect(const lazy_pointcloud& other) const {
		if (cloud() != other.cloud()) {
			throw std::exception("Lazy clouds can only be intersected if they share the same source");
		}

		if (!indexed()) {
			return other;
		}
		if (!other.indexed()) {
			return *this;
		}

		std::vector<int> x((std::min)(size(), other.size()));
		auto it = std::set_intersection(indices().begin(), indices().end(), other.indices().begin(), other.indices().end(), x.begin());
		x.resize(it - x.begin());

		// std::sort() necessary?

		indices_type indices;
		indices.insert(boost::container::ordered_unique_range, x.begin(), x.end());

		return lazy_pointcloud(cloud(), tree_, indices);
	}

	lazy_pointcloud subtract(const lazy_pointcloud& other) const {
		if (cloud() != other.cloud()) {
			throw std::exception("Lazy clouds can only be intersected if they share the same source");
		}

		if (!other.indexed()) {
			return create_empty();
		}

		if (!indexed()) {
			std::vector<int> idxs;
			idxs.reserve(size() - other.size());
			for (int i = 0; i < size(); ++i) {
				if (other.indices().find(i) == other.indices().end()) {
					idxs.push_back(i);
				}
			}
			// sorted by construction
			indices_type indices;
			indices.insert(boost::container::ordered_unique_range, idxs.begin(), idxs.end());
			return lazy_pointcloud(cloud(), tree_, indices);
		}

		std::vector<int> x(indices().size() - other.indices().size());
		auto it = std::set_difference(indices().begin(), indices().end(), other.indices().begin(), other.indices().end(), x.begin());
		x.resize(it - x.begin());

		// std::sort() necessary?

		indices_type indices;
		indices.insert(boost::container::ordered_unique_range, x.begin(), x.end());

		return lazy_pointcloud(cloud(), tree_, indices);
	}

	static lazy_pointcloud create_empty() {
		pcl::PointCloud<point_t>::Ptr cloud(new pcl::PointCloud<point_t>);
		return lazy_pointcloud(cloud);
	}

private:
	pcl::PointCloud<point_t>::Ptr cloud_;
	pcl::search::KdTree<point_t>::Ptr tree_;
	indices_type indices_;
	bool indexed_ = false;
};

/*
lazy_pointcloud operator+(const lazy_pointcloud& a, const lazy_pointcloud& b) {
	if (a.cloud() == b.cloud()) {
		if (a.indexed_ && b.indexed_) {
			lazy_pointcloud::indices_type indices;
			indices.insert(a.indices_.begin(), a.indices_.end());
			indices.insert(b.indices_.begin(), b.indices_.end());
			return lazy_pointcloud(a.cloud(), a.tree_, indices);
		} else {
			return lazy_pointcloud(a.cloud(), a.tree_);
		}
	} else {
		lazy_pointcloud A = a.flatten();
		lazy_pointcloud B = b.flatten();
		pcl::PointCloud<point_t>::Ptr cloud(new pcl::PointCloud<point_t>);
		*cloud += *A.cloud();
		*cloud += *B.cloud();
		return lazy_pointcloud(cloud);
	}
}

lazy_pointcloud operator-(const lazy_pointcloud& a, const lazy_pointcloud& b) {
	if (a.cloud() == b.cloud()) {
		if (a.indexed_ && b.indexed_) {
			lazy_pointcloud::indices_type indices;
			indices.insert(a.indices_.begin(), a.indices_.end());
			indices.erase(b.indices_.begin(), b.indices_.end());
			return lazy_pointcloud(a.cloud(), a.tree_, indices);
		} else if (!a.indexed_ && !b.indexed_) {
			pcl::PointCloud<point_t>::Ptr cloud(new pcl::PointCloud<point_t>);
			return lazy_pointcloud(cloud);
		} else {
			lazy_pointcloud::indices_type indices;
			for (int i = 0; i < a.cloud()->size(); ++i) {
				indices.insert(i);
			}
			if (a.indexed_) {
				indices.erase(a.indices_.begin(), a.indices_.end());
			} else {
				indices.erase(b.indices_.begin(), b.indices_.end());
			}
			return lazy_pointcloud(a.cloud(), a.tree_);
		}
	} else {
		throw std::exception("Lazy clouds can only be subtracted if they share the same source");
	}
}
*/

class stop_iteration : public std::exception {};

double face_area(const TopoDS_Face& f) {
	GProp_GProps prop;
	BRepGProp::SurfaceProperties(f, prop);
	return prop.Mass();
}

template <typename T>
T read_from_stream(std::istream& stream) {
	T t;
	stream >> t;
	if (stream.eof()) {
		throw stop_iteration();
	}
	return t;
}

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>

template <>
TopoDS_Face read_from_stream(std::istream& stream) {
	BRep_Builder B;
	BRepTools_ShapeSet ss(B);
	ss.Read(stream);

	// int i = ss.Index(ss.Shape(ss.NbShapes()));
	// ss.Locations().Dump(std::cerr);

	return TopoDS::Face(TopoDS_Iterator(ss.Shape(ss.NbShapes())).Value());

	/* while (true) {
		TopoDS_Shape S;
		ss.Read(S, stream);
		if (S.ShapeType() == TopAbs_FACE) {
			return TopoDS::Face(S);
		}
	} */
	// return TopoDS::Face(ss.Shape(ss.NbShapes())); //  .Moved(ss.Locations().Location(ss.NbShapes())));
	// BRep_Builder builder;
	// TopoDS_Face face;
	// BRepTools::Read(face, stream, builder);
	// return face;
}

TopoDS_Face string_to_face(const std::string& s) {
	std::istringstream ss(s);
	return read_from_stream<TopoDS_Face>(ss);
}

std::string face_to_string(const TopoDS_Face& f) {
	std::ostringstream sss;
	BRepTools_ShapeSet ss;
	ss.Add(f);
	ss.Write(f, sss);
	return sss.str();
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
			pos = Eigen::Vector3f(-coef[3] / coef[0], 0., 0.);
			z = Eigen::Vector3f(coef[0], coef[1], coef[2]);
			x = Eigen::Vector3f(-coef[2], 0., coef[0]);
		} else {
			pos = Eigen::Vector3f(0., 0., -coef[3] / coef[2]);
			z = Eigen::Vector3f(coef[0], coef[1], coef[2]);
			x = Eigen::Vector3f(coef[2], 0., -coef[0]);
		}
	} else if (abc[0] <= abc[1] && abc[0] <= abc[2]) {
		if (abc[1] > abc[2]) {
			pos = Eigen::Vector3f(0., -coef[3] / coef[1], 0.);
			z = Eigen::Vector3f(coef[0], coef[1], coef[2]);
			x = Eigen::Vector3f(0., -coef[2], coef[1]);
		} else {
			pos = Eigen::Vector3f(0., 0., -coef[3] / coef[2]);
			z = Eigen::Vector3f(coef[0], coef[1], coef[2]);
			x = Eigen::Vector3f(0., coef[2], -coef[1]);
		}
	} else {
		if (abc[0] > abc[1]) {
			pos = Eigen::Vector3f(-coef[3] / coef[0], 0., 0.);
			z = Eigen::Vector3f(coef[0], coef[1], coef[2]);
			x = Eigen::Vector3f(-coef[1], coef[0], 0.);
		} else {
			pos = Eigen::Vector3f(0., -coef[3] / coef[1], 0.);
			z = Eigen::Vector3f(coef[0], coef[1], coef[2]);
			x = Eigen::Vector3f(coef[1], -coef[0], 0.);
		}
	}

	z.normalize();
	x.normalize();

	y = x.cross(z).normalized();
}

void xyz_to_vec(const gp_XYZ& xyz, Eigen::Vector3f& vec) {
	vec.data()[0] = (float)xyz.X();
	vec.data()[1] = (float)xyz.Y();
	vec.data()[2] = (float)xyz.Z();
}

template <typename T>
void dump(const T& t);

template <>
void dump(const gp_XYZ& p) {
	std::cerr << "(" << p.X() << " " << p.Y() << " " << p.Z() << ")" << std::endl;
}

template <>
void dump(const gp_Pnt& p) {
	dump(p.XYZ());
}

template <>
void dump(const gp_Dir& d) {
	dump(d.XYZ());
}

template <>
void dump(const Eigen::Vector3f& v) {
	std::cerr << "(" << v.x() << " " << v.y() << " " << v.z() << ")" << std::endl;
}

template <>
void dump(const Eigen::Vector4f& v) {
	std::cerr << "(" << v.x() << " " << v.y() << " " << v.z() << " " << v.w() << ")" << std::endl;
}

template <typename T>
ptrdiff_t calc_distance(const T& a, const T& b) {
	return std::distance(a, b);
}

template <>
ptrdiff_t calc_distance(const lazy_pointcloud::const_iterator& a, const lazy_pointcloud::const_iterator& b) {
	// original_index() returns offset from complete cloud begin()
	return b.original_index();
};

class rectangular_trimmed_plane {
public:
	rectangular_trimmed_plane() {}

	rectangular_trimmed_plane(const Eigen::Vector4f& plane_coefficients)
		: plane_coefficients(plane_coefficients)
	{
		const float* data = plane_coefficients.data();
		std::vector<float> coef(data, data + 4);
		plane_to_axis_placement(coef, pos, x, y, z);

		init_bounds();
	}

	rectangular_trimmed_plane(const gp_Pln& p)
	{
		double a, b, c, d;
		p.Coefficients(a, b, c, d);
		plane_coefficients.data()[0] = a;
		plane_coefficients.data()[1] = b;
		plane_coefficients.data()[2] = c;
		plane_coefficients.data()[3] = d;

		xyz_to_vec(p.Position().Location().XYZ(), pos);
		xyz_to_vec(p.Position().XDirection().XYZ(), x);
		xyz_to_vec(p.Position().YDirection().XYZ(), y);
		xyz_to_vec(p.Position().Direction().XYZ(), z);

		init_bounds();
	}

	void fill(float u, float v) {
		uv_min.x() = (std::min)(uv_min.x(), u);
		uv_min.y() = (std::min)(uv_min.y(), v);

		uv_max.x() = (std::max)(uv_max.x(), u);
		uv_max.y() = (std::max)(uv_max.y(), v);
	}

	void fill(const TopoDS_Face& face) {
		TopExp_Explorer exp(face, TopAbs_VERTEX);
		for (; exp.More(); exp.Next()) {
			const TopoDS_Vertex& V = TopoDS::Vertex(exp.Current());
			gp_Pnt pnt = BRep_Tool::Pnt(V);
			float xyz[] = { (float)pnt.X(), (float)pnt.Y(), (float)pnt.Z() };

			float u, v;
			Eigen::Vector3f p(xyz);
			project(p, u, v);

			fill(u, v);
		}
	}

	void fill(const pcl::PointCloud<point_t>& cloud) {
		for (auto it = cloud.begin(); it != cloud.end(); ++it) {
			float u, v;
			Eigen::Vector3f p(it->data);
			project(p, u, v);

			fill(u, v);
		}
	}

	template <typename T>
	void segment(float treshold, const T& cloud, std::vector<int>& inliers, std::vector<float>& distances) {
		// for (auto it = guess.begin(); it != guess.end(); ++it) {
		// 	const int I = *it;
		// 	const point_t& P = cloud[*it];

		for (auto it = cloud.begin(); it != cloud.end(); ++it) {
			// const int I = std::distance(cloud.begin(), it);
			const int I = calc_distance(cloud.begin(), it);
			const point_t& P = *it;

			const float d = distance(P);

			if (d <= treshold) {
				float u, v;
				Eigen::Vector3f p(P.data);
				project(p, u, v);

				if (u >= uv_min[0] && u <= uv_max[0]) {
					if (v >= uv_min[1] && v <= uv_max[1]) {
						inliers.push_back(I);
						distances.push_back(d);
					}
				}
			}
		}
	}

	void project(const Eigen::Vector3f& p, float& u, float& v) {
		Eigen::Vector3f op = p - pos;
		u = op.dot(x);
		v = op.dot(y);
	}

	void project(const Eigen::Vector3f& p, float& u, float& v, float& w) {
		Eigen::Vector3f op = p - pos;
		u = op.dot(x);
		v = op.dot(y);
		w = op.dot(z);
	}

	void write_transformed(std::ostream& stream, const pcl::PointCloud<point_t>& cloud) {
		unsigned int size = cloud.size();
		stream.write((char*)&size, sizeof(unsigned int));

		for (auto it = cloud.begin(); it != cloud.end(); ++it) {
			Eigen::Vector3f p(it->data);
			// const float x = it->x;
			// const float d = x - (-138.84408569);
			// const bool m = fabs(d < 1e-5);
			// stream.write((char*)&m, sizeof(bool));
			float uvw[3];
			project(p, uvw[0], uvw[1], uvw[2]);
			stream.write((char*)uvw, 3 * sizeof(float));
		}
	}

	float distance(const point_t& p) {
		return fabs(plane_coefficients[0] * p.x + plane_coefficients[1] * p.y + plane_coefficients[2] * p.z + plane_coefficients[3]);
	}

	friend std::ostream& operator<< (std::ostream& stream, const rectangular_trimmed_plane& plane) {
		float zero = 0.f, one = 1.f;
		stream.write((char*)plane.x.data(), sizeof(float) * 3);
		stream.write((char*)&zero, sizeof(float));
		stream.write((char*)plane.y.data(), sizeof(float) * 3);
		stream.write((char*)&zero, sizeof(float));
		stream.write((char*)plane.z.data(), sizeof(float) * 3);
		stream.write((char*)&zero, sizeof(float));
		stream.write((char*)plane.pos.data(), sizeof(float) * 3);
		stream.write((char*)&one, sizeof(float));

		return stream;
	}

	Eigen::Vector4f plane_coefficients;
	Eigen::Vector3f pos, x, y, z;
	Eigen::Vector2f uv_min, uv_max;
private:

	void init_bounds() {
		uv_min.x() = +std::numeric_limits<float>::infinity();
		uv_min.y() = +std::numeric_limits<float>::infinity();

		uv_max.x() = -std::numeric_limits<float>::infinity();
		uv_max.y() = -std::numeric_limits<float>::infinity();
	}

};

class face_def {
public:
	int hashcode;

	face_def(IfcParse::IfcFile* file, IfcSchema::IfcProduct* product, int face_id, const TopoDS_Face& face)
		: valid_(false)
		, product_(product)
		, face_id_(face_id)
		, face_(face)
		, area_(face_area(face_))
	{
		hashcode = face.HashCode(1 << 24);
		calculate(file);
	}

	face_def(IfcParse::IfcFile* file, std::istream& stream)
		: valid_(false)
		, product_(file->entityById(read_from_stream<int>(stream))->as<IfcSchema::IfcProduct>())
		, face_id_(read_from_stream<int>(stream))
		, face_(read_from_stream<TopoDS_Face>(stream))
		, area_(face_area(face_))
	{
		calculate();
	}

	/* face_def& operator=(const face_def & other) {
		product_ = other.product_;
		face_id_ = other.face_id_;
		face_ = other.face_;
		area_ = other.area_;
		return *this;
	} */

	bool operator<(const face_def& other) const {
		return area_ < other.area_;
	}

	friend std::ostream& operator<< (std::ostream& stream, const face_def& face_def) {
		BRep_Builder b;
		TopoDS_Compound c;
		b.MakeCompound(c);
		b.Add(c, face_def.face_);

		stream << face_def.product_->entity->id() << " " << face_def.face_id_ << std::endl;
		BRepTools_ShapeSet ss;
		ss.Add(c);
		ss.Write(stream);
		return stream;
	}

	// pcl::PointCloud<point_t>::Ptr
	void segment(float treshold, const lazy_pointcloud& cloud, std::vector<int>& inliers, std::vector<float>& distances) {
		
		/*
		std::cerr << "Face:" << std::endl;
		TopExp_Explorer exp(face_, TopAbs_VERTEX);
		for (; exp.More(); exp.Next()) {
			gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(exp.Current()));
			std::cerr << p.X() << " " << p.Y() << " " << p.Z() << std::endl;
		}
		*/

		Eigen::Vector2f uv_center = (pln_.uv_min + pln_.uv_max) / 2.;
		Eigen::Vector2f uv_extent = (pln_.uv_max - pln_.uv_min);
		Eigen::Vector3f xyz_center = pln_.pos + pln_.x * uv_center.x() + pln_.y * uv_center.y();

		// std::cerr << "Center:" << std::endl;
		// std::cerr << xyz_center.x() << " " << xyz_center.y() << " " << xyz_center.z() << std::endl;

		point_t point_center(xyz_center.x(), xyz_center.y(), xyz_center.z());
		float r = uv_extent.norm() / 2.;
		r = sqrt(r * r + treshold * treshold) + 1e-5;
		
		std::vector<int> all_inliers;
		std::vector<float> all_distances;
		if (r < 12) {
			auto Q = cloud.radius_search(point_center, r);
			// std::cerr << "Tree selection with radius " << r << ":" << Q.size() << "/" << cloud.size() << std::endl;
			pln_.segment(treshold, Q, all_inliers, all_distances);
		} else {
			pln_.segment(treshold, cloud, all_inliers, all_distances);
		}

		inliers.reserve(all_inliers.size());
		distances.reserve(all_inliers.size());

		ShapeAnalysis_Surface sas(BRep_Tool::Surface(face_));
		BRepTopAdaptor_FClass2d cls(face_, 1e-9);

		auto it = all_inliers.begin();
		auto it2 = all_distances.begin();

		for (; it != all_inliers.end();	++it, ++it2) {
			float u, v;
			const point_t& p = (*cloud.cloud())[*it];
			Eigen::Vector3f V(p.data);

			pln_.project(V, u, v);
			
			if (cls.Perform(gp_Pnt2d(u, v)) != TopAbs_OUT) {
				inliers.push_back(*it);
				distances.push_back(*it2);
			}
		}

		// std::cerr << "Topology selection: " << inliers.size() << "/" << all_inliers.size() << std::endl;
	}

	bool valid() const { return valid_; }
	rectangular_trimmed_plane pln() const { return pln_; }
	IfcSchema::IfcProduct* product() const { return product_; }
	double area() const { return area_; }
	int face_id() const { return face_id_; }

private:
	void calculate(IfcParse::IfcFile* file = 0) {
		Handle_Geom_Surface surf = BRep_Tool::Surface(face_);
		Handle_Geom_Plane plane = Handle_Geom_Plane::DownCast(surf);
		if (!plane.IsNull()) {
			gp_Pln pln = plane->Pln();
			// if (face_.Orientation() == TopAbs_REVERSED) {
			// 	pln.Mirror(pln.Location());
			// }
			pln_ = rectangular_trimmed_plane(pln);
			pln_.fill(face_);
			valid_ = true;
		}
	}
	
	bool valid_;
	IfcSchema::IfcProduct* product_;
	int face_id_;
	TopoDS_Face face_;
	double area_;
	rectangular_trimmed_plane pln_;
};

void load_compressed_pcd(const std::string& fn, pcl::PointCloud<point_t>::Ptr& cloud) {
	char filename[L_tmpnam];
	std::tmpnam(filename);

	std::stringstream ss;
	std::ifstream file(fn.c_str(), std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
	in.push(boost::iostreams::gzip_decompressor());
	in.push(file);
	try {
		ofstream f(filename);
		boost::iostreams::copy(in, f);
	} catch (const boost::exception& e) {
		std::cerr << boost::diagnostic_information(e) << std::endl;
	}

	pcl::io::loadPCDFile<point_t>(filename, *cloud);

	std::remove(filename);
}

Eigen::Matrix4f load_matrix(const std::string& fn) {
	std::ifstream ifs(fn.c_str());
	std::stringstream buffer;
	buffer << ifs.rdbuf();
	std::string s = buffer.str();

	boost::replace_all(s, "[", " ");
	boost::replace_all(s, "]", " ");
	boost::replace_all(s, ",", " ");

	std::istringstream ss(s);
	float data[16];
	for (int i = 0; i < 16; ++i) {
		ss >> data[i];
	}

	Eigen::Matrix4f M(data);
	M.transposeInPlace();

	return M;
}

void dump_binary_file(std::ostream& fs, const pcl::PointCloud<point_t>& cloud) {
	for (auto it = cloud.begin(); it != cloud.end(); ++it) {
		fs.write((const char*)it->data, sizeof(float) * 3);
	}
}

void dump_binary_file(const std::string& fn, const pcl::PointCloud<point_t>& cloud, bool append=false) {
	int mode = std::ios_base::binary;
	if (append) {
		mode = mode | std::ios_base::app;
	}
	std::ofstream fs(fn.c_str(), mode);
	dump_binary_file(fs, cloud);
}

static int SEGMENTATION_SUBSAMPLING = 1;

void load_binary_file(const std::string& fn, pcl::PointCloud<point_t>& cloud) {
	std::ifstream f(fn.c_str(), std::ios_base::binary);
	f.seekg(0, std::ios::end);
	size_t s1 = cloud.size();
	size_t s2 = f.tellg() / sizeof(float) / 3 / SEGMENTATION_SUBSAMPLING;
	f.seekg(0, std::ios::beg);
	cloud.resize(s1 + s2);
	auto it = cloud.begin() + s1;
	for (; it < cloud.end(); it += SEGMENTATION_SUBSAMPLING) {
		f.read((char*)it->data, sizeof(float) * 3);
	}
}

int main_x(int argc, char** argv) {
	const std::string fn = argv[1];
	const std::string dump = fn + ".bin";
	
	pcl::PointCloud<point_t>::Ptr cloud(new pcl::PointCloud<point_t>);
	pcl::io::loadPCDFile(fn, *cloud);
	dump_binary_file(dump, *cloud);

	return 0;
}

class labelled_index_remover : public std::unary_function<int, bool> {
public:
	labelled_index_remover(const std::vector<bool>& labeled)
		: labeled(labeled)
	{}
	bool operator ()(int i) { return labeled[i]; }
private:
	const std::vector<bool>& labeled;
};

class smaller_than : public std::unary_function<const pcl::PointIndices&, bool> {
public:
	smaller_than(int n) : n(n) {}
	bool operator ()(const pcl::PointIndices& pi) { return pi.indices.size() < n; }
private:
	const int n;
};

template <typename T>
size_t copy_if_not_in(const std::vector<T>& source, std::vector<T>& dest, const std::vector<T>& cmp) {
	size_t n = 0;
	T p = T();
	for (typename std::vector<T>::const_iterator it = source.begin(); it != source.end(); ++it) {
		if (it != source.begin() && *it == p) continue;
		p = *it;
		if (!std::binary_search(cmp.begin(), cmp.end(), *it)) {
			dest.push_back(*it);
			++n;
		}
	}
	return n;
}

std::vector<pcl::PointIndices> euclidian_cluster(pcl::PointCloud<point_t>::Ptr& cloud, double tolerance, int min_size) {
	pcl::search::KdTree<point_t>::Ptr tree(new pcl::search::KdTree<point_t>);
	tree->setInputCloud(cloud);

	std::vector<bool> labeled(cloud->size(), false);

	std::vector<pcl::PointIndices> clusters;
	
	for (auto it = cloud->begin(); it != cloud->end(); ++it) {
		auto i = std::distance(cloud->begin(), it);
		if (labeled[i]) continue;

		pcl::PointIndices indices;
		std::vector<float> distances;

		size_t indices_offset = 0;
		indices.indices.push_back(i);

		while (true) {
			std::vector<int> augmented_indices;
			for (auto jt = indices.indices.begin() + indices_offset; jt != indices.indices.end(); ++jt) {
				std::vector<int> tmp; // radiusSearch() clears vector (why??)
				tree->radiusSearch(*jt, tolerance, tmp, distances);
				augmented_indices.reserve(augmented_indices.size() + tmp.size());
				// radiusSearch contains query point as well?
				// std::copy(tmp.begin(), tmp.end(), std::back_inserter(augmented_indices));
				for (auto it = tmp.begin(); it != tmp.end(); ++it) {
					if (*it != *jt) {
						augmented_indices.push_back(*it);
					}
				}
			}
			indices_offset = indices.indices.size();

			auto truncated = std::remove_if(augmented_indices.begin(), augmented_indices.end(), labelled_index_remover(labeled));
			augmented_indices.erase(truncated, augmented_indices.end());

			std::sort(indices.indices.begin(), indices.indices.end());
			std::sort(augmented_indices.begin(), augmented_indices.end());

			if (copy_if_not_in(augmented_indices, indices.indices, indices.indices) == 0) {
				break;
			}
		}

		for (auto jt = indices.indices.begin(); jt != indices.indices.end(); ++jt) {
			labeled[*jt] = true;
		}

		clusters.push_back(indices);
		
	}

	std::cerr << "n " << min_size << "c " << clusters.size() << " ";
	auto truncated = std::remove_if(clusters.begin(), clusters.end(), smaller_than(min_size));
	clusters.erase(truncated, clusters.end());
	std::cerr << clusters.size() << std::endl;

	return clusters;
}

pointcloud_pair_t planar_segmentation(const std::string& model_name, pcl::PointCloud<point_t>::Ptr& cloud_complete) {
	const std::string segmented_path = "intermediate_files\\" + model_name + "\\planar_segmentation.bin";
	const std::string unsegmented_path = "intermediate_files\\" + model_name + "\\unsegmented";
	const std::string voxel_clustering_path = "intermediate_files\\" + model_name + "\\voxel_clustering.bin";

	std::cerr << "Subsampling..." << std::endl;

	voxel_mapping voxelizer(cloud_complete, 0.1);
	auto cluster_indices = voxelizer.components();
	int min_size = 0;
	{
		pcl::PointCloud<point_t>::Ptr cloud = voxelizer.applyFilter();
		std::cerr << cloud_complete->size() << " -> " << cloud->size() << " points. With " << cluster_indices.size() << " clusters." << std::endl;
		min_size = (int)pow((double)cloud->size(), 0.5);
	}

	// size_t subset_size_sum = 0;
	// for (auto it = cluster_indices.begin(); it != cluster_indices.end(); ++it) {
	// 	voxel_subset subset(voxelizer, *it);
	// 	subset_size_sum += subset.size();
	// }
	// std::cerr << cloud_complete->size() << " " << subset_size_sum << std::endl;

	// Haus 30
	// const double tolerance = 20. / 100.;
	// Byg 72
	const double tolerance = 20. / 100.;
	const int iterations = 1000;
	const double distance = 0.15;

	/*
	std::cerr << "Building tree..." << std::endl;

	pcl::search::KdTree<point_t>::Ptr tree(new pcl::search::KdTree<point_t>);
	tree->setInputCloud(cloud);

	std::cerr << cloud->size() << std::endl;
	
	// auto t0 = time(0);
	
	std::cerr << "Finding clusters..." << std::endl;

	std::vector<pcl::PointIndices> cluster_indices;
	pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
	ec.setClusterTolerance(tolerance);
	ec.setMinClusterSize(min_size);
	ec.setSearchMethod(tree);
	ec.setInputCloud(cloud);
	ec.extract(cluster_indices);
	
	// auto cluster_indices = euclidian_cluster(cloud, tolerance, min_size);

	// auto t1 = time(0);

	// std::cerr << "Found " << cluster_indices.size() << " clusters in " << (t1 - t0) << " // " << cloud->size() << std::endl;
	*/
	/*
	size_t S = 0;
	std::vector<int> reg(cloud_complete->size(), 0);
	for (auto it = cluster_indices.begin(); it != cluster_indices.end(); ++it) {
		auto I = std::distance(cluster_indices.begin(), it);
		voxel_subset subset(voxelizer, *it);
		S += subset.size();
		size_t s = 0;
		for (auto jt = subset.begin(); jt != subset.end(); ++jt, ++s) {
			reg[*jt] += 1;
		}
		if (s != subset.size() || s > 10000) {
			std::cerr << I << " " << s << " " << subset.size() << std::endl;
		}
	}

	std::cerr << cloud_complete->size() << " " << S << std::endl;

	for (auto it = reg.begin(); it != reg.end(); ++it) {
		if (*it != 1) {
			std::cerr << std::distance(reg.begin(), it) << ":" << (*it) << std::endl;
		}
	}

	return pointcloud_pair_t(nullptr, nullptr);
	*/
	std::ofstream ofs(segmented_path.c_str(), std::ios_base::binary);

	std::ofstream voxel_clustering(voxel_clustering_path.c_str(), std::ios_base::binary);

	unsigned threshold = min_size;

	size_t could_have_written = 0;
	size_t has_written = 0;
	size_t segmentated_count = 0;
	size_t unsegmentated_count = 0;

	{
		// Truncate the file
		const std::string fn = unsegmented_path + ".bin";
		std::ofstream fs(fn.c_str(), std::ios_base::binary);
	}	

	int n = 0;
	for (auto it = cluster_indices.begin(); it != cluster_indices.end(); ++it, ++n) {

		voxel_subset subset(voxelizer, *it);

		std::cerr << "Written: " << has_written << " / " << could_have_written << std::endl;

		could_have_written += subset.size();

		if (subset.size() < threshold) {
			// Don't split explicitly, but rather write streaming based on indices
			const std::string fn = unsegmented_path + ".bin";
			int mode = std::ios_base::binary | std::ios_base::app;
			std::ofstream fs(fn.c_str(), mode);
			for (auto jt = subset.begin(); jt != subset.end(); ++jt) {
				fs.write((const char*)(*cloud_complete)[*jt].data, sizeof(float) * 3);
				has_written += 1;
				unsegmentated_count += 1;
			}
			if (fs.tellp() / 12 > unsegmentated_count) {
				std::cerr << "fp: " << fs.tellp() / 12 << " " << unsegmentated_count << std::endl;
			}
			continue;
		}

		std::vector<int> indices_contiguous;
		indices_contiguous.reserve(subset.size());
		std::copy(subset.begin(), subset.end(), std::back_inserter(indices_contiguous));

		auto split_by_cluster = split_pointcloud(cloud_complete, indices_contiguous, true);
		auto cluster = split_by_cluster.first;
		auto initial_size = cluster->size();

		unsigned int size = indices_contiguous.size();
		voxel_clustering.write((char*)&size, sizeof(unsigned int));
		for (auto jt = cluster->begin(); jt != cluster->end(); ++jt) {
			voxel_clustering.write((char*)jt->data, sizeof(float) * 3);
		}

		while (cluster->size() > threshold) {

			pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients());
			pcl::PointIndices::Ptr inliers(new pcl::PointIndices());
			pcl::SACSegmentation<point_t> seg;
			seg.setModelType(pcl::SACMODEL_PLANE);
			seg.setMethodType(0);
			seg.setMaxIterations(iterations);
			seg.setDistanceThreshold(distance);
			seg.setInputCloud(cluster);
			seg.segment(*inliers, *coefficients);

			if (inliers->indices.size() < threshold) {
				break;
			}

			auto split_by_plane = split_pointcloud(cluster, inliers->indices);
			auto cluster_plane = split_by_plane.first;

			rectangular_trimmed_plane plane(Eigen::Vector4f(coefficients->values.data()));
			plane.fill(*cluster_plane);
			ofs << plane;

			plane.write_transformed(ofs, *cluster_plane);

			has_written += cluster_plane->size();
			segmentated_count += cluster_plane->size();

			std::cerr << "\rCluster " << n << ": " << initial_size << ". Remaining: " << cluster->size() << "      ";

			cluster = split_by_plane.second;
		}

		dump_binary_file(unsegmented_path + ".bin", *cluster, true);

		has_written += cluster->size();
		unsegmentated_count += cluster->size();
	}

	std::cerr << "Unseg: " << unsegmentated_count << " seg: " << segmentated_count << std::endl;
	
	return pointcloud_pair_t(nullptr, nullptr);
}

int cluster_sizes(const std::string& model_name, pcl::PointCloud<point_t>::Ptr& cloud) {
	pcl::search::KdTree<point_t>::Ptr tree(new pcl::search::KdTree<point_t>);
	tree->setInputCloud(cloud);

	std::ofstream fs("intermediate_files\\" + model_name + "\\clusters.csv");

	for (int j = 0; j < 3; ++j) {
		// const int min_sizes[3] = { 100, 1000, 10000 };
		// const int min_size = cloud->size() / min_sizes[j];
		// const int min_sizes[4] = { 1, 10, 100, 1000 };
		// const int min_size = min_sizes[j];
		const int min_size = (int) pow((double)cloud->size(), j * 0.25);

		fs << min_size;

		// for (int i = 1; i <= 100; i += ((i<10) ? 1 : ((i<50) ? 5 : 10))) {
		for (int i = 1; i <= 50; i += (i < 5 ? 1 : (i < 30 ? 5 : 10))) {
			const double tolerance = (double)i / 100.;
			
			std::vector<pcl::PointIndices> cluster_indices;
			pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
			ec.setClusterTolerance(tolerance);
			ec.setMinClusterSize(min_size);
			// ec.setMaxClusterSize(cloud->size() / 10);
			ec.setSearchMethod(tree);
			ec.setInputCloud(cloud);
			ec.extract(cluster_indices);
			
			// cluster_indices = euclidian_cluster(cloud, tolerance, min_size);

			std::cerr << i << ": " << cluster_indices.size() << " ";
			fs << "," << i << "," << cluster_indices.size();
		}
		fs << std::endl;
	}
	std::cin.get();
}

void load_ifc_from_model(IfcParse::IfcFile& file, std::vector<face_def>& faces) {
	IfcGeom::IteratorSettings settings;
	settings.set(IfcGeom::IteratorSettings::DISABLE_TRIANGULATION, true);
	IfcGeom::Iterator<double> context(settings, &file);

	std::set<std::string> es;
	es.insert("IfcWall");
	es.insert("IfcSlab");
	es.insert("IfcStairFlight");
	es.insert("IfcStair");
	context.includeEntities(es);

	if (!context.initialize()) {
		return;
	}

	do {
		IfcGeom::BRepElement<double>* elem = (IfcGeom::BRepElement<double>*) context.get();
		IfcSchema::IfcProduct* product = context.getFile()->entityById(elem->id())->as<IfcSchema::IfcProduct>();

		int face_id = 0;
		for (auto it = elem->geometry().begin(); it != elem->geometry().end(); ++it) {
			TopoDS_Shape shape;
			if (it->Placement().Form() == gp_Other) {
				shape = BRepBuilderAPI_GTransform(it->Shape(), it->Placement(), true);
			} else {
				shape = BRepBuilderAPI_Transform(it->Shape(), it->Placement().Trsf(), true);
			}
			shape = BRepBuilderAPI_Transform(shape, elem->transformation().data(), true);

			TopExp_Explorer exp(shape, TopAbs_FACE);
			for (; exp.More(); exp.Next(), ++face_id) {
				const TopoDS_Face& face = TopoDS::Face(exp.Current());
				faces.push_back(face_def(&file, product, face_id, face));
			}
		}
	} while (context.next());

	std::sort(faces.begin(), faces.end());
	std::reverse(faces.begin(), faces.end());
}

void load_ifc_from_cache(IfcParse::IfcFile& file, std::ifstream& stream, std::vector<face_def>& faces) {
	std::cerr << "Loading IFC from cache" << std::endl;

	while (true) {
		try {
			face_def fs(&file, stream);
			faces.push_back(fs);
		} catch (stop_iteration&) {
			break;
		}
	}
}

std::vector<face_def> load_ifc(IfcParse::IfcFile& file, const std::string& ifc_cache_path, int num_faces) {
	std::ifstream stream(ifc_cache_path.c_str());

	std::vector<face_def> faces;
	faces.reserve(file.entitiesByType<IfcSchema::IfcProduct>()->size() * 8);

	if (stream.good()) {
		load_ifc_from_cache(file, stream, faces);
	} else {
		load_ifc_from_model(file, faces);

		std::ofstream ofs(ifc_cache_path.c_str());
		for (auto it = faces.begin(); it != faces.end(); ++it) {
			ofs << *it;
		}
	}

	if (faces.size() > num_faces) {
		faces.erase(faces.begin() + num_faces, faces.end());
	}

	return faces;
}

int load_ascii_scans(const std::string& path, const Eigen::Matrix4f& matrix, pcl::PointCloud<point_t>::Ptr& cloud, int offset=0) {
	int stopped_at = -1;
	for (int N = offset;; ++N) {
		// TODO: How much, this is single scan?
		if (cloud->size() > 20 * 1000 * 1000) {
			stopped_at = N;
			break;
		}

		std::cerr << "\rLoading scan: " << N << " Points: " << cloud->size();
		std::cerr.flush();

		std::stringstream ss;
		ss << path << N << ".pcd";
		std::string fn = ss.str();

		pcl::PointCloud<point_t>::Ptr cloud_part(new pcl::PointCloud<pcl::PointXYZ>);
		
		if (std::ifstream(fn.c_str()).good()) {
			pcl::io::loadPCDFile<point_t>(fn, *cloud_part);
		} else {
			ss << ".gz";
			fn = ss.str();
			if (std::ifstream(fn.c_str()).good()) {
				load_compressed_pcd(fn, cloud_part);
			} else {
				// Finished loading scans.
				stopped_at = -1;
				break;
			}
		}

		(*cloud) += *cloud_part;
	}

	std::cerr << "\rLoaded " << cloud->size() << " points                 " << std::endl;
	std::cerr << "Applying transformation" << std::endl;

	float d[4] = { 0., 0., 0., 1. };
	for (auto it = cloud->begin(); it != cloud->end(); ++it) {
		memcpy(d, it->data, sizeof(float) * 3);
		Eigen::Vector4f V = matrix * Eigen::Vector4f(d);
		memcpy(it->data, V.data(), sizeof(float) * 3);
	}

	return stopped_at;
}

pointcloud_pair_t associate_points(const std::string& model_name) {
	// File paths
	const std::string ifc_path = "input_data\\" + model_name + ".ifc";
	const std::string matrix_path = "input_data\\" + model_name + "_alignment_matrix.json";
	const std::string ifc_cache_path = "intermediate_files\\" + model_name + ".ifc.cache.txt";
	const std::string ascii_pcd_path = "input_data\\" + model_name + "\\";
	const std::string binary_pcd_path = ascii_pcd_path + "pointcloud.pcdx";
	const std::string association_results_path = "intermediate_files\\" + model_name + "\\association_results.bin";
	const std::string distances_path = "intermediate_files\\" + model_name + "\\";
	const std::string associated_points_path = "intermediate_files\\" + model_name + "\\associated";
	const std::string unassociated_points_path = "intermediate_files\\" + model_name + "\\unassociated";

	// Constants
	const int num_faces_to_associate = 1000;
	const unsigned NUM_TRESHOLDS = 3;
	float association_thresholds[NUM_TRESHOLDS] = { 0.01f, 0.05f, 0.15f };
	
	// float association_thresholds[NUM_TRESHOLDS] = { 0.1f, 0.2f, 0.5f };
	// Used for distance density plot:
	// const unsigned NUM_TRESHOLDS = 4;
	// float association_thresholds[NUM_TRESHOLDS] = { 0.01f, 0.03f, 0.09f, 0.27f };

	// Initialisation
	Eigen::Matrix4f M = load_matrix(matrix_path);

	IfcParse::IfcFile file;
	file.Init(ifc_path);

	std::vector<face_def> faces = load_ifc(file, ifc_cache_path, num_faces_to_associate);

	pcl::PointCloud<point_t>::Ptr cloud_(new pcl::PointCloud<point_t>);
	pcl::PointCloud<point_t>::Ptr all_associated(new pcl::PointCloud<point_t>);

	std::set<IfcSchema::Type::Enum> written_distances;

	const bool binary_pcd_fn_exists = std::ifstream(binary_pcd_path.c_str()).good();

	std::ofstream association_results(association_results_path.c_str(), std::ios_base::binary);

	int scan_number = 0;
	bool multiple_batches = false;
	bool first_batch = true;

next_batch:
	if (binary_pcd_fn_exists) {
		pcl::io::loadPCDFile(binary_pcd_path, *cloud_);
		std::cerr << "Loaded " << cloud_->size() << " points" << std::endl;
	} else {
		scan_number = load_ascii_scans(ascii_pcd_path, M, cloud_, scan_number);
		if (scan_number != -1) {
			multiple_batches = true;
		} else if (!multiple_batches) {
			pcl::io::savePCDFileBinaryCompressed(binary_pcd_path, *cloud_);
		}
	}

	if (!cloud_->empty()) {
		lazy_pointcloud cloud(cloud_);

		for (int i = 0; i < NUM_TRESHOLDS; ++i) {

			for (auto it = faces.begin(); it != faces.end(); ++it) {

				if (it->valid()) {
					std::cerr << "\rAssociation pass: " << i << " Remaining: " << cloud.size() << "     ";
					std::cerr.flush();

					std::vector<float> distances;
					pcl::IndicesPtr inliers(new std::vector<int>);
					it->segment(association_thresholds[i], cloud, *inliers, distances);

					if (inliers->size() > 32) {

						// auto splitted = split_pointcloud(cloud, *inliers);
						auto splitted = cloud.split(inliers);

						{
							int open_mode = std::ios_base::binary;
							if (written_distances.find(it->product()->type()) == written_distances.end()) {
								written_distances.insert(it->product()->type());
							} else {
								open_mode |= std::iostream::app;
							}
							const std::string fn = distances_path + IfcSchema::Type::ToString(it->product()->type()) + ".bin";
							std::ofstream fs(fn.c_str(), open_mode);
							fs.write((const char*)distances.data(), sizeof(float) * distances.size());
						}

						cloud = splitted.second;
						auto just_associated = splitted.first.extract();
						*all_associated += *just_associated;

						uint32_t temp;
						temp = it->product()->entity->id();
						association_results.write((char*)&temp, sizeof(uint32_t));

						temp = it->face_id();
						association_results.write((char*)&temp, sizeof(uint32_t));

						association_results << it->pln();
						it->pln().write_transformed(association_results, *just_associated);
					}
				}
			}

			std::cerr << std::endl;

			if (i < NUM_TRESHOLDS - 1) {
				// Force to rebuild the tree to speed up lookup
				cloud = lazy_pointcloud(cloud.extract());
			}
		}

		if (multiple_batches) {
			dump_binary_file(unassociated_points_path + ".bin", *cloud.extract(), !first_batch);

			cloud_->clear();
			all_associated->clear();

			if (scan_number != -1) {
				// Yay, evil
				first_batch = false;
				goto next_batch;
			}
		} else {
			dump_binary_file(unassociated_points_path + ".bin", *cloud.extract());
			pcl::io::savePCDFileBinaryCompressed(unassociated_points_path + ".pcdx", *cloud.extract());

			dump_binary_file(associated_points_path + ".bin", *all_associated);
			pcl::io::savePCDFileBinaryCompressed(associated_points_path + ".pcdx", *all_associated);
		}
	}
	
	return pointcloud_pair_t(all_associated, cloud_);
}

int main(int argc, char** argv) {
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <model>" << std::endl;
		return 1;
	}

	const std::string model_name = argv[1];
	// associate_points(model_name);
	
	const std::string unassociated_points_path = "intermediate_files\\" + model_name + "\\unassociated";
	pcl::PointCloud<point_t>::Ptr unassociated(new pcl::PointCloud<point_t>);
	load_binary_file(unassociated_points_path + ".bin", *unassociated);
	planar_segmentation(model_name, unassociated);
	
	// cluster_sizes(model_name, unassociated);
	// cluster_sizes(model_name, association.second);
	// planar_segmentation(model_name, association.second);

	std::cin.get();
}