#pragma once

#include <deque>
#include <iterator>
#include <limits>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Vector_3 Vector;
typedef std::array<K::FT, 3> Point;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor HalfEdgeDescriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor EdgeDescriptor;

struct HalfEdge2Edge {
    HalfEdge2Edge(const Mesh& m, std::vector<EdgeDescriptor>& edges) : m_mesh(m), m_edges(edges) {
    }

    void operator()(const HalfEdgeDescriptor& h) const {
        m_edges.push_back(edge(h, m_mesh));
    }

    const Mesh& m_mesh;
    std::vector<EdgeDescriptor>& m_edges;
};

struct PointArray_traits {
    struct Equal_3 {
        bool operator()(const Point& p, const Point& q) const {
          return (p == q);
        }
    };

    struct Less_xyz_3 {
        bool operator()(const Point& p, const Point& q) const {
          return std::lexicographical_compare(p.begin(), p.end(), q.begin(), q.end());
        }
    };

    Equal_3 equal_3_object() const { return Equal_3(); }

    Less_xyz_3 less_xyz_3_object() const { return Less_xyz_3(); }
};

struct Point_traits {
    struct Equal_3 {
        bool operator()(const K::Point_3& p, const K::Point_3& q) const {
          return (p == q);
        }
    };

    struct Less_xyz_3 {
        bool operator()(const K::Point_3& p, const K::Point_3& q) const {
          return std::lexicographical_compare(p.cartesian_begin(), p.cartesian_end(),
                                              q.cartesian_begin(), q.cartesian_end());
        }
    };

    Equal_3 equal_3_object() const { return Equal_3(); }

    Less_xyz_3 less_xyz_3_object() const { return Less_xyz_3(); }
};

template <typename Traits, typename PointRange, typename PolygonRange>
std::vector<std::size_t> find_invalid_polygons_in_polygon_soup(const PointRange& points,
                                                               const PolygonRange& polygons,
                                                               const Traits& traits = Traits())
{
  typedef typename PMP::internal::Polygon_types<PointRange, PolygonRange>::Polygon_3 Polygon_3;

  std::vector<std::size_t> invalid;
  const std::size_t ini_polygons_size = polygons.size();
  for(std::size_t polygon_index=0; polygon_index!=ini_polygons_size; ++polygon_index) {
    const Polygon_3& polygon = polygons[polygon_index];
    const std::size_t N = polygon.size(), last = N-1;
    if (N < 3) {
      invalid.push_back(polygon_index);
      continue;
    }

    for(std::size_t i=0; i<N; ++i) {
      const std::size_t next_i = (i == last) ? 0 : i+1;
      if(polygon[i] == polygon[next_i] || // combinatorial equality
         traits.equal_3_object()(points[polygon[i]], points[polygon[next_i]])) // geometric equality
      {
        invalid.push_back(polygon_index);
        break;
      }
    }
  }

  return invalid;
}

template <typename PointRange, typename PolygonRange, typename NamedParameters = CGAL::parameters::Default_named_parameters>
std::vector<std::size_t> find_duplicate_polygons_in_polygon_soup(const PointRange& points,
                                                                 const PolygonRange& polygons,
                                                                 const NamedParameters& np = CGAL::parameters::default_values())
{
  typedef typename PMP::internal::GetPolygonGeomTraits<PointRange, PolygonRange, NamedParameters>::type Traits;
  Traits traits = CGAL::parameters::choose_parameter<Traits>(CGAL::parameters::get_parameter(np, CGAL::internal_np::geom_traits));

  typedef typename PMP::internal::Polygon_types<PointRange, PolygonRange>::P_ID P_ID;

  std::deque<std::vector<P_ID> > all_duplicate_polygons;
  PMP::internal::collect_duplicate_polygons(points, polygons, std::back_inserter(all_duplicate_polygons), traits, false);

  if (all_duplicate_polygons.empty()) {
    return {};
  }

  std::vector<std::size_t> duplicates;
  for (const auto& duplicate_set: all_duplicate_polygons) {
    if (duplicate_set.size() > 1) {
      duplicates.insert(duplicates.end(), duplicate_set.begin() + 1, duplicate_set.end());
    }
  }

  std::sort(duplicates.begin(), duplicates.end());

  return duplicates;
}

template <typename T>
void remove_indices(std::vector<T>& data, const std::vector<size_t>& to_remove) {
  if (to_remove.empty()) {
      return;
  }

  auto ascending_remove_indices = to_remove;
  std::sort(ascending_remove_indices.begin(), ascending_remove_indices.end());
  while (ascending_remove_indices.back() >= data.size() && !ascending_remove_indices.empty()) {
    ascending_remove_indices.pop_back();
  }

  std::vector<T> output;
  output.reserve(data.size());

  size_t idx = 0;
  for (size_t next_remove_idx: ascending_remove_indices) {
    if (idx != next_remove_idx) {
      output.insert(output.end(), data.begin() + idx, data.begin() + next_remove_idx);
    }
    idx = next_remove_idx + 1;
  }

  if (idx < data.size()) {
    output.insert(output.end(), data.begin() + idx, data.end());
  }

  data = output;
}

bool is_potentially_non_manifold(Mesh& mesh, Mesh::Vertex_index v1, Mesh::Vertex_index v2) {
    // Count halfedges from v1 to v2 and from v2 to v1
    int count_v1_to_v2 = 0;
    int count_v2_to_v1 = 0;

    // Iterate over halfedges around v1
    for (auto h : CGAL::halfedges_around_target(mesh.halfedge(v1), mesh)) {
        if (mesh.source(h) == v2) {
            count_v1_to_v2++;
        }
    }

    // Iterate over halfedges around v2
    for (auto h : CGAL::halfedges_around_target(mesh.halfedge(v2), mesh)) {
        if (mesh.source(h) == v1) {
            count_v2_to_v1++;
        }
    }

    // Non-manifold condition if any count exceeds 1
    return (count_v1_to_v2 > 1 || count_v2_to_v1 > 1);
}

template <typename Traits>
Mesh merge_meshes(const std::vector<Mesh>& meshes, const Traits& traits = Traits()) {
    Mesh merged;

    if (meshes.empty()) {
      return merged;
    }

    merged = meshes[0];

    if (meshes.size() == 1) {
        return merged;
    }

    // construct vertex lookup
    typedef typename Traits::Less_xyz_3 Less_xyz_3;
    typedef std::map<Mesh::Point, Mesh::Vertex_index, Less_xyz_3> VertexIndex;
    VertexIndex vertex_lookup(traits.less_xyz_3_object());
    for (const auto& v : meshes[0].vertices()) {
        vertex_lookup[meshes[0].point(v)] = v;
    }

    for (size_t i = 1; i < meshes.size(); i++) {
        const auto& mesh = meshes[i];
        for (const auto& f : mesh.faces()) {
            std::vector<Mesh::Vertex_index> face_vertices;
            for (const auto& v : mesh.vertices_around_face(mesh.halfedge(f))) {

                auto result = vertex_lookup.find(mesh.point(v));
                if (result == vertex_lookup.end()) {
                    auto index = merged.add_vertex(mesh.point(v));
                    vertex_lookup[mesh.point(v)] = index;
                    face_vertices.push_back(index);
                }
                else {
                  face_vertices.push_back(result->second);
                }
            }
            auto face_index = merged.add_face(face_vertices);
            // TODO(malban): check than face was added
        }
    }

    return merged;
}
