#ifndef H3_CORE_H
#define H3_CORE_H

#include <godot_cpp/classes/ref_counted.hpp>
#include <godot_cpp/variant/array.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/string.hpp>
#include <godot_cpp/variant/vector2.hpp>
#include <godot_cpp/variant/vector4.hpp>
#include "h3api.h"

namespace godot {
class H3Core : public RefCounted {
    GDCLASS(H3Core, RefCounted)

protected:
    static void _bind_methods();

public:
    H3Core() = default;
    ~H3Core() = default;

    // Indexing functions
    H3Index lat_lng_to_cell(double lat, double lon, int res) const;
    Vector2 cell_to_lat_lng(H3Index h3_index) const;
    Array cell_to_boundary(H3Index h3_index) const;
    PackedVector2Array cell_to_boundary_packed(H3Index h3_index) const;

    // Index inspection functions
    int get_resolution(H3Index h3_index) const;
    int get_base_cell_number(H3Index h3_index) const;
    H3Index string_to_h3(const String& str) const;
    String h3_to_string(const H3Index& h3_index) const;
    Array strings_to_h3(const Array& strs) const;
    Array h3s_to_string(const Array& h3_indices) const;
    PackedInt64Array strings_to_h3_packed(const PackedStringArray& strs) const;
    PackedStringArray h3s_to_string_packed(const PackedInt64Array& h3_indices) const;
    bool is_valid_cell(H3Index h3_index) const;
    bool is_res_class_iii(H3Index h3_index) const;
    bool is_pentagon(H3Index h3_index) const;
    Array get_icosahedron_faces(H3Index h3_index) const;
    PackedInt64Array get_icosahedron_faces_packed(H3Index h3_index) const;
    int max_face_count(H3Index h3_index) const;

    // Grid traversal functions
    int64_t grid_distance(H3Index origin, H3Index h3) const;
    Array grid_ring_unsafe(H3Index h3_index, int k) const;
    PackedInt64Array grid_ring_unsafe_packed(H3Index h3_index, int k) const;
    Array grid_disk(H3Index h3_index, int k) const;
    PackedInt64Array grid_disk_packed(H3Index h3_index, int k) const;
    int64_t max_grid_disk_size(int k) const;
    Array grid_disk_distances(H3Index h3_index, int k) const;
    Array grid_disk_unsafe(H3Index h3_index, int k) const;
    PackedInt64Array grid_disk_unsafe_packed(H3Index h3_index, int k) const;
    Array grid_disk_distances_unsafe(H3Index h3_index, int k) const;
    // Array grid_disks_unsafe(const Array& h3_set, int k) const;
    Array grid_path_cells(H3Index start, H3Index end) const;
    PackedInt64Array grid_path_cells_packed(H3Index start, H3Index end) const;
    int64_t grid_path_cells_size(H3Index start, H3Index end) const;
    Vector2i cell_to_local_ij(H3Index origin, H3Index h3_index) const;
    H3Index local_ij_to_cell(H3Index origin, const Vector2i& ij) const;

    // Hierarchical grid functions
    H3Index cell_to_parent(H3Index h3_index, int parent_res) const;
    Array cell_to_children(H3Index h3_index, int child_res) const;
    PackedInt64Array cell_to_children_packed(H3Index h3_index, int child_res) const;
    int64_t cell_to_child_size(H3Index h3_index, int child_res) const;
    H3Index cell_to_center_child(H3Index h3_index, int child_res) const;
    int64_t cell_to_child_pos(H3Index child, int parent_res) const;
    H3Index child_pos_to_cell(int64_t child_pos, H3Index parent, int child_res) const;
    Array compact_cells(const Array& h3_set) const;
    PackedInt64Array compact_cells_packed(const PackedInt64Array& h3_set) const;
    Array uncompact_cells(const Array& compacted_set, int res) const;
    PackedInt64Array uncompact_cells_packed(const PackedInt64Array& compacted_set, int res) const;
    int64_t uncompact_cells_size(const Array& compacted_set, int res) const;
    int64_t uncompact_cells_size_packed(const PackedInt64Array& compacted_set, int res) const;

    // Region functions
    Array polygon_to_cells(const Array& polygon, int res, uint32_t flags, bool is_geojson = false) const;
    PackedInt64Array polygon_to_cells_packed(const Array& polygon, int res, uint32_t flags, bool is_geojson = false) const;
    int64_t max_polygon_to_cells_size(const Array& polygon, int res, uint32_t flags, bool is_geojson = false) const;
    Array polygon_to_cells_experimental(const Array& polygon, int res, uint32_t flags, bool is_geojson = false) const;
    PackedInt64Array polygon_to_cells_experimental_packed(const Array& polygon, int res, uint32_t flags, bool is_geojson = false) const;
    int64_t max_polygon_to_cells_experimental_size(const Array& polygon, int res, uint32_t flags, bool is_geojson = false) const;
    Array cells_to_linked_multi_polygon(const Array& h3_set, bool is_geojson = false) const;

    // Directed edge functions
    bool are_neighbor_cells(H3Index origin, H3Index destination) const;
    H3Index cells_to_directed_edge(H3Index origin, H3Index destination) const;
    bool is_valid_directed_edge(const H3Index edge) const;
    H3Index get_directed_edge_origin(const H3Index edge) const;
    H3Index get_directed_edge_destination(const H3Index edge) const;
    Array directed_edge_to_cells(const H3Index edge) const;
    PackedInt64Array directed_edge_to_cells_packed(const H3Index edge) const;
    Array origin_to_directed_edges(H3Index origin) const;
    PackedInt64Array origin_to_directed_edges_packed(H3Index origin) const;
    Array directed_edge_to_boundary(const H3Index edge) const;
    PackedVector2Array directed_edge_to_boundary_packed(const H3Index edge) const;

    // Vertex functions
    H3Index cell_to_vertex(H3Index origin, int vertex_num) const;
    Array cell_to_vertices(H3Index origin) const;
    PackedInt64Array cell_to_vertices_packed(H3Index origin) const;
    Vector2 vertex_to_lat_lng(const H3Index vertex) const;
    bool is_valid_vertex(const H3Index vertex) const;

    // Miscellaneous H3 functions
    double degs_to_rads(double degrees) const;
    double rads_to_degs(double radians) const;
    double get_hexagon_area_avg_km2(int res) const;
    double get_hexagon_area_avg_m2(int res) const;
    double cell_area_rads2(H3Index h3_index) const;
    double cell_area_km2(H3Index h3_index) const;
    double cell_area_m2(H3Index h3_index) const;
    double get_hexagon_edge_length_avg_km(int res) const;
    double get_hexagon_edge_length_avg_m(int res) const;
    double edge_length_km(const H3Index edge) const;
    double edge_length_m(const H3Index edge) const;
    double edge_length_rads(const H3Index edge) const;
    int64_t get_num_cells(int res) const;
    Array get_res0_cells() const;
    PackedInt64Array get_res0_cells_packed() const;
    int res0_cell_count() const;
    Array get_pentagons(int res) const;
    PackedInt64Array get_pentagons_packed(int res) const;
    int pentagon_count() const;
    double great_circle_distance_km(const Vector2 point1, const Vector2 point2) const;
    double great_circle_distance_m(const Vector2 point1, const Vector2 point2) const;
    double great_circle_distance_rads(const Vector2 point1, const Vector2 point2) const;

    // Godot functions
    Vector4 h3_to_vector4(H3Index h3_index) const;
    H3Index vector4_to_h3(const Vector4& vec) const;
};

} // namespace godot

#endif // H3_CORE_H