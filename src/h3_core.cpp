#include "h3_core.h"
#include "h3api.h"
#include <algorithm>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/error_macros.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <math.h>

#define EARTH_RADIUS_KM 6371.0
#define M_PI 3.14159265358979323846

using namespace godot;

static bool handle_h3_error(H3Error error) {
    if (error != E_SUCCESS) {
        const char* error_str = H3_EXPORT(describeH3Error)(error);
        ERR_PRINT(String("H3 Error: ") + error_str);
        return false;
    }
    return true;
}

static GeoPolygon array_to_geopolygon(const Array& polygon_array, bool is_geojson) {
    GeoPolygon polygon;
    Array normalized_array;

    if (polygon_array.size() > 0) {
        if (polygon_array[0].get_type() == Variant::VECTOR2) {
            normalized_array.append(polygon_array);
        } else {
            normalized_array = polygon_array;
        }
    }

    if (normalized_array.is_empty()) {
        polygon.geoloop.numVerts = 0;
        polygon.geoloop.verts = nullptr;
        polygon.numHoles = 0;
        polygon.holes = nullptr;
        return polygon;
    }

    Array perimeter = normalized_array[0];
    polygon.geoloop.numVerts = perimeter.size();
    polygon.geoloop.verts = new LatLng[perimeter.size()];

    for (int i = 0; i < perimeter.size(); i++) {
        Vector2 vertex = perimeter[i];
        if (is_geojson) {
            polygon.geoloop.verts[i].lat = H3_EXPORT(degsToRads)(vertex.y);
            polygon.geoloop.verts[i].lng = H3_EXPORT(degsToRads)(vertex.x);
        } else {
            polygon.geoloop.verts[i].lat = H3_EXPORT(degsToRads)(vertex.x);
            polygon.geoloop.verts[i].lng = H3_EXPORT(degsToRads)(vertex.y);
        }
    }

    polygon.numHoles = normalized_array.size() - 1;
    if (polygon.numHoles > 0) {
        polygon.holes = new GeoLoop[polygon.numHoles];

        for (int h = 0; h < polygon.numHoles; h++) {
            Array hole = normalized_array[h + 1];
            polygon.holes[h].numVerts = hole.size();
            polygon.holes[h].verts = new LatLng[hole.size()];

            for (int i = 0; i < hole.size(); i++) {
                Vector2 vertex = hole[i];
                if (is_geojson) {
                    polygon.holes[h].verts[i].lat = H3_EXPORT(degsToRads)(vertex.y);
                    polygon.holes[h].verts[i].lng = H3_EXPORT(degsToRads)(vertex.x);
                } else {
                    polygon.holes[h].verts[i].lat = H3_EXPORT(degsToRads)(vertex.x);
                    polygon.holes[h].verts[i].lng = H3_EXPORT(degsToRads)(vertex.y);
                }
            }
        }
    } else {
        polygon.holes = nullptr;
    }

    return polygon;
}

static void free_geopolygon(GeoPolygon& polygon) {
    delete[] polygon.geoloop.verts;

    if (polygon.holes) {
        for (int h = 0; h < polygon.numHoles; h++) {
            delete[] polygon.holes[h].verts;
        }
        delete[] polygon.holes;
    }
}

static Array linked_geo_polygon_to_dict(const LinkedGeoPolygon* polygon, bool is_geojson) {
    Array result;

    const LinkedGeoPolygon* currentPolygon = polygon;
    while (currentPolygon != nullptr) {
        Array loop;
        const LinkedGeoLoop* currentLoop = currentPolygon->first;

        while (currentLoop != nullptr) {
            Array vertices;
            const LinkedLatLng* currentVertex = currentLoop->first;

            while (currentVertex != nullptr) {
                Vector2 vertex;
                if (is_geojson) {
                    vertex.x = H3_EXPORT(radsToDegs)(currentVertex->vertex.lng);
                    vertex.y = H3_EXPORT(radsToDegs)(currentVertex->vertex.lat);
                } else {
                    vertex.x = H3_EXPORT(radsToDegs)(currentVertex->vertex.lat);
                    vertex.y = H3_EXPORT(radsToDegs)(currentVertex->vertex.lng);
                }
                vertices.append(vertex);
                currentVertex = currentVertex != currentLoop->last ? currentVertex->next : nullptr;
            }

            loop.append(vertices);
            currentLoop = currentLoop != currentPolygon->last ? currentLoop->next : nullptr;
        }

        result.append(loop);
        currentPolygon = currentPolygon->next;
    }

    return result;
}

void H3Core::_bind_methods() {
    // Indexing Functions
    ClassDB::bind_method(D_METHOD("lat_lng_to_cell", "lat", "lon", "res"), &H3Core::lat_lng_to_cell);
    ClassDB::bind_method(D_METHOD("cell_to_lat_lng", "h3_index"), &H3Core::cell_to_lat_lng);
    ClassDB::bind_method(D_METHOD("cell_to_boundary", "h3_index"), &H3Core::cell_to_boundary);
    ClassDB::bind_method(D_METHOD("cell_to_boundary_packed", "h3_index"), &H3Core::cell_to_boundary_packed);

    // Index inspection functions
    ClassDB::bind_method(D_METHOD("get_resolution", "h3_index"), &H3Core::get_resolution);
    ClassDB::bind_method(D_METHOD("get_base_cell_number", "h3_index"), &H3Core::get_base_cell_number);
    ClassDB::bind_method(D_METHOD("string_to_h3", "h3_str"), &H3Core::string_to_h3);
    ClassDB::bind_method(D_METHOD("h3_to_string", "h3_index"), &H3Core::h3_to_string);
    ClassDB::bind_method(D_METHOD("strings_to_h3", "strs"), &H3Core::strings_to_h3);
    ClassDB::bind_method(D_METHOD("h3s_to_string", "h3_indices"), &H3Core::h3s_to_string);
    ClassDB::bind_method(D_METHOD("strings_to_h3_packed", "strs"), &H3Core::strings_to_h3_packed);
    ClassDB::bind_method(D_METHOD("h3s_to_string_packed", "h3_indices"), &H3Core::h3s_to_string_packed);
    ClassDB::bind_method(D_METHOD("is_valid_cell", "h3_index"), &H3Core::is_valid_cell);
    ClassDB::bind_method(D_METHOD("is_res_class_iii", "h3_index"), &H3Core::is_res_class_iii);
    ClassDB::bind_method(D_METHOD("is_pentagon", "h3_index"), &H3Core::is_pentagon);
    ClassDB::bind_method(D_METHOD("get_icosahedron_faces", "h3_index"), &H3Core::get_icosahedron_faces);
    ClassDB::bind_method(D_METHOD("get_icosahedron_faces_packed", "h3_index"), &H3Core::get_icosahedron_faces_packed);
    ClassDB::bind_method(D_METHOD("max_face_count", "h3_index"), &H3Core::max_face_count);

    // Grid traversal functions
    ClassDB::bind_method(D_METHOD("grid_distance", "origin", "h3"), &H3Core::grid_distance);
    ClassDB::bind_method(D_METHOD("grid_ring_unsafe", "h3_index", "k"), &H3Core::grid_ring_unsafe);
    ClassDB::bind_method(D_METHOD("grid_ring_unsafe_packed", "h3_index", "k"), &H3Core::grid_ring_unsafe_packed);
    ClassDB::bind_method(D_METHOD("grid_disk", "h3_index", "k"), &H3Core::grid_disk);
    ClassDB::bind_method(D_METHOD("grid_disk_packed", "h3_index", "k"), &H3Core::grid_disk_packed);
    ClassDB::bind_method(D_METHOD("max_grid_disk_size", "k"), &H3Core::max_grid_disk_size);
    ClassDB::bind_method(D_METHOD("grid_disk_distances", "h3_index", "k"), &H3Core::grid_disk_distances);
    ClassDB::bind_method(D_METHOD("grid_disk_unsafe", "h3_index", "k"), &H3Core::grid_disk_unsafe);
    ClassDB::bind_method(D_METHOD("grid_disk_unsafe_packed", "h3_index", "k"), &H3Core::grid_disk_unsafe_packed);
    ClassDB::bind_method(D_METHOD("grid_disk_distances_unsafe", "h3_index", "k"), &H3Core::grid_disk_distances_unsafe);
    // ClassDB::bind_method(D_METHOD("grid_disks_unsafe", "h3_set", "k"), &H3Core::grid_disks_unsafe); // not provided
    ClassDB::bind_method(D_METHOD("grid_path_cells", "start", "end"), &H3Core::grid_path_cells);
    ClassDB::bind_method(D_METHOD("grid_path_cells_packed", "start", "end"), &H3Core::grid_path_cells_packed);
    ClassDB::bind_method(D_METHOD("grid_path_cells_size", "start", "end"), &H3Core::grid_path_cells_size);
    ClassDB::bind_method(D_METHOD("cell_to_local_ij", "origin", "h3_index"), &H3Core::cell_to_local_ij);
    ClassDB::bind_method(D_METHOD("local_ij_to_cell", "origin", "ij"), &H3Core::local_ij_to_cell);

    // Hierarchical grid functions
    ClassDB::bind_method(D_METHOD("cell_to_parent", "h3_index", "parent_res"), &H3Core::cell_to_parent);
    ClassDB::bind_method(D_METHOD("cell_to_children", "h3_index", "child_res"), &H3Core::cell_to_children);
    ClassDB::bind_method(D_METHOD("cell_to_children_packed", "h3_index", "child_res"), &H3Core::cell_to_children_packed);
    ClassDB::bind_method(D_METHOD("cell_to_child_size", "h3_index", "child_res"), &H3Core::cell_to_child_size);
    ClassDB::bind_method(D_METHOD("cell_to_center_child", "h3_index", "child_res"), &H3Core::cell_to_center_child);
    ClassDB::bind_method(D_METHOD("cell_to_child_pos", "child", "parent_res"), &H3Core::cell_to_child_pos);
    ClassDB::bind_method(D_METHOD("child_pos_to_cell", "child_pos", "parent", "child_res"), &H3Core::child_pos_to_cell);
    ClassDB::bind_method(D_METHOD("compact_cells", "h3_set"), &H3Core::compact_cells);
    ClassDB::bind_method(D_METHOD("compact_cells_packed", "h3_set"), &H3Core::compact_cells_packed);
    ClassDB::bind_method(D_METHOD("uncompact_cells", "compacted_set", "res"), &H3Core::uncompact_cells);
    ClassDB::bind_method(D_METHOD("uncompact_cells_packed", "compacted_set", "res"), &H3Core::uncompact_cells_packed);
    ClassDB::bind_method(D_METHOD("uncompact_cells_size", "compacted_set", "res"), &H3Core::uncompact_cells_size);
    ClassDB::bind_method(D_METHOD("uncompact_cells_size_packed", "compacted_set", "res"), &H3Core::uncompact_cells_size_packed);

    // Region functions
    ClassDB::bind_method(D_METHOD("polygon_to_cells", "polygon", "res", "flags", "is_geojson"), &H3Core::polygon_to_cells, DEFVAL(false));
    ClassDB::bind_method(D_METHOD("polygon_to_cells_packed", "polygon", "res", "flags", "is_geojson"), &H3Core::polygon_to_cells_packed, DEFVAL(false));
    ClassDB::bind_method(D_METHOD("max_polygon_to_cells_size", "polygon", "res", "flags", "is_geojson"), &H3Core::max_polygon_to_cells_size, DEFVAL(false));
    ClassDB::bind_method(D_METHOD("polygon_to_cells_experimental", "polygon", "res", "flags", "is_geojson"), &H3Core::polygon_to_cells_experimental, DEFVAL(false));
    ClassDB::bind_method(D_METHOD("polygon_to_cells_experimental_packed", "polygon", "res", "flags", "is_geojson"), &H3Core::polygon_to_cells_experimental_packed, DEFVAL(false));
    ClassDB::bind_method(D_METHOD("max_polygon_to_cells_experimental_size", "polygon", "res", "flags", "is_geojson"), &H3Core::max_polygon_to_cells_experimental_size, DEFVAL(false));
    ClassDB::bind_method(D_METHOD("cells_to_linked_multi_polygon", "h3_set", "is_geojson"), &H3Core::cells_to_linked_multi_polygon, DEFVAL(false));

    // Directed edge functions
    ClassDB::bind_method(D_METHOD("are_neighbor_cells", "origin", "destination"), &H3Core::are_neighbor_cells);
    ClassDB::bind_method(D_METHOD("cells_to_directed_edge", "origin", "destination"), &H3Core::cells_to_directed_edge);
    ClassDB::bind_method(D_METHOD("is_valid_directed_edge", "edge"), &H3Core::is_valid_directed_edge);
    ClassDB::bind_method(D_METHOD("get_directed_edge_origin", "edge"), &H3Core::get_directed_edge_origin);
    ClassDB::bind_method(D_METHOD("get_directed_edge_destination", "edge"), &H3Core::get_directed_edge_destination);
    ClassDB::bind_method(D_METHOD("directed_edge_to_cells", "edge"), &H3Core::directed_edge_to_cells);
    ClassDB::bind_method(D_METHOD("directed_edge_to_cells_packed", "edge"), &H3Core::directed_edge_to_cells_packed);
    ClassDB::bind_method(D_METHOD("origin_to_directed_edges", "origin"), &H3Core::origin_to_directed_edges);
    ClassDB::bind_method(D_METHOD("origin_to_directed_edges_packed", "origin"), &H3Core::origin_to_directed_edges_packed);
    ClassDB::bind_method(D_METHOD("directed_edge_to_boundary", "edge"), &H3Core::directed_edge_to_boundary);
    ClassDB::bind_method(D_METHOD("directed_edge_to_boundary_packed", "edge"), &H3Core::directed_edge_to_boundary_packed);

    // Vertex functions
    ClassDB::bind_method(D_METHOD("cell_to_vertex", "origin", "vertex_num"), &H3Core::cell_to_vertex);
    ClassDB::bind_method(D_METHOD("cell_to_vertices", "origin"), &H3Core::cell_to_vertices);
    ClassDB::bind_method(D_METHOD("cell_to_vertices_packed", "origin"), &H3Core::cell_to_vertices_packed);
    ClassDB::bind_method(D_METHOD("vertex_to_lat_lng", "vertex"), &H3Core::vertex_to_lat_lng);
    ClassDB::bind_method(D_METHOD("is_valid_vertex", "vertex"), &H3Core::is_valid_vertex);

    // Miscellaneous H3 functions
    ClassDB::bind_method(D_METHOD("degs_to_rads", "degrees"), &H3Core::degs_to_rads);
    ClassDB::bind_method(D_METHOD("rads_to_degs", "radians"), &H3Core::rads_to_degs);
    ClassDB::bind_method(D_METHOD("get_hexagon_area_avg_km2", "res"), &H3Core::get_hexagon_area_avg_km2);
    ClassDB::bind_method(D_METHOD("get_hexagon_area_avg_m2", "res"), &H3Core::get_hexagon_area_avg_m2);
    ClassDB::bind_method(D_METHOD("cell_area_rads2", "h3_index"), &H3Core::cell_area_rads2);
    ClassDB::bind_method(D_METHOD("cell_area_km2", "h3_index"), &H3Core::cell_area_km2);
    ClassDB::bind_method(D_METHOD("cell_area_m2", "h3_index"), &H3Core::cell_area_m2);
    ClassDB::bind_method(D_METHOD("get_hexagon_edge_length_avg_km", "res"), &H3Core::get_hexagon_edge_length_avg_km);
    ClassDB::bind_method(D_METHOD("get_hexagon_edge_length_avg_m", "res"), &H3Core::get_hexagon_edge_length_avg_m);
    ClassDB::bind_method(D_METHOD("edge_length_km", "edge"), &H3Core::edge_length_km);
    ClassDB::bind_method(D_METHOD("edge_length_m", "edge"), &H3Core::edge_length_m);
    ClassDB::bind_method(D_METHOD("edge_length_rads", "edge"), &H3Core::edge_length_rads);
    ClassDB::bind_method(D_METHOD("get_num_cells", "res"), &H3Core::get_num_cells);
    ClassDB::bind_method(D_METHOD("get_res0_cells"), &H3Core::get_res0_cells);
    ClassDB::bind_method(D_METHOD("get_res0_cells_packed"), &H3Core::get_res0_cells_packed);
    ClassDB::bind_method(D_METHOD("res0_cell_count"), &H3Core::res0_cell_count);
    ClassDB::bind_method(D_METHOD("get_pentagons", "res"), &H3Core::get_pentagons);
    ClassDB::bind_method(D_METHOD("get_pentagons_packed", "res"), &H3Core::get_pentagons_packed);
    ClassDB::bind_method(D_METHOD("pentagon_count"), &H3Core::pentagon_count);
    ClassDB::bind_method(D_METHOD("great_circle_distance_km", "point1", "point2"), &H3Core::great_circle_distance_km);
    ClassDB::bind_method(D_METHOD("great_circle_distance_m", "point1", "point2"), &H3Core::great_circle_distance_m);
    ClassDB::bind_method(D_METHOD("great_circle_distance_rads", "point1", "point2"), &H3Core::great_circle_distance_rads);

    // Godot functions
    ClassDB::bind_method(D_METHOD("h3_to_vector4", "h3_index"), &H3Core::h3_to_vector4);
    ClassDB::bind_method(D_METHOD("vector4_to_h3", "vec"), &H3Core::vector4_to_h3);
}

// Indexing Functions
H3Index H3Core::lat_lng_to_cell(double lat, double lon, int res) const {
    LatLng coord = {H3_EXPORT(degsToRads)(lat), H3_EXPORT(degsToRads)(lon)};
    H3Index h3Index = 0;
    return handle_h3_error(H3_EXPORT(latLngToCell)(&coord, res, &h3Index)) ? h3Index : H3_NULL;
}

Vector2 H3Core::cell_to_lat_lng(H3Index h3_index) const {
    LatLng coord;
    return handle_h3_error(H3_EXPORT(cellToLatLng)(h3_index, &coord)) ? Vector2(H3_EXPORT(radsToDegs)(coord.lat), H3_EXPORT(radsToDegs)(coord.lng)) : Vector2();
}

Array H3Core::cell_to_boundary(H3Index h3_index) const {
    CellBoundary boundary;
    if (!handle_h3_error(H3_EXPORT(cellToBoundary)(h3_index, &boundary))) {
        return Array();
    }

    Array result;
    for (int i = 0; i < boundary.numVerts; i++) {
        result.append(Vector2(H3_EXPORT(radsToDegs)(boundary.verts[i].lat), H3_EXPORT(radsToDegs)(boundary.verts[i].lng)));
    }
    return result;
}

PackedVector2Array H3Core::cell_to_boundary_packed(H3Index h3_index) const {
    CellBoundary boundary;
    if (!handle_h3_error(H3_EXPORT(cellToBoundary)(h3_index, &boundary))) {
        return PackedVector2Array();
    }

    PackedVector2Array result;
    for (int i = 0; i < boundary.numVerts; i++) {
        result.append(Vector2(H3_EXPORT(radsToDegs)(boundary.verts[i].lat), H3_EXPORT(radsToDegs)(boundary.verts[i].lng)));
    }
    return result;
}

// Index inspection functions
int H3Core::get_resolution(H3Index h3_index) const {
    return H3_EXPORT(getResolution)(h3_index);
}

int H3Core::get_base_cell_number(H3Index h3_index) const {
    return H3_EXPORT(getBaseCellNumber)(h3_index);
}

H3Index H3Core::string_to_h3(const String& h3_str) const {
    H3Index h3Index;
    return handle_h3_error(H3_EXPORT(stringToH3)(h3_str.utf8().get_data(), &h3Index)) ? h3Index : 0;
}

String H3Core::h3_to_string(const H3Index& h3_index) const {
    char str[17];
    return handle_h3_error(H3_EXPORT(h3ToString)(h3_index, str, sizeof(str))) ? String(str) : String();
}

Array H3Core::strings_to_h3(const Array& strs) const {
    Array result;
    result.resize(strs.size());
    for (int i = 0; i < strs.size(); i++) {
        result[i] = string_to_h3(strs[i]);
    }
    return result;
}

Array H3Core::h3s_to_string(const Array& h3_indices) const {
    Array result;
    result.resize(h3_indices.size());
    for (int i = 0; i < h3_indices.size(); i++) {
        result[i] = h3_to_string(h3_indices[i]);
    }
    return result;
}

PackedInt64Array H3Core::strings_to_h3_packed(const PackedStringArray& strs) const {
    PackedInt64Array result;
    result.resize(strs.size());
    for (int i = 0; i < strs.size(); i++) {
        result[i] = string_to_h3(strs[i]);
    }
    return result;
}

PackedStringArray H3Core::h3s_to_string_packed(const PackedInt64Array& h3_indices) const {
    PackedStringArray result;
    result.resize(h3_indices.size());
    for (int i = 0; i < h3_indices.size(); i++) {
        result[i] = h3_to_string(h3_indices[i]);
    }
    return result;
}

bool H3Core::is_valid_cell(H3Index h3_index) const {
    return H3_EXPORT(isValidCell)(h3_index);
}

bool H3Core::is_res_class_iii(H3Index h3_index) const {
    return H3_EXPORT(isResClassIII)(h3_index);
}

bool H3Core::is_pentagon(H3Index h3_index) const {
    return H3_EXPORT(isPentagon)(h3_index);
}

Array H3Core::get_icosahedron_faces(H3Index h3_index) const {
    int maxFaces;
    if (!handle_h3_error(H3_EXPORT(maxFaceCount)(h3_index, &maxFaces))) {
        return Array();
    }

    if (maxFaces <= 0) {
        return Array();
    }

    std::vector<int> faceNums(maxFaces);
    if (!handle_h3_error(H3_EXPORT(getIcosahedronFaces)(h3_index, faceNums.data()))) {
        return Array();
    }

    int count = 0;
    for (size_t i = 0; i < faceNums.size(); i++) {
        if (faceNums[i] >= 0) {
            count++;
        }
    }
    Array result;
    result.resize(count);
    for (size_t i = 0; i < faceNums.size(); i++) {
        if (faceNums[i] >= 0) {
            result[i] = faceNums[i];
        }
    }
    return result;
}

PackedInt64Array H3Core::get_icosahedron_faces_packed(H3Index h3_index) const {
    int maxFaces;
    if (!handle_h3_error(H3_EXPORT(maxFaceCount)(h3_index, &maxFaces))) {
        return PackedInt64Array();
    }

    if (maxFaces <= 0) {
        return PackedInt64Array();
    }

    std::vector<int> faceNums(maxFaces);
    if (!handle_h3_error(H3_EXPORT(getIcosahedronFaces)(h3_index, faceNums.data()))) {
        return PackedInt64Array();
    }

    int count = 0;
    for (size_t i = 0; i < faceNums.size(); i++) {
        if (faceNums[i] >= 0) {
            count++;
        }
    }

    PackedInt64Array result;
    result.resize(count);
    for (size_t i = 0; i < faceNums.size(); i++) {
        if (faceNums[i] >= 0) {
            result[i] = faceNums[i];
        }
    }
    return result;
}

int H3Core::max_face_count(H3Index h3_index) const {
    int maxFaces;
    if (!handle_h3_error(H3_EXPORT(maxFaceCount)(h3_index, &maxFaces))) {
        return 0;
    }
    if (maxFaces <= 0) {
        return 0;
    }

    std::vector<int> faceNums(maxFaces);
    if (!handle_h3_error(H3_EXPORT(getIcosahedronFaces)(h3_index, faceNums.data()))) {
        return 0;
    }

    int count = 0;
    for (size_t i = 0; i < faceNums.size(); i++) {
        if (faceNums[i] >= 0) {
            count++;
        }
    }
    return count;
}

// Grid traversal functions
int64_t H3Core::grid_distance(H3Index origin, H3Index h3) const {
    int64_t distance;
    if (!handle_h3_error(H3_EXPORT(gridDistance)(origin, h3, &distance))) {
        return 0;
    }
    return distance;
}

Array H3Core::grid_ring_unsafe(H3Index h3_index, int k) const {
    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxGridDiskSize)(k, &maxSize))) {
        return Array();
    }

    if (maxSize <= 0) {
        return Array();
    }

    std::vector<H3Index> indices(maxSize);
    if (!handle_h3_error(H3_EXPORT(gridRingUnsafe)(h3_index, k, indices.data()))) {
        return Array();
    }

    Array result;
    for (const auto& index : indices) {
        if (index != 0) {
            result.append(index);
        }
    }
    return result;
}

PackedInt64Array H3Core::grid_ring_unsafe_packed(H3Index h3_index, int k) const {
    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxGridDiskSize)(k, &maxSize))) {
        return PackedInt64Array();
    }

    if (maxSize <= 0) {
        return PackedInt64Array();
    }

    std::vector<H3Index> indices(maxSize);
    if (!handle_h3_error(H3_EXPORT(gridRingUnsafe)(h3_index, k, indices.data()))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    int count = 0;
    result.resize(indices.size());
    for (size_t i = 0; i < indices.size(); i++) {
        if (indices[i] != 0) {
            result[i] = indices[i];
            count++;
        }
    }

    result.resize(count);
    return result;
}

Array H3Core::grid_disk(H3Index h3_index, int k) const {
    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxGridDiskSize)(k, &maxSize))) {
        return Array();
    }

    if (maxSize <= 0) {
        return Array();
    }

    std::vector<H3Index> indices(maxSize);
    if (!handle_h3_error(H3_EXPORT(gridDisk)(h3_index, k, indices.data()))) {
        return Array();
    }

    Array result;
    for (const auto& index : indices) {
        if (index != 0) {
            result.append(index);
        }
    }
    return result;
}

PackedInt64Array H3Core::grid_disk_packed(H3Index h3_index, int k) const {
    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxGridDiskSize)(k, &maxSize))) {
        return PackedInt64Array();
    }

    if (maxSize <= 0) {
        return PackedInt64Array();
    }

    std::vector<H3Index> indices(maxSize);
    if (!handle_h3_error(H3_EXPORT(gridDisk)(h3_index, k, indices.data()))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(indices.size());
    int count = 0;
    for (size_t i = 0; i < indices.size(); i++) {
        if (indices[i] != 0) {
            result[i] = indices[i];
            count++;
        }
    }
    result.resize(count);
    return result;
}

int64_t H3Core::max_grid_disk_size(int k) const {
    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxGridDiskSize)(k, &maxSize))) {
        return 0;
    }
    return maxSize;
}

Array H3Core::grid_disk_distances(H3Index h3_index, int k) const {
    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxGridDiskSize)(k, &maxSize))) {
        return Array();
    }

    if (maxSize <= 0) {
        return Array();
    }

    std::vector<H3Index> indices(maxSize);
    std::vector<int> distances(maxSize);
    if (!handle_h3_error(H3_EXPORT(gridDiskDistances)(h3_index, k, indices.data(), distances.data()))) {
        return Array();
    }

    Array result;
    result.resize(k + 1);

    for (int i = 0; i <= k; i++) {
        result[i] = Array();
    }

    for (size_t i = 0; i < maxSize; i++) {
        if (indices[i] != 0) {
            Array current = result[distances[i]];
            current.append(indices[i]);
            result[distances[i]] = current;
        }
    }

    return result;
}

Array H3Core::grid_disk_unsafe(H3Index h3_index, int k) const {
    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxGridDiskSize)(k, &maxSize))) {
        return Array();
    }

    if (maxSize <= 0) {
        return Array();
    }

    std::vector<H3Index> indices(maxSize);
    if (!handle_h3_error(H3_EXPORT(gridDiskUnsafe)(h3_index, k, indices.data()))) {
        return Array();
    }

    Array result;
    for (const auto& index : indices) {
        if (index != 0) {
            result.append(index);
        }
    }
    return result;
}

PackedInt64Array H3Core::grid_disk_unsafe_packed(H3Index h3_index, int k) const {
    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxGridDiskSize)(k, &maxSize))) {
        return PackedInt64Array();
    }

    if (maxSize <= 0) {
        return PackedInt64Array();
    }

    std::vector<H3Index> indices(maxSize);
    if (!handle_h3_error(H3_EXPORT(gridDiskUnsafe)(h3_index, k, indices.data()))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(indices.size());
    int count = 0;
    for (size_t i = 0; i < indices.size(); i++) {
        if (indices[i] != 0) {
            result[i] = indices[i];
            count++;
        }
    }
    result.resize(count);
    return result;
}

Array H3Core::grid_disk_distances_unsafe(H3Index h3_index, int k) const {
    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxGridDiskSize)(k, &maxSize))) {
        return Array();
    }

    if (maxSize <= 0) {
        return Array();
    }

    std::vector<H3Index> indices(maxSize);
    std::vector<int> distances(maxSize);
    if (!handle_h3_error(H3_EXPORT(gridDiskDistancesUnsafe)(h3_index, k, indices.data(), distances.data()))) {
        return Array();
    }

    Array result;
    result.resize(k + 1);

    for (int i = 0; i <= k; i++) {
        result[i] = Array();
    }

    for (size_t i = 0; i < maxSize; i++) {
        if (indices[i] != 0) {
            Array current = result[distances[i]];
            current.append(indices[i]);
            result[distances[i]] = current;
        }
    }

    return result;
}

// Array H3Core::grid_disks_unsafe(const Array& h3_set, int k) const {
//     std::vector<H3Index> input;
//     for (int i = 0; i < h3_set.size(); i++) {
//         input.push_back(h3_set[i].operator H3Index());
//     }

//     int64_t maxSize;
//     H3Error error = H3_EXPORT(maxGridDiskSize)(k, &maxSize);
//     handle_h3_error(error);
//     maxSize *= input.size();

//     if (maxSize <= 0) {
//         return Array();
//     }

//     std::vector<H3Index> out(maxSize);
//     error = H3_EXPORT(gridDisksUnsafe)(input.data(), input.size(), k, out.data());
//     handle_h3_error(error);

//     Array result;
//     for (const auto& index : out) {
//         if (index != 0) {
//             result.append(index);
//         }
//     }
//     return result;
// }

Array H3Core::grid_path_cells(H3Index start, H3Index end) const {
    int64_t size;
    if (!handle_h3_error(H3_EXPORT(gridPathCellsSize)(start, end, &size))) {
        return Array();
    }

    std::vector<H3Index> path(size);
    if (!handle_h3_error(H3_EXPORT(gridPathCells)(start, end, path.data()))) {
        return Array();
    }

    Array result;
    for (const auto& cell : path) {
        if (cell != 0) {
            result.append(cell);
        }
    }
    return result;
}

PackedInt64Array H3Core::grid_path_cells_packed(H3Index start, H3Index end) const {
    int64_t size;
    if (!handle_h3_error(H3_EXPORT(gridPathCellsSize)(start, end, &size))) {
        return PackedInt64Array();
    }

    std::vector<H3Index> path(size);
    if (!handle_h3_error(H3_EXPORT(gridPathCells)(start, end, path.data()))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(path.size());
    int count = 0;
    for (size_t i = 0; i < path.size(); i++) {
        if (path[i] != 0) {
            result[i] = path[i];
            count++;
        }
    }
    result.resize(count);
    return result;
}

int64_t H3Core::grid_path_cells_size(H3Index start, H3Index end) const {
    int64_t size;
    if (!handle_h3_error(H3_EXPORT(gridPathCellsSize)(start, end, &size))) {
        return 0;
    }
    return size;
}

Vector2i H3Core::cell_to_local_ij(H3Index origin, H3Index h3_index) const {
    CoordIJ ij;
    if (!handle_h3_error(H3_EXPORT(cellToLocalIj)(origin, h3_index, 0, &ij))) {
        return Vector2i();
    }

    Vector2i result;
    result.x = ij.i;
    result.y = ij.j;
    return result;
}

H3Index H3Core::local_ij_to_cell(H3Index origin, const Vector2i& ij) const {
    CoordIJ coord;
    coord.i = ij.x;
    coord.j = ij.y;

    H3Index h3_index;
    if (!handle_h3_error(H3_EXPORT(localIjToCell)(origin, &coord, 0, &h3_index))) {
        return 0;
    }
    return h3_index;
}

// Hierarchical grid functions
H3Index H3Core::cell_to_parent(H3Index h3_index, int parent_res) const {
    H3Index parent;
    if (!handle_h3_error(H3_EXPORT(cellToParent)(h3_index, parent_res, &parent))) {
        return 0;
    }
    return parent;
}

Array H3Core::cell_to_children(H3Index h3_index, int child_res) const {
    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(cellToChildrenSize)(h3_index, child_res, &maxSize))) {
        return Array();
    }

    if (maxSize <= 0) {
        return Array();
    }

    std::vector<H3Index> children(maxSize);
    if (!handle_h3_error(H3_EXPORT(cellToChildren)(h3_index, child_res, children.data()))) {
        return Array();
    }

    Array result;
    for (const auto& child : children) {
        if (child != 0) {
            result.append(child);
        }
    }
    return result;
}

PackedInt64Array H3Core::cell_to_children_packed(H3Index h3_index, int child_res) const {
    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(cellToChildrenSize)(h3_index, child_res, &maxSize))) {
        return PackedInt64Array();
    }

    if (maxSize <= 0) {
        return PackedInt64Array();
    }

    std::vector<H3Index> children(maxSize);
    if (!handle_h3_error(H3_EXPORT(cellToChildren)(h3_index, child_res, children.data()))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(children.size());
    int count = 0;
    for (size_t i = 0; i < children.size(); i++) {
        if (children[i] != 0) {
            result[i] = children[i];
            count++;
        }
    }
    result.resize(count);
    return result;
}

int64_t H3Core::cell_to_child_size(H3Index h3_index, int child_res) const {
    int64_t size;
    if (!handle_h3_error(H3_EXPORT(cellToChildrenSize)(h3_index, child_res, &size))) {
        return 0;
    }
    return size;
}

H3Index H3Core::cell_to_center_child(H3Index h3_index, int child_res) const {
    H3Index center;
    if (!handle_h3_error(H3_EXPORT(cellToCenterChild)(h3_index, child_res, &center))) {
        return 0;
    }
    return center;
}

int64_t H3Core::cell_to_child_pos(H3Index child, int parent_res) const {
    int64_t child_pos;
    if (!handle_h3_error(H3_EXPORT(cellToChildPos)(child, parent_res, &child_pos))) {
        return 0;
    }
    return child_pos;
}

H3Index H3Core::child_pos_to_cell(int64_t child_pos, H3Index parent, int child_res) const {
    H3Index child;
    if (!handle_h3_error(H3_EXPORT(childPosToCell)(child_pos, parent, child_res, &child))) {
        return 0;
    }
    return child;
}

Array H3Core::compact_cells(const Array& h3_set) const {
    std::vector<H3Index> input;
    std::vector<H3Index> compacted;

    for (int i = 0; i < h3_set.size(); i++) {
        input.push_back(h3_set[i].operator H3Index());
    }

    compacted.resize(input.size());
    if (!handle_h3_error(H3_EXPORT(compactCells)(input.data(), compacted.data(), input.size()))) {
        return Array();
    }

    Array result;
    for (const auto& cell : compacted) {
        if (cell != 0) {
            result.append(cell);
        }
    }
    return result;
}

PackedInt64Array H3Core::compact_cells_packed(const PackedInt64Array& h3_set) const {
    std::vector<H3Index> input;
    std::vector<H3Index> compacted;

    for (int i = 0; i < h3_set.size(); i++) {
        input.push_back(h3_set[i]);
    }

    compacted.resize(input.size());
    if (!handle_h3_error(H3_EXPORT(compactCells)(input.data(), compacted.data(), input.size()))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(compacted.size());
    int count = 0;
    for (size_t i = 0; i < compacted.size(); i++) {
        if (compacted[i] != 0) {
            result[i] = compacted[i];
            count++;
        }
    }
    result.resize(count);
    return result;
}

Array H3Core::uncompact_cells(const Array& compacted_set, int res) const {
    std::vector<H3Index> input;

    for (int i = 0; i < compacted_set.size(); i++) {
        input.push_back(compacted_set[i].operator H3Index());
    }

    int64_t size;
    if (!handle_h3_error(H3_EXPORT(uncompactCellsSize)(input.data(), input.size(), res, &size))) {
        return Array();
    }

    std::vector<H3Index> uncompacted(size);
    if (!handle_h3_error(H3_EXPORT(uncompactCells)(input.data(), input.size(), uncompacted.data(), size, res))) {
        return Array();
    }

    Array result;
    for (const auto& cell : uncompacted) {
        if (cell != 0) {
            result.append(cell);
        }
    }
    return result;
}

PackedInt64Array H3Core::uncompact_cells_packed(const PackedInt64Array& compacted_set, int res) const {
    std::vector<H3Index> input;
    std::vector<H3Index> uncompacted;

    for (int i = 0; i < compacted_set.size(); i++) {
        input.push_back(compacted_set[i]);
    }

    int64_t size;
    if (!handle_h3_error(H3_EXPORT(uncompactCellsSize)(input.data(), input.size(), res, &size))) {
        return PackedInt64Array();
    }

    uncompacted.resize(size);
    if (!handle_h3_error(H3_EXPORT(uncompactCells)(input.data(), input.size(), uncompacted.data(), size, res))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(uncompacted.size());
    int count = 0;
    for (size_t i = 0; i < uncompacted.size(); i++) {
        if (uncompacted[i] != 0) {
            result[i] = uncompacted[i];
            count++;
        }
    }
    result.resize(count);
    return result;
}

int64_t H3Core::uncompact_cells_size(const Array& compacted_set, int res) const {
    std::vector<H3Index> input;
    for (int i = 0; i < compacted_set.size(); i++) {
        input.push_back(compacted_set[i].operator H3Index());
    }
    int64_t size;
    if (!handle_h3_error(H3_EXPORT(uncompactCellsSize)(input.data(), input.size(), res, &size))) {
        return 0;
    }
    return size;
}

int64_t H3Core::uncompact_cells_size_packed(const PackedInt64Array& compacted_set, int res) const {
    std::vector<H3Index> input;
    for (int i = 0; i < compacted_set.size(); i++) {
        input.push_back(compacted_set[i]);
    }
    int64_t size;
    if (!handle_h3_error(H3_EXPORT(uncompactCellsSize)(input.data(), input.size(), res, &size))) {
        return 0;
    }
    return size;
}

// Region functions
Array H3Core::polygon_to_cells(const Array& polygon_array, int res, uint32_t flags, bool is_geojson) const {
    GeoPolygon polygon = array_to_geopolygon(polygon_array, is_geojson);

    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxPolygonToCellsSize)(&polygon, res, flags, &maxSize))) {
        return Array();
    }

    if (maxSize <= 0) {
        free_geopolygon(polygon);
        return Array();
    }

    std::vector<H3Index> out(maxSize);
    if (!handle_h3_error(H3_EXPORT(polygonToCells)(&polygon, res, flags, out.data()))) {
        return Array();
    }

    Array result;
    for (const auto& cell : out) {
        if (cell != 0) {
            result.append(cell);
        }
    }

    free_geopolygon(polygon);
    return result;
}

PackedInt64Array H3Core::polygon_to_cells_packed(const Array& polygon_array, int res, uint32_t flags, bool is_geojson) const {
    GeoPolygon polygon = array_to_geopolygon(polygon_array, is_geojson);

    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxPolygonToCellsSize)(&polygon, res, flags, &maxSize))) {
        return PackedInt64Array();
    }

    if (maxSize <= 0) {
        free_geopolygon(polygon);
        return Array();
    }

    std::vector<H3Index> out(maxSize);
    if (!handle_h3_error(H3_EXPORT(polygonToCells)(&polygon, res, flags, out.data()))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(out.size());
    int count = 0;
    for (size_t i = 0; i < out.size(); i++) {
        if (out[i] != 0) {
            result[count] = out[i];
            count++;
        }
    }

    free_geopolygon(polygon);
    result.resize(count);
    return result;
}

int64_t H3Core::max_polygon_to_cells_size(const Array& polygon_array, int res, uint32_t flags, bool is_geojson) const {
    GeoPolygon polygon = array_to_geopolygon(polygon_array, is_geojson);

    int64_t out;
    if (!handle_h3_error(H3_EXPORT(maxPolygonToCellsSize)(&polygon, res, flags, &out))) {
        return 0;
    }

    free_geopolygon(polygon);
    return out;
}

Array H3Core::polygon_to_cells_experimental(const Array& polygon_array, int res, uint32_t flags, bool is_geojson) const {
    GeoPolygon polygon = array_to_geopolygon(polygon_array, is_geojson);

    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxPolygonToCellsSizeExperimental)(&polygon, res, flags, &maxSize))) {
        return Array();
    }

    if (maxSize <= 0) {
        free_geopolygon(polygon);
        return Array();
    }

    std::vector<H3Index> out(maxSize);
    if (!handle_h3_error(H3_EXPORT(polygonToCellsExperimental)(&polygon, res, flags, maxSize, out.data()))) {
        return Array();
    }

    Array result;
    for (const auto& cell : out) {
        if (cell != 0) {
            result.append(cell);
        }
    }

    free_geopolygon(polygon);
    return result;
}

PackedInt64Array H3Core::polygon_to_cells_experimental_packed(const Array& polygon_array, int res, uint32_t flags, bool is_geojson) const {
    GeoPolygon polygon = array_to_geopolygon(polygon_array, is_geojson);

    int64_t maxSize;
    if (!handle_h3_error(H3_EXPORT(maxPolygonToCellsSizeExperimental)(&polygon, res, flags, &maxSize))) {
        return PackedInt64Array();
    }

    if (maxSize <= 0) {
        free_geopolygon(polygon);
        return PackedInt64Array();
    }

    std::vector<H3Index> out(maxSize);
    if (!handle_h3_error(H3_EXPORT(polygonToCellsExperimental)(&polygon, res, flags, maxSize, out.data()))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(out.size());
    int count = 0;
    for (size_t i = 0; i < out.size(); i++) {
        if (out[i] != 0) {
            result[i] = out[i];
            count++;
        }
    }

    result.resize(count);
    return result;
}

int64_t H3Core::max_polygon_to_cells_experimental_size(const Array& polygon_array, int res, uint32_t flags, bool is_geojson) const {
    GeoPolygon polygon = array_to_geopolygon(polygon_array, is_geojson);

    int64_t out;
    if (!handle_h3_error(H3_EXPORT(maxPolygonToCellsSizeExperimental)(&polygon, res, flags, &out))) {
        return 0;
    }

    free_geopolygon(polygon);
    return out;
}

Array H3Core::cells_to_linked_multi_polygon(const Array& h3_set, bool is_geojson) const {
    std::vector<H3Index> indices;
    for (int i = 0; i < h3_set.size(); i++) {
        indices.push_back(h3_set[i].operator H3Index());
    }

    LinkedGeoPolygon polygon;
    if (!handle_h3_error(H3_EXPORT(cellsToLinkedMultiPolygon)(indices.data(), indices.size(), &polygon))) {
        return Array();
    }

    Array result = linked_geo_polygon_to_dict(&polygon, is_geojson);
    H3_EXPORT(destroyLinkedMultiPolygon)(&polygon);

    return result;
}

// Directed edge functions
bool H3Core::are_neighbor_cells(H3Index origin, H3Index destination) const {
    int result;
    return handle_h3_error(H3_EXPORT(areNeighborCells)(origin, destination, &result)) ? result != 0 : false;
}

H3Index H3Core::cells_to_directed_edge(H3Index origin, H3Index destination) const {
    H3Index edge;
    return handle_h3_error(H3_EXPORT(cellsToDirectedEdge)(origin, destination, &edge)) ? edge : 0;
}

bool H3Core::is_valid_directed_edge(const H3Index edge) const {
    return H3_EXPORT(isValidDirectedEdge)(edge) != 0;
}

H3Index H3Core::get_directed_edge_origin(const H3Index edge) const {
    H3Index origin;
    return handle_h3_error(H3_EXPORT(getDirectedEdgeOrigin)(edge, &origin)) ? origin : 0;
}

H3Index H3Core::get_directed_edge_destination(const H3Index edge) const {
    H3Index destination;
    return handle_h3_error(H3_EXPORT(getDirectedEdgeDestination)(edge, &destination)) ? destination : 0;
}

Array H3Core::directed_edge_to_cells(const H3Index edge) const {
    H3Index originDestination[2];
    if (!handle_h3_error(H3_EXPORT(directedEdgeToCells)(edge, originDestination))) {
        return Array();
    }

    Array result;
    result.append(originDestination[0]);
    result.append(originDestination[1]);
    return result;
}

PackedInt64Array H3Core::directed_edge_to_cells_packed(const H3Index edge) const {
    H3Index originDestination[2];
    if (!handle_h3_error(H3_EXPORT(directedEdgeToCells)(edge, originDestination))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(2);
    result[0] = originDestination[0];
    result[1] = originDestination[1];
    return result;
}

Array H3Core::origin_to_directed_edges(H3Index origin) const {
    H3Index edges[6];
    if (!handle_h3_error(H3_EXPORT(originToDirectedEdges)(origin, edges))) {
        return Array();
    }

    Array result;
    for (int i = 0; i < 6; i++) {
        if (edges[i] != 0) {
            result.append(edges[i]);
        }
    }
    return result;
}

PackedInt64Array H3Core::origin_to_directed_edges_packed(H3Index origin) const {
    H3Index edges[6];
    if (!handle_h3_error(H3_EXPORT(originToDirectedEdges)(origin, edges))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(6);
    int count = 0;
    for (int i = 0; i < 6; i++) {
        if (edges[i] != 0) {
            result[count] = edges[i];
            count++;
        }
    }

    result.resize(count);
    return result;
}

Array H3Core::directed_edge_to_boundary(const H3Index edge) const {
    CellBoundary boundary;
    if (!handle_h3_error(H3_EXPORT(directedEdgeToBoundary)(edge, &boundary))) {
        return Array();
    }

    Array result;

    for (int i = 0; i < boundary.numVerts; i++) {
        Vector2 vertex;
        vertex.x = H3_EXPORT(radsToDegs)(boundary.verts[i].lat);
        vertex.y = H3_EXPORT(radsToDegs)(boundary.verts[i].lng);
        result.append(vertex);
    }

    return result;
}

PackedVector2Array H3Core::directed_edge_to_boundary_packed(const H3Index edge) const {
    CellBoundary boundary;
    if (!handle_h3_error(H3_EXPORT(directedEdgeToBoundary)(edge, &boundary))) {
        return PackedVector2Array();
    }

    PackedVector2Array result;
    result.resize(boundary.numVerts);
    for (int i = 0; i < boundary.numVerts; i++) {
        result[i] = Vector2(H3_EXPORT(radsToDegs)(boundary.verts[i].lat), H3_EXPORT(radsToDegs)(boundary.verts[i].lng));
    }
    return result;
}

// Vertex functions
H3Index H3Core::cell_to_vertex(H3Index origin, int vertex_num) const {
    H3Index vertex;
    return handle_h3_error(H3_EXPORT(cellToVertex)(origin, vertex_num, &vertex)) ? vertex : 0;
}

Array H3Core::cell_to_vertices(H3Index origin) const {
    H3Index vertices[6];
    if (!handle_h3_error(H3_EXPORT(cellToVertexes)(origin, vertices))) {
        return Array();
    }

    Array result;
    int count = H3_EXPORT(isPentagon)(origin) == 0 ? 6 : 5;
    for (int i = 0; i < count; i++) {
        result.append(vertices[i]);
    }
    return result;
}

PackedInt64Array H3Core::cell_to_vertices_packed(H3Index origin) const {
    H3Index vertices[6];
    if (!handle_h3_error(H3_EXPORT(cellToVertexes)(origin, vertices))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    int count = H3_EXPORT(isPentagon)(origin) == 0 ? 6 : 5;
    result.resize(count);
    for (int i = 0; i < count; i++) {
        result[i] = vertices[i];
    }

    return result;
}

Vector2 H3Core::vertex_to_lat_lng(const H3Index vertex) const {
    LatLng point;
    if (!handle_h3_error(H3_EXPORT(vertexToLatLng)(vertex, &point))) {
        return Vector2();
    }

    Vector2 result;
    result.x = H3_EXPORT(radsToDegs)(point.lat);
    result.y = H3_EXPORT(radsToDegs)(point.lng);
    return result;
}

bool H3Core::is_valid_vertex(const H3Index vertex) const {
    return H3_EXPORT(isValidVertex)(vertex) != 0;
}

// Miscellaneous H3 functions
double H3Core::degs_to_rads(double degrees) const {
    return H3_EXPORT(degsToRads)(degrees);
}

double H3Core::rads_to_degs(double radians) const {
    return H3_EXPORT(radsToDegs)(radians);
}

double H3Core::get_hexagon_area_avg_km2(int res) const {
    double area;
    return handle_h3_error(H3_EXPORT(getHexagonAreaAvgKm2)(res, &area)) ? area : 0;
}

double H3Core::get_hexagon_area_avg_m2(int res) const {
    double area;
    return handle_h3_error(H3_EXPORT(getHexagonAreaAvgM2)(res, &area)) ? area : 0;
}

double H3Core::cell_area_rads2(H3Index h3_index) const {
    double area;
    return handle_h3_error(H3_EXPORT(cellAreaRads2)(h3_index, &area)) ? area : 0;
}

double H3Core::cell_area_km2(H3Index h3_index) const {
    double area;
    return handle_h3_error(H3_EXPORT(cellAreaKm2)(h3_index, &area)) ? area : 0;
}

double H3Core::cell_area_m2(H3Index h3_index) const {
    double area;
    return handle_h3_error(H3_EXPORT(cellAreaM2)(h3_index, &area)) ? area : 0;
}

double H3Core::get_hexagon_edge_length_avg_km(int res) const {
    double length;
    return handle_h3_error(H3_EXPORT(getHexagonEdgeLengthAvgKm)(res, &length)) ? length : 0;
}

double H3Core::get_hexagon_edge_length_avg_m(int res) const {
    double length;
    return handle_h3_error(H3_EXPORT(getHexagonEdgeLengthAvgM)(res, &length)) ? length : 0;
}

double H3Core::edge_length_km(const H3Index edge) const {
    double length;
    return handle_h3_error(H3_EXPORT(edgeLengthKm)(edge, &length)) ? length : 0;
}

double H3Core::edge_length_m(const H3Index edge) const {
    double length;
    return handle_h3_error(H3_EXPORT(edgeLengthM)(edge, &length)) ? length : 0;
}

double H3Core::edge_length_rads(const H3Index edge) const {
    double length;
    return handle_h3_error(H3_EXPORT(edgeLengthRads)(edge, &length)) ? length : 0;
}

int64_t H3Core::get_num_cells(int res) const {
    int64_t num;
    return handle_h3_error(H3_EXPORT(getNumCells)(res, &num)) ? num : 0;
}

Array H3Core::get_res0_cells() const {
    std::vector<H3Index> cells(H3_EXPORT(res0CellCount)());
    if (!handle_h3_error(H3_EXPORT(getRes0Cells)(cells.data()))) {
        return Array();
    }

    Array result;
    for (const auto& cell : cells) {
        result.append(cell);
    }
    return result;
}

PackedInt64Array H3Core::get_res0_cells_packed() const {
    int count = H3_EXPORT(res0CellCount)();
    std::vector<H3Index> cells(count);
    if (!handle_h3_error(H3_EXPORT(getRes0Cells)(cells.data()))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(count);
    for (int i = 0; i < count; i++) {
        result[i] = cells[i];
    }
    return result;
}

int H3Core::res0_cell_count() const {
    return H3_EXPORT(res0CellCount)();
}

Array H3Core::get_pentagons(int res) const {
    std::vector<H3Index> pentagons(12);
    if (!handle_h3_error(H3_EXPORT(getPentagons)(res, pentagons.data()))) {
        return Array();
    }

    Array result;
    for (const auto& pentagon : pentagons) {
        if (pentagon != 0) {
            result.append(pentagon);
        }
    }
    return result;
}

PackedInt64Array H3Core::get_pentagons_packed(int res) const {
    std::vector<H3Index> pentagons(12);
    if (!handle_h3_error(H3_EXPORT(getPentagons)(res, pentagons.data()))) {
        return PackedInt64Array();
    }

    PackedInt64Array result;
    result.resize(12);
    for (int i = 0; i < 12; i++) {
        result[i] = pentagons[i];
    }
    return result;
}

int H3Core::pentagon_count() const {
    return H3_EXPORT(pentagonCount)();
}

double H3Core::great_circle_distance_km(const Vector2 point1, const Vector2 point2) const {
    LatLng p1, p2;
    p1.lat = H3_EXPORT(degsToRads)(static_cast<double>(point1.x));
    p1.lng = H3_EXPORT(degsToRads)(static_cast<double>(point1.y));
    p2.lat = H3_EXPORT(degsToRads)(static_cast<double>(point2.x));
    p2.lng = H3_EXPORT(degsToRads)(static_cast<double>(point2.y));

    return H3_EXPORT(greatCircleDistanceKm)(&p1, &p2);
}

double H3Core::great_circle_distance_m(const Vector2 point1, const Vector2 point2) const {
    LatLng p1, p2;
    p1.lat = H3_EXPORT(degsToRads)(static_cast<double>(point1.x));
    p1.lng = H3_EXPORT(degsToRads)(static_cast<double>(point1.y));
    p2.lat = H3_EXPORT(degsToRads)(static_cast<double>(point2.x));
    p2.lng = H3_EXPORT(degsToRads)(static_cast<double>(point2.y));

    return H3_EXPORT(greatCircleDistanceM)(&p1, &p2);
}

double H3Core::great_circle_distance_rads(const Vector2 point1, const Vector2 point2) const {
    LatLng p1, p2;
    p1.lat = H3_EXPORT(degsToRads)(static_cast<double>(point1.x));
    p1.lng = H3_EXPORT(degsToRads)(static_cast<double>(point1.y));
    p2.lat = H3_EXPORT(degsToRads)(static_cast<double>(point2.x));
    p2.lng = H3_EXPORT(degsToRads)(static_cast<double>(point2.y));

    return H3_EXPORT(greatCircleDistanceRads)(&p1, &p2);
}

// Godot functions
Vector4 H3Core::h3_to_vector4(H3Index h3_index) const {
    char str[17];
    if (!handle_h3_error(H3_EXPORT(h3ToString)(h3_index, str, sizeof(str)))) {
        return Vector4();
    }

    size_t len = strlen(str);
    size_t padding = (16 - len);
    char padded_str[17];
    memset(padded_str, '0', padding);
    strcpy(padded_str + padding, str);

    uint16_t segments[4] = {0, 0, 0, 0};
    for (int i = 0; i < 16; i += 4) {
        char segment[5] = {0};
        strncpy(segment, padded_str + i, 4);
        unsigned long value = strtoul(segment, nullptr, 16);
        segments[i / 4] = static_cast<uint16_t>(value);
    }

    return Vector4(segments[0], segments[1], segments[2], segments[3]);
}

H3Index H3Core::vector4_to_h3(const Vector4& vec) const {
    char str[17];
    snprintf(str, sizeof(str), "%04x%04x%04x%04x", static_cast<uint16_t>(vec.x), static_cast<uint16_t>(vec.y), static_cast<uint16_t>(vec.z), static_cast<uint16_t>(vec.w));

    int start_pos = 0;
    while (start_pos < 15 && str[start_pos] == '0') {
        start_pos++;
    }

    H3Index h3_index;
    if (!handle_h3_error(H3_EXPORT(stringToH3)(str + start_pos, &h3_index))) {
        return 0;
    }

    return h3_index;
}
