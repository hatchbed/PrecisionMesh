#include <deque>
#include <filesystem>
#include <limits>
#include <map>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include <boost/algorithm/string/case_conv.hpp>
#include <CLI/CLI.hpp>
#include <spdlog/fmt/fmt.h>
#include <spdlog/spdlog.h>
#include <tbb/parallel_for.h>
#include <tbb/global_control.h>

// CGAL
#define CGAL_NO_PRECONDITIONS
#define CGAL_NO_ASSERTIONS
#define CGAL_NO_WARNINGS
#include <CGAL/Polygon_mesh_processing/Adaptive_sizing_field.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Surface_mesh.h>

// OCCT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <BRep_Tool.hxx>
#include <BRepGProp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <GProp_GProps.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Message.hxx>
#include <Message_Alert.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <StlAPI_Writer.hxx>
#include <TColStd_SequenceOfAsciiString.hxx>
#include <TDataStd_Name.hxx>
#include <TDF_Label.hxx>
#include <TDF_LabelSequence.hxx>
#include <TDocStd_Document.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopExp_Explorer.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#pragma GCC diagnostic pop

#include <precision_mesh/mesh_util.h>
#include <precision_mesh/ply.h>
#include <precision_mesh/step_util.h>
#include <precision_mesh/stl.h>

namespace PMP = CGAL::Polygon_mesh_processing;

std::unordered_set<std::string> mesh_formats = {".obj", ".off", ".ply", ".stl", ".ts", ".vtp"};

struct Component {
    Component(TDF_Label label, const std::string& name, size_t index,
              const std::string& reference_name, size_t depth) :
        label(label), name(name), index(index), reference_name(reference_name), depth(depth) {}

    TDF_Label label;
    std::string name;
    int32_t index;
    std::string reference_name;
    size_t depth;
    TopoDS_Shape shape;
    double surface_area = 0.0;
    std::vector<std::shared_ptr<Component>> children;
};

std::string getName(const TDF_Label& label) {
    std::string name = "<unknown>";
    Handle(TDataStd_Name) name_attribute;
    if (label.FindAttribute(TDataStd_Name::GetID(), name_attribute)) {
        TCollection_AsciiString utf8String(name_attribute->Get().ToExtString(), Standard_False);
        name = utf8String.ToCString();
    }
    return name;
}

void saveOutput(const std::vector<std::string>& outputs, const std::vector<Mesh>& meshes) {
    for (const auto& output: outputs) {

        std::filesystem::path output_path(output);
        std::string extension = boost::algorithm::to_lower_copy(output_path.extension().string());

        if (extension == ".ply") {
            spdlog::info("saving mesh to: {}", output);
            saveComponentsToPly<Point_traits>(output, meshes);
        }
        else if (extension == ".stl") {
            spdlog::info("saving mesh to: {}", output);
            saveComponentsToStl<Point_traits>(output, meshes);
        }
    }
}

void exploreAssembly(Handle(XCAFDoc_ShapeTool) &assembly, Component& parent, std::map<std::string, size_t>& counts) {
    TDF_LabelSequence children;
    assembly->GetComponents(parent.label, children);

    for (Standard_Integer i = 1; i <= children.Length(); i++) {
        TDF_Label label = children.Value(i);

        std::string name = getName(label);
        std::string reference_name;

        size_t index = ++counts[name];

        TDF_Label reference;
        if (assembly->GetReferredShape(label, reference)) {
            reference_name = getName(reference);
            label = reference;
            spdlog::debug("has reference: {}", reference_name);
        }

        spdlog::debug("    {}{}", std::string(2 * (parent.depth + 1), ' '), name);

        parent.children.push_back(
            std::make_shared<Component>(label, name, index, reference_name, parent.depth + 1));

        exploreAssembly(assembly, *parent.children.back(), counts);
    }
}

std::string normalizeUnit(const std::string& unit) {
    if (unit.empty()) {
        return "mesh units";
    }

    auto lower = boost::algorithm::to_lower_copy(unit);
    if (lower == "millimetre" || lower == "millimetres" || lower == "millimeter" ||
        lower == "millimeters" ||  lower == "mm") {
        return "mm";
    }
    else if (lower == "centimetre" || lower == "centimeters" || lower == "centimeter" ||
             lower == "centimeters" || lower == "cm") {
        return "cm";
    }
    else if (lower == "metre" || lower == "metres" || lower == "meter" || lower == "meters" ||
             lower == "m") {
        return "m";
    }
    else if (lower == "inch" || lower == "inches" || lower == "in") {
        return "in";
    }
    else if (lower == "foot" || lower == "feet" || lower == "ft") {
        return "ft";
    }
    else if (lower == "yard" || lower == "yards" || lower == "yd") {
        return "yd";
    }

    return lower;
}

int main(int argc, char **argv) {

    CLI::App app{"A flexible STEP to mesh conversion tool and general adaptive isotropic remesher.", "precision_mesh"};
    argv = app.ensure_utf8(argv);

    std::string input = "";
    app.add_option("-i,--input", input, "Input file (.obj|.off|.ply|.step|.stl|.ts|.vtp)")
        ->check(CLI::ExistingFile)
        ->required();

    std::string shape_name;
    app.add_option("-n,--name", shape_name, "STEP file shape name.");

    int shape_instance_index = 1;
    app.add_option("--instance", shape_instance_index, "STEP file shape instance index.")
        ->check(CLI::PositiveNumber);

    std::vector<std::string> outputs;
    app.add_option("-o,--output", outputs, "Output file (.obj|.off|.ply|.stl|.ts|.vtp)")->take_all();

    bool list_step_components = false;
    app.add_flag("--list",  list_step_components, "List STEP components.");

    bool raw_step_mesh = false;
    app.add_flag("-r,--raw-step-mesh",  raw_step_mesh, "Generate raw mesh from STEP component without remeshing.");

    int iterations = 0;
    app.add_option("--iterations", iterations, "Iterations")->check(CLI::PositiveNumber);

    bool no_projection = false;
    app.add_flag("--no-projection",  no_projection, "Skip vertex reprojection.");

    bool unfreeze_step_boundaries = false;
    app.add_flag("--unfreeze-step-boundaries", unfreeze_step_boundaries,
                 "Unfreeze STEP surface boundaries when remeshing");

    double crease_angle = std::numeric_limits<double>::quiet_NaN();;
    app.add_option("-a,--crease-angle", crease_angle,
                   "Minimum threshold in degrees of the dihidral angle of edges to freeze when "
                   "remeshing")
        ->check(CLI::PositiveNumber)
        ->check(CLI::Range(0.0, 180.0));

    double max_surface_error_percent = std::numeric_limits<double>::quiet_NaN();
    auto max_surface_error_percent_opt = app.add_option("--max-surface-error-percent",
        max_surface_error_percent,
        "Target maximum surface error when remeshing as percent of sqrt of surface area")
        ->check(CLI::PositiveNumber)
        ->check(CLI::Range(0.0, 100.0));

    double max_boundary_surface_error_percent = std::numeric_limits<double>::quiet_NaN();
    auto max_boundary_surface_error_percent_opt = app.add_option(
        "--max-boundary-surface-error-percent",
        "Target maximum STEP boundary surface error when as percent of sqrt of surface area")
        ->check(CLI::PositiveNumber)
        ->check(CLI::Range(0.0, 100.0));

    double max_edge_length_percent = std::numeric_limits<double>::quiet_NaN();;
    auto max_edge_length_percent_opt = app.add_option("--max-edge-length-percent",
        max_edge_length_percent, "Target maximum edge length as percent of sqrt of surface area")
        ->check(CLI::PositiveNumber)
        ->check(CLI::Range(0.0, 100.0));


    double min_edge_length_percent = std::numeric_limits<double>::quiet_NaN();;
    auto min_edge_length_percent_opt = app.add_option("--min-edge-length-percent",
        min_edge_length_percent, "Target minimum edge length as percent of sqrt of surface area")
        ->check(CLI::PositiveNumber)
        ->check(CLI::Range(0.0, 100.0));

    double max_surface_error = std::numeric_limits<double>::quiet_NaN();
    app.add_option("-s,--max-surface-error", max_surface_error,
                   "Target maximum surface error when remeshing")
        ->excludes(max_surface_error_percent_opt)
        ->check(CLI::PositiveNumber);

    double max_boundary_surface_error = std::numeric_limits<double>::quiet_NaN();
    app.add_option("-b,--max-boundary-surface-error", max_boundary_surface_error,
                   "Target maximum STEP boundary surface error")
        ->excludes(max_boundary_surface_error_percent_opt)
        ->check(CLI::PositiveNumber);

    double max_edge_length = std::numeric_limits<double>::quiet_NaN();
    app.add_option("--max-edge-length", max_edge_length,"Target maximum edge length")
        ->excludes(max_edge_length_percent_opt);

    double min_edge_length = std::numeric_limits<double>::quiet_NaN();
    app.add_option("--min-edge-length", min_edge_length, "Target minimum edge length")
        ->excludes(min_edge_length_percent_opt);

    std::vector<std::pair<std::string, spdlog::level::level_enum>> level_map{
        {"trace", spdlog::level::trace},
        {"debug", spdlog::level::debug},
        {"info", spdlog::level::info},
        {"warn", spdlog::level::warn},
        {"error", spdlog::level::err},
        {"critical", spdlog::level::critical},
        {"off", spdlog::level::off}};
    auto log_level = spdlog::level::info;
    app.add_option("-l,--level", log_level, "Log level")->transform(CLI::CheckedTransformer(level_map));

    CLI11_PARSE(app, argc, argv);

    spdlog::set_level(log_level);
    spdlog::set_pattern("[%^%l%$] %v");

    std::filesystem::path input_path(input);
    std::string extension = boost::algorithm::to_lower_copy(input_path.extension().string());

    bool is_step = extension == ".stp" || extension == ".step";

    spdlog::info("parameters:");
    spdlog::info("  input                    = {}", input);

    if (is_step) {
        if (!shape_name.empty()) {
            spdlog::info("  shape                    = {}", shape_name);
        }
        else {
            spdlog::info("  shape                    = <largest> (default)");
        }
        if (shape_instance_index > 1) {
            spdlog::info("  shape instance           = {}", shape_instance_index);
        }
        else {
            spdlog::info("  shape instance           = {} (default)", shape_instance_index);
        }
    }

    if (iterations > 0) {
        spdlog::info("  iterations               = {}", iterations);
    }
    else {
        iterations = 10;
        spdlog::info("  iterations               = {} (default)", iterations);
    }

    if (std::isfinite(crease_angle)) {
        spdlog::info("  crease angle             = {} degrees", crease_angle);
    }
    else {
        crease_angle = 30.0;
        spdlog::info("  crease angle             = {} degrees (default)", crease_angle);
    }

    if (is_step) {
        if (unfreeze_step_boundaries) {
            spdlog::info("  unfreeze step boundaries = {}", unfreeze_step_boundaries);
        }
        else {
            spdlog::info("  unfreeze step boundaries = {} (default)", unfreeze_step_boundaries);
        }
    }

    std::string unit = "mesh_units";

    bool max_surface_error_from_percent = !std::isfinite(max_surface_error);
    bool max_edge_length_from_percent = !std::isfinite(max_edge_length);
    bool min_edge_length_from_percent = !std::isfinite(min_edge_length);
    bool max_boundary_surface_error_from_percent = !std::isfinite(max_boundary_surface_error);

    if (max_surface_error_from_percent) {
        if (std::isfinite(max_surface_error_percent)) {
            spdlog::info("  max surface error        = {:.2f} %", max_surface_error_percent);
        }
        else {
            max_surface_error_percent = 0.05;
            spdlog::info("  max surface error        = {:.2f} % (default)", max_surface_error_percent);
        }
    }
    else {
        spdlog::info("  max surface error        = {:.2f} {}", max_surface_error, unit);
    }

    if (is_step) {
        if (max_boundary_surface_error_from_percent) {
            if (std::isfinite(max_boundary_surface_error_percent)) {
                spdlog::info("  max boundary surface error = {:.2f} %", max_boundary_surface_error_percent);
            }
            else if (max_surface_error_from_percent) {
                max_boundary_surface_error_percent = max_surface_error_percent;
                spdlog::info("  max boundary surface error = {:.2f} % (max_surface_error)", max_surface_error_percent);
            }
            else {
                max_boundary_surface_error_from_percent = false;
                max_boundary_surface_error = max_surface_error;
                spdlog::info("  max boundary surface error = {:.2f} {} (max_surface_error)", max_boundary_surface_error, unit);
            }
        }
        else {
            spdlog::info("  max boundary surface error = {:.2f} {}", max_boundary_surface_error, unit);
        }
    }

    if (max_edge_length_from_percent) {
        if (std::isfinite(max_edge_length_percent)) {
            spdlog::info("  max edge length          = {:.2f} %", max_edge_length_percent);
        }
        else {
            max_edge_length_percent = 1.0;
            spdlog::info("  max edge length          = {:.2f} % (default)", max_edge_length_percent);
        }
    }
    else {
        spdlog::info("  max edge length          = {:.2f} {}", max_edge_length, unit);
    }
    if (min_edge_length_from_percent) {
        if (std::isfinite(min_edge_length_percent)) {
            spdlog::info("  min edge length          = {:.2f} %", min_edge_length_percent);
        }
        else {
            min_edge_length_percent = 0.1;
            spdlog::info("  min edge length          = {:.2f} % (default)", min_edge_length_percent);
        }
    }
    else {
        spdlog::info("  min edge length          = {:.2f} {}", min_edge_length, unit);
    }

    //tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, 1);

    std::mutex mutex;

    double surface_area = 0.0;
    std::vector<Mesh> meshes;
    std::shared_ptr<Component> selected_component;
    if (is_step) {

        /*
        for (auto& printer: Message::DefaultMessenger()->Printers()) {
            printer->SetTraceLevel(Message_Alarm);
        }
        */

        Handle(TDocStd_Document) doc;
        XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);

        spdlog::info("reading STEP file...");
        STEPCAFControl_Reader reader;
        auto status = reader.ReadFile(input.c_str());
        if (status != IFSelect_ReturnStatus::IFSelect_RetDone) {
            spdlog::critical("Failed to open {}.  Invalid STEP file.", input);
            return 1;
        }

        reader.Transfer(doc);

        TColStd_SequenceOfAsciiString unit_length_names, unit_angle_names, unit_solid_angle_names;
        reader.ChangeReader().FileUnits(unit_length_names, unit_angle_names, unit_solid_angle_names);
        if (!unit_length_names.IsEmpty()) {
            unit = normalizeUnit(unit_length_names.First().ToCString());
            spdlog::info("  length unit: {}", unit);
        }

        auto assembly = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
        TDF_LabelSequence labels;
        assembly->GetFreeShapes(labels);

        std::vector<std::shared_ptr<Component>> components;
        std::map<std::string, size_t> counts;
        for (auto i = 1; i <= labels.Length(); i++) {
            TDF_Label label = labels.Value(i);

            std::string name = getName(label);
            std::string reference_name;

            size_t index = ++counts[name];

            TDF_Label reference;
            if (assembly->GetReferredShape(label, reference)) {
                reference_name = getName(reference);
                label = reference;
            }

            components.push_back(
                std::make_shared<Component>(label, name, index, reference_name, 0));

            exploreAssembly(assembly, *components.back(), counts);
        }

        spdlog::info("  components:");
        std::deque<std::shared_ptr<Component>> queue;
        queue.insert(queue.end(), components.begin(), components.end());
        std::shared_ptr<Component> largest_component;
        double largest_surface_area = 0.0;
        while (!queue.empty()) {
            auto& component = queue.back();
            queue.pop_back();
            queue.insert(queue.end(), component->children.rbegin(), component->children.rend());

            bool matched = component->name == shape_name && component->index == shape_instance_index;
            if (matched || shape_name.empty()) {
                if (XCAFDoc_ShapeTool::IsSimpleShape(component->label)) {
                    component->shape = XCAFDoc_ShapeTool::GetShape(component->label);
                    if (!component->shape.IsNull()) {
                        GProp_GProps surface_props;
                        BRepGProp::SurfaceProperties(component->shape, surface_props);
                        component->surface_area = surface_props.Mass();

                        if (component->surface_area > largest_surface_area) {
                            largest_surface_area = component->surface_area;
                            largest_component = component;
                        }
                    }
                }
            }

            std::string matched_str;
            if (matched) {
                selected_component = component;
                matched_str = " <-----";
            }

            std::string surface_area_str = "";
            if (component->surface_area > 0.0) {
                surface_area_str = fmt::format(": {:.2f} sq {}", component->surface_area, unit);
            }

            spdlog::info("{}{} [instance = {}]{}{}", std::string(2 * (component->depth + 1), ' '),
                         component->name, component->index, surface_area_str, matched_str);
        }

        if (list_step_components) {
            return 0;
        }

        if (shape_name.empty()) {
            selected_component = largest_component;
        }


    }
    else if (mesh_formats.count(extension) != 0) {
        spdlog::info("reading mesh file...");
        Mesh mesh;
        if (!PMP::IO::read_polygon_mesh(input, mesh) || !CGAL::is_triangle_mesh(mesh)) {
            spdlog::critical("Not a valid input file: {}", input);
            return 1;
        }

        surface_area = PMP::area(mesh);
        spdlog::info("  surface area: {:.2f}", surface_area);
        spdlog::info("  loaded {} triangles", mesh.number_of_faces());
        meshes.push_back(mesh);
    }

    if (is_step) {
        if (!selected_component) {
            spdlog::error("No STEP component selected.");
            return 1;
        }
        spdlog::info("  shape: {}[instance = {}]", selected_component->name, selected_component->index);
        surface_area = selected_component->surface_area;
    }

    spdlog::info("  surface area: {:.2f} sq {}", surface_area, unit);

    if (max_surface_error_from_percent) {
        max_surface_error = std::sqrt(surface_area) * max_surface_error_percent / 100.0;
        spdlog::info("  max surface error: {:.2f} {} (from %)", max_surface_error, unit);
    }
    else {
        spdlog::info("  max surface error: {:.2f} {}", max_surface_error, unit);
    }

    if (is_step) {
        if (max_boundary_surface_error_from_percent) {
            max_boundary_surface_error = std::sqrt(surface_area) * max_boundary_surface_error_percent / 100.0;
            spdlog::info("  max boundary surface error: {:.2f} {} (from %)", max_boundary_surface_error, unit);
        }
        else {
            spdlog::info("  max bounary surface error: {:.2f} {}", max_boundary_surface_error, unit);
        }
    }

    if (max_edge_length_from_percent) {
        max_edge_length = std::sqrt(surface_area) * max_edge_length_percent / 100.0;
        spdlog::info("  max edge length: {:.2f} {} (from %)", max_edge_length, unit);
    }
    else {
        spdlog::info("  max edge length: {:.2f} {}", max_edge_length, unit);
    }
    if (min_edge_length_from_percent) {
        min_edge_length = std::sqrt(surface_area) * min_edge_length_percent / 100.0;
        spdlog::info("  min edge length: {:.2f} {} (from %)", min_edge_length, unit);
    }
    else {
        spdlog::info("  min edge length: {:.2f} {}", min_edge_length, unit);
    }

    std::vector<TopoDS_Face> segments;
    if (is_step) {


        double max_surface_error_auto = find_surface_error_param<Mesh>(selected_component->shape,
                                                                       min_edge_length, 10);
        max_surface_error = max_surface_error_auto;


        spdlog::info("  subdividing faces ...");


        selected_component->shape = subdivide_step_shape(selected_component->shape, min_edge_length,
                                                         max_edge_length, max_surface_error);

        save_shape_as_step("subdivided.step", selected_component->shape);

        spdlog::info("  tessalating ...");



        auto tesselation = tessalate_shape<Mesh>(selected_component->shape, max_surface_error);
        size_t total_faces = 0;
        for (const auto& [mesh, face]: tesselation) {
            meshes.push_back(mesh);
            segments.push_back(face);
            total_faces += mesh.number_of_faces();
        }


        spdlog::info("  tesselated component into {} faces over {} segments.", total_faces, meshes.size());
    }

    size_t total_faces_init = 0;
    for (auto& mesh: meshes) {
        total_faces_init += mesh.number_of_faces();
    }

    spdlog::info("    faces: {}", total_faces_init);


    if (is_step && raw_step_mesh) {
        saveOutput(outputs, meshes);
        return 0;
    }



    spdlog::info("  splitting long border edges ...");
    // TODO(malban): ensure that border splits are consistent on shared edges
    size_t border_num_before = 0;
    size_t border_num_after = 0;

    struct Less_xyz_3 {
        bool operator()(const typename Mesh::Point& p, const typename Mesh::Point& q) const {
          return std::lexicographical_compare(p.cartesian_begin(), p.cartesian_end(),
                                              q.cartesian_begin(), q.cartesian_end());
        }
    };

    typedef std::map<typename Mesh::Point, typename Mesh::Point, Less_xyz_3> PointMap;
    PointMap projection_map{Less_xyz_3()};

    for (size_t i = 0; i < meshes.size(); i++) {
        auto& mesh = meshes[i];
        std::vector<EdgeDescriptor> border_edges;
        PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(
            HalfEdge2Edge(mesh, border_edges)));
        border_num_before += border_edges.size();
        PMP::split_long_edges(border_edges, max_edge_length, mesh);
        border_edges.clear();
        PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(
            HalfEdge2Edge(mesh, border_edges)));

        std::vector<typename Mesh::Point> vertices;
        for (auto v: mesh.vertices()) {
            vertices.push_back(mesh.point(v));
        }

/*
        if (!no_projection && meshes.size() == segments.size()) {
            // TODO: reproject border vertices back to STEP border
            //project_to_step_border<Mesh>(segments[i], mesh);
        }

        std::vector<typename Mesh::Point> projected;
        for (auto v: mesh.vertices()) {
            projected.push_back(mesh.point(v));
        }

        for (size_t j = 0; j < vertices.size(); j++) {
            auto it = projection_map.find(vertices[j]);
            if (it == projection_map.end()) {
                projection_map[vertices[j]] = projected[j];
            }
            else {
                double dx = it->second.x() - projected[j].x();
                double dy = it->second.y() - projected[j].y();
                double dz = it->second.z() - projected[j].z();
                double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                double d1 = get_distance_to_face(segments[i], projected[j].x(), projected[j].y(), projected[j].z());

                if (d1 > 0.001) {
                    spdlog::error("d1 = {}", d1);
                }


                if (dist > 0.001) {
                    spdlog::warn("large projection delta: {}", dist);

                    double d2 = get_distance_to_face(segments[i], it->second.x(), it->second.y(), it->second.z());

                    spdlog::warn("d1 = {}, d2 = {}", d1, d2);
                }
            }
        }
        */

        border_num_after += border_edges.size();
    }



    if (is_step && raw_step_mesh) {
        saveOutput(outputs, meshes);
        return 0;
    }


    size_t total_faces_2 = 0;
    for (auto& mesh: meshes) {
        total_faces_2 += mesh.number_of_faces();
    }
    spdlog::info("    border edges: {} -> {}", border_num_before, border_num_after);
    spdlog::info("    faces: {} -> {}", total_faces_init, total_faces_2);

    if (!is_step) {
        spdlog::info("  splitting long crease edges ...");
        size_t crease_num_before = 0;
        size_t crease_num_after = 0;
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
                size_t num_before = 0;
                size_t num_after = 0;
                for (size_t i=r.begin(); i!=r.end(); ++i) {
                    Mesh& mesh = meshes[i];
                    Mesh::Property_map<Mesh::Edge_index, bool> crease_features =
                        mesh.add_property_map<Mesh::Edge_index, bool>("crease", false).first;
                    PMP::detect_sharp_edges(mesh, crease_angle, crease_features);
                    std::vector<EdgeDescriptor> crease_edges;
                    for(const auto& edge: edges(mesh)) {
                        if (get(crease_features, edge)) {
                            crease_edges.push_back(edge);
                        }
                    }
                    num_before += crease_edges.size();
                    PMP::split_long_edges(crease_edges, max_edge_length, mesh);
                    mesh.remove_property_map(crease_features);
                    crease_edges.clear();

                    Mesh::Property_map<Mesh::Edge_index, bool> crease_features2 =
                        mesh.add_property_map<Mesh::Edge_index, bool>("crease", false).first;
                    PMP::detect_sharp_edges(mesh, crease_angle, crease_features2);
                    for(const auto& edge: edges(mesh)) {
                        if (get(crease_features2, edge)) {
                            crease_edges.push_back(edge);
                        }
                    }

                    num_after += crease_edges.size();
                }
                std::scoped_lock<std::mutex> lock(mutex);
                crease_num_before += num_before;
                crease_num_after += num_after;
            });
        size_t total_faces_3 = 0;
        for (auto& mesh: meshes) {
            total_faces_3 += mesh.number_of_faces();
        }
        spdlog::info("    crease edges: {} -> {}", crease_num_before, crease_num_after);
        spdlog::info("    faces: {} -> {}", total_faces_2, total_faces_3);
    }

    double max_remeshing_surface_error = std::min(max_surface_error, min_edge_length * 0.1);

    spdlog::info("  adaptive isotropic remeshing ...");
    for (int i = 0; i < iterations; i++) {
        spdlog::info("    iteration {}", i + 1);
        spdlog::info("      remeshing ...");
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
                for (size_t m=r.begin(); m!=r.end(); ++m) {
                    Mesh& mesh = meshes[m];

                    const std::pair edge_min_max{min_edge_length, max_edge_length};
                    PMP::Adaptive_sizing_field<Mesh> sizing_field(max_remeshing_surface_error,
                                                                  edge_min_max, faces(mesh), mesh);
                    auto crease_map = mesh.property_map<Mesh::Edge_index, bool>("crease").first;

                    try {
                        if (is_step) {
                            PMP::isotropic_remeshing(faces(mesh), sizing_field, mesh,
                                CGAL::parameters::number_of_iterations(1)
                                                 .number_of_relaxation_steps(3)
                                                 .protect_constraints(true));
                        }
                        else {
                            PMP::isotropic_remeshing(faces(mesh), sizing_field, mesh,
                                CGAL::parameters::number_of_iterations(1)
                                                 .number_of_relaxation_steps(3)
                                                 .edge_is_constrained_map(crease_map)
                                                 .protect_constraints(true));
                        }
                    }
                    catch (const std::out_of_range& ex) {
                        spdlog::warn("Out of range exception when remeshing segment {}.", m);
                    }
                }});

        if (is_step) {
            spdlog::info("      projecting ...");
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
                    for (size_t m=r.begin(); m!=r.end(); ++m) {
                        double weight = 1.0 / (iterations - i);
                        project_to_step<Mesh>(selected_component->shape, segments[m], meshes[m], weight);
            }});
        }
    }

    size_t total_faces_remeshed = 0;
    for (auto& mesh: meshes) {
        total_faces_remeshed += mesh.number_of_faces();
    }
    //spdlog::info("    faces: {} -> {}", total_faces_3, total_faces_remeshed);

    // auto merged = merge_meshes(meshes, Point_traits());
    // spdlog::info("  merged faces: {}",  merged.number_of_faces());
    saveOutput(outputs, meshes);


    return 0;
}
