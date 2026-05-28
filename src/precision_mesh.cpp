#include <filesystem>
#include <limits>
#include <map>
#include <string>
#include <thread>
#include <unordered_map>
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
#include <gp_Trsf.hxx>
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
#include <precision_mesh/remeshing.h>
#include <precision_mesh/step_face_util.h>
#include <precision_mesh/step_projection.h>
#include <precision_mesh/step_reader.h>
#include <precision_mesh/step_subdivision.h>
#include <precision_mesh/unit_conversion.h>
#include <precision_mesh/step_tessellation.h>
#include <precision_mesh/stl.h>

namespace PMP = CGAL::Polygon_mesh_processing;

std::unordered_set<std::string> mesh_formats = {".obj", ".off", ".ply", ".stl", ".ts", ".vtp"};

bool saveOutput(const std::vector<std::string>& outputs, const std::vector<Mesh>& meshes,
                const std::vector<TopoDS_Face>& faces,
                const std::unordered_map<size_t, int>& component_map, float scale=1)
{
    bool ok = true;
    for (const auto& output: outputs) {

        std::filesystem::path output_path(output);
        std::string extension = boost::algorithm::to_lower_copy(output_path.extension().string());

        if (extension == ".ply") {
            spdlog::info("saving mesh to: {}", output);
            ok = saveComponentsToPly<Point_traits>(output, meshes, faces, component_map, scale) && ok;
        }
        else if (extension == ".stl") {
            spdlog::info("saving mesh to: {}", output);
            ok = saveComponentsToStl<Point_traits>(output, meshes, scale) && ok;
        }
        else {
            spdlog::error("Unsupported output format: {}", extension);
            ok = false;
        }
    }
    return ok;
}

int main(int argc, char **argv) {

    CLI::App app{"A flexible STEP to mesh conversion tool and general adaptive isotropic remesher.", "precision_mesh"};
    argv = app.ensure_utf8(argv);

    std::string input;
    app.add_option("-i,--input", input, "Input file (.obj|.off|.p21|.ply|.step|.stl|.stp|.ts|.vtp)")
        ->check(CLI::ExistingFile)
        ->required();

    std::string shape_name;
    auto shape_name_opt = app.add_option("-n,--name", shape_name, "STEP file shape name.");

    int shape_index = -1;
    app.add_option("--index", shape_index, "STEP file shape index.")
        ->check(CLI::NonNegativeNumber)
        ->excludes(shape_name_opt);

    std::vector<std::string> outputs;
    app.add_option("-o,--output", outputs, "Output file (.obj|.off|.ply|.stl|.ts|.vtp)")->take_all();

    std::string output_unit;
    app.add_option("--output-units", output_unit, "Output units (mm|cm|m|in|ft|yd)");

    bool list_step_components = false;
    app.add_flag("--list",  list_step_components, "List STEP components.");

    bool raw_step_mesh = false;
    app.add_flag("-r,--raw-step-mesh",  raw_step_mesh, "Generate raw mesh from STEP component without remeshing.");

    int iterations = 0;
    app.add_option("--iterations", iterations, "Iterations")->check(CLI::PositiveNumber);

    bool no_subdivision = false;
    app.add_flag("--no-subdivision",  no_subdivision, "Skip STEP shape subdivision.");

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
        "--max-boundary-surface-error-percent", max_boundary_surface_error_percent,
        "Target maximum STEP boundary surface error as percent of sqrt of surface area")
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

    bool is_step = extension == ".stp" || extension == ".step" || extension == ".p21";

    spdlog::info("parameters:");
    spdlog::info("  input                    = {}", input);

    if (is_step) {
        if (!shape_name.empty()) {
            spdlog::info("  shape                    = {}", shape_name);
        }
        else if (shape_index >= 0) {
            spdlog::info("  shape index              = {}", shape_index);
        }
        else {
            spdlog::info("  shape                    = <largest> (default)");
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
    if (!output_unit.empty()) {
        unit = output_unit;
    }

    float conversion_scale = 1.0;

    bool max_surface_error_from_percent = !std::isfinite(max_surface_error);
    bool max_edge_length_from_percent = !std::isfinite(max_edge_length);
    bool min_edge_length_from_percent = !std::isfinite(min_edge_length);
    bool max_boundary_surface_error_from_percent = !std::isfinite(max_boundary_surface_error);

    bool max_surface_error_from_default = false;

    if (max_surface_error_from_percent) {
        if (std::isfinite(max_surface_error_percent)) {
            spdlog::info("  max surface error        = {:.4f} %", max_surface_error_percent);
        }
        else {
            max_surface_error_percent = 0.05;
            spdlog::info("  max surface error        = {:.4f} % (default)", max_surface_error_percent);
            max_surface_error_from_default = true;
        }
    }
    else {
        spdlog::info("  max surface error        = {:.4f} {}", max_surface_error, unit);
    }

    if (is_step) {
        if (max_boundary_surface_error_from_percent) {
            if (std::isfinite(max_boundary_surface_error_percent)) {
                spdlog::info("  max boundary surface error = {:.2f} %", max_boundary_surface_error_percent);
            }
            else if (max_surface_error_from_percent) {
                max_boundary_surface_error_percent = max_surface_error_percent;
                spdlog::info("  max boundary surface error = {:.4f} % (max_surface_error)", max_surface_error_percent);
            }
            else {
                max_boundary_surface_error_from_percent = false;
                max_boundary_surface_error = max_surface_error;
                spdlog::info("  max boundary surface error = {:.4f} {} (max_surface_error)", max_boundary_surface_error, unit);
            }
        }
        else {
            spdlog::info("  max boundary surface error = {:.4f} {}", max_boundary_surface_error, unit);
        }
    }

    if (max_edge_length_from_percent) {
        if (std::isfinite(max_edge_length_percent)) {
            spdlog::info("  max edge length          = {:.4f} %", max_edge_length_percent);
        }
        else {
            max_edge_length_percent = 1.0;
            spdlog::info("  max edge length          = {:.4f} % (default)", max_edge_length_percent);
        }
    }
    else {
        spdlog::info("  max edge length          = {:.4f} {}", max_edge_length, unit);
    }
    if (min_edge_length_from_percent) {
        if (std::isfinite(min_edge_length_percent)) {
            spdlog::info("  min edge length          = {:.4f} %", min_edge_length_percent);
        }
        else {
            min_edge_length_percent = 0.1;
            spdlog::info("  min edge length          = {:.4f} % (default)", min_edge_length_percent);
        }
    }
    else {
        spdlog::info("  min edge length          = {:.4f} {}", min_edge_length, unit);
    }

    if (!is_step && !output_unit.empty()) {
        spdlog::error("Output units only valid for a STEP input.");
        return 1;
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

        if (output_unit.empty()) {
            output_unit = unit;
            spdlog::info("  output unit: {} (STEP)", output_unit);
        }
        else {
            std::string normalized_output_unit = normalizeUnit(output_unit);
            if (!isKnownUnit(normalized_output_unit)) {
                spdlog::error("Output unit {} is not supported.", output_unit);
                return 1;
            }
            output_unit = normalized_output_unit;
            spdlog::info("  output unit: {} ", output_unit);
        }

        if (unit != output_unit) {
            if (!isKnownUnit(unit)) {
                spdlog::error("Unable to convert between input unit: {} and output unit: {}.",
                              unit, output_unit);
                return 1;
            }

            conversion_scale = getUnitConversionScale(unit, output_unit);
        }
        spdlog::info("  unit conversion scale: {} ", conversion_scale);

        auto assembly = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
        TDF_LabelSequence labels;
        assembly->GetFreeShapes(labels);

        std::vector<std::shared_ptr<Component>> components;
        collectComponents(assembly, labels, "", components);

        spdlog::info("  components:");
        std::shared_ptr<Component> largest_component;
        double largest_surface_area = 0.0;
        for (size_t i = 0; i < components.size(); i++) {
            const auto& component = components[i];

            GProp_GProps surface_props;
            BRepGProp::SurfaceProperties(component->shape, surface_props);
            component->surface_area = surface_props.Mass();

            if (component->surface_area > largest_surface_area) {
                largest_surface_area = component->surface_area;
                largest_component = component;
            }

            std::string matched_str;
            bool matched = component->qualified_name == shape_name ||
                           shape_index == static_cast<int>(i);
            if (matched) {
                selected_component = component;
                matched_str = " <-----";
            }

            std::string surface_area_str = "";
            if (component->surface_area > 0.0) {
                surface_area_str = fmt::format(": {:.2f} sq {}", component->surface_area * conversion_scale * conversion_scale, output_unit);
            }

            spdlog::info("    [{}] {}{}{}", i, component->qualified_name, surface_area_str,
                         matched_str);
        }

        if (list_step_components) {
            return 0;
        }

        if (shape_name.empty() && shape_index < 0) {
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
        surface_area = selected_component->surface_area * conversion_scale * conversion_scale;
    }

    spdlog::info("  surface area: {:.4f} sq {}", surface_area, output_unit);

    if (max_surface_error_from_percent) {
        max_surface_error = std::sqrt(surface_area) * max_surface_error_percent / 100.0;
        spdlog::info("  max surface error: {:.4f} {} (from %)", max_surface_error, output_unit);
    }
    else {
        spdlog::info("  max surface error: {:.4f} {}", max_surface_error, output_unit);
    }

    if (is_step) {
        if (max_boundary_surface_error_from_percent) {
            max_boundary_surface_error = std::sqrt(surface_area) * max_boundary_surface_error_percent / 100.0;
            spdlog::info("  max boundary surface error: {:.4f} {} (from %)", max_boundary_surface_error, output_unit);
        }
        else {
            spdlog::info("  max boundary surface error: {:.4f} {}", max_boundary_surface_error, output_unit);
        }
    }

    if (max_edge_length_from_percent) {
        max_edge_length = std::sqrt(surface_area) * max_edge_length_percent / 100.0;
        spdlog::info("  max edge length: {:.4f} {} (from %)", max_edge_length, output_unit);
    }
    else {
        spdlog::info("  max edge length: {:.4f} {}", max_edge_length, output_unit);
    }
    if (min_edge_length_from_percent) {
        min_edge_length = std::sqrt(surface_area) * min_edge_length_percent / 100.0;
        spdlog::info("  min edge length: {:.4f} {} (from %)", min_edge_length, output_unit);
    }
    else {
        spdlog::info("  min edge length: {:.4f} {}", min_edge_length, output_unit);
    }

    std::vector<TopoDS_Face> segments;
    std::vector<TopoDS_Face> original_faces;
    std::unordered_map<size_t, int> component_map;

    min_edge_length /= conversion_scale;
    max_edge_length /= conversion_scale;
    max_surface_error /= conversion_scale;
    max_boundary_surface_error /= conversion_scale;

    if (is_step) {

        if (max_surface_error_from_default) {
            double max_surface_error_auto = find_surface_error_param(selected_component->shape,
                                                                           min_edge_length, 10, 
                                                                           conversion_scale);
            max_surface_error = max_surface_error_auto;
        }

        original_faces = get_faces(selected_component->shape);

        FaceMap face_map;
        if (!no_subdivision && !raw_step_mesh) {
            spdlog::info("  subdividing faces ...");


            std::tie(selected_component->shape, face_map) = 
                subdivide_step_shape(selected_component->shape, min_edge_length, max_edge_length, 
                                     max_surface_error);
        }

        spdlog::info("  tessellating ...");

        auto tessellation = tessellate_shape(selected_component->shape, max_surface_error);
        size_t total_faces = 0;
        for (const auto& [mesh, face]: tessellation) {
            meshes.push_back(mesh);
            segments.push_back(face);
            total_faces += mesh.number_of_faces();
        }

        spdlog::info("  tessellated component into {} faces over {} segments.", total_faces, meshes.size());

        // Create mapping of subdivided components to the original components prior to subdivision
        if (!face_map.empty()) {
            for (size_t i = 0; i < segments.size(); i++) {
                component_map[i] = face_map[segments[i]];
            }
        }
    }



    size_t total_faces_init = 0;
    for (auto& mesh: meshes) {
        total_faces_init += mesh.number_of_faces();
    }

    spdlog::info("    faces: {}", total_faces_init);

    if (is_step && raw_step_mesh) {
        return saveOutput(outputs, meshes, segments, {}, conversion_scale) ? 0 : 1;
    }

    split_border_edges(meshes, max_edge_length);

    if (!is_step) {
        split_crease_edges(meshes, crease_angle, max_edge_length);
    }

    spdlog::info("  adaptive isotropic remeshing ...");

    WireProjectorCachePtr wire_projectors;
    std::vector<std::unique_ptr<StepProjector>> surface_projectors;
    std::vector<std::unique_ptr<StepBorderProjector>> border_projectors;
    if (is_step) {
        spdlog::info("    creating edge projectors ...");
        wire_projectors = get_edge_vertex_wire_projectors(selected_component->shape);

        if (!no_projection) {
            spdlog::info("    initializing surface projectors ...");
            for (const auto& segment : segments) {
                surface_projectors.push_back(std::make_unique<StepProjector>(segment));
                border_projectors.push_back(std::make_unique<StepBorderProjector>(segment));
            }
        }
    }

    RemeshParams rparams{min_edge_length, max_edge_length, max_surface_error,
                         iterations, is_step, no_projection};
    remesh_and_project(meshes, segments, wire_projectors, surface_projectors, border_projectors,
                       rparams);

    return saveOutput(outputs, meshes, original_faces, component_map, conversion_scale) ? 0 : 1;
}
