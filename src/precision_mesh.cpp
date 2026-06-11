#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
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
#include <Geom_Curve.hxx>
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
#include <Precision.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
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

#ifdef PRECISION_MESH_HAS_VIEWER
#include <precision_mesh/viewer.h>
#endif

namespace PMP = CGAL::Polygon_mesh_processing;

std::unordered_set<std::string> mesh_formats = {".obj", ".off", ".ply", ".stl", ".ts", ".vtp"};

#ifdef PRECISION_MESH_HAS_VIEWER
static void push_to_viewer(Viewer* viewer,
                            const std::vector<Mesh>& meshes,
                            const std::vector<TopoDS_Face>& segments,
                            const std::vector<TopoDS_Face>& original_faces,
                            const std::unordered_map<size_t, int>& component_map,
                            bool is_step,
                            const std::string& label,
                            double wire_deflection,
                            const std::vector<float>& free_edge_lines = {},
                            const std::vector<float>& healed_edge_lines = {})
{
    if (!viewer) return;

    // Count how many tessellation segments belong to each original face
    std::unordered_map<int, int> subdiv_counts;
    for (size_t i = 0; i < meshes.size(); i++) {
        int orig = component_map.count(i) ? component_map.at(i) : (int)i;
        subdiv_counts[orig]++;
    }

    // Pre-sample original face wires once per unique original face index
    std::unordered_map<int, std::vector<float>> orig_wire_cache;
    if (is_step) {
        for (size_t i = 0; i < meshes.size(); i++) {
            int orig_idx = component_map.count(i) ? component_map.at(i) : (int)i;
            if (!orig_wire_cache.count(orig_idx) && orig_idx < (int)original_faces.size())
                orig_wire_cache[orig_idx] = sample_face_wire(original_faces[orig_idx], wire_deflection);
        }
    }

    std::vector<ViewerSegment> vsegs;
    vsegs.reserve(meshes.size());
    for (size_t i = 0; i < meshes.size(); i++) {
        int orig_idx = component_map.count(i) ? component_map.at(i) : (int)i;
        int ft       = 0;
        float area   = 0.0f;
        std::string desc = "Mesh";
        std::vector<float> subseg_wire, orig_face_wire;
        if (is_step && i < segments.size()) {
            ft   = get_face_type(segments[i]);
            area = get_face_area(segments[i]);
            desc = get_face_description(segments[i]);
            subseg_wire = sample_face_wire(segments[i], wire_deflection);
            if (orig_wire_cache.count(orig_idx))
                orig_face_wire = orig_wire_cache[orig_idx];
        }
        int sdivs = subdiv_counts.count(orig_idx) ? subdiv_counts[orig_idx] : 1;
        vsegs.push_back({&meshes[i], ft, (int)i, orig_idx, sdivs, area, desc,
                         std::move(subseg_wire), std::move(orig_face_wire)});
    }
    viewer->update(vsegs, label, free_edge_lines, healed_edge_lines);
}
#endif

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

    bool validate = false;
    app.add_flag("--validate", validate,
                 "Validate the final tessellation against the BREP (vertex placement, "
                 "watertightness).  Slow; for diagnostics.");

    bool validate_surface_error = false;
    app.add_flag("--validate-surface-error", validate_surface_error,
                 "During --validate, also sample triangle interiors (centroid + edge "
                 "midpoints) against the BREP to measure chord deviation.  Expensive.");

    std::string validate_report;
    app.add_option("--validate-report", validate_report,
                   "Write validation metrics, health counters, and all run parameters to "
                   "a YAML report file at the given path (implies --validate).");

    bool no_tess_repair = false;
    app.add_flag("--no-tess-repair", no_tess_repair,
                 "Disable the inverted/holed/dropped-corner tessellation repair (diagnostic).");

    bool no_uv_tess = false;
    app.add_flag("--no-uv-tess", no_uv_tess,
                 "Disable UV-grid CDT tessellation for regularly-curved faces; fall back to "
                 "BRepMesh + isotropic remeshing (diagnostic).");

    std::string dump_cdt_dir;
    app.add_option("--dump-cdt-obj", dump_cdt_dir,
                   "Write each UV-grid CDT face's 2D triangulation to Face_<idx>_CDT.obj in "
                   "the given directory, created if missing (diagnostic).");

#ifdef PRECISION_MESH_HAS_VIEWER
    bool enable_display = false;
    app.add_flag("--display", enable_display, "Open interactive 3D viewer during processing.");
#endif

    int iterations = -1;
    app.add_option("--iterations", iterations, "Iterations (0 = skip remeshing)")->check(CLI::NonNegativeNumber);

    bool no_subdivision = false;
    app.add_flag("--no-subdivision",  no_subdivision, "Skip STEP shape subdivision.");

    bool boundary_mesh_mode = false;
    app.add_flag("--boundary-mesh", boundary_mesh_mode,
                 "Visualise raw sub-face boundaries (fan-polygon meshes, no BRepMesh).");

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

    if (!validate_report.empty()) {
        validate = true;
    }

    if (!dump_cdt_dir.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(dump_cdt_dir, ec);
        if (ec) {
            spdlog::error("Failed to create --dump-cdt-obj directory '{}': {}",
                          dump_cdt_dir, ec.message());
            return 1;
        }
    }

    spdlog::set_level(log_level);
    spdlog::set_pattern("[%^%l%$] %v");

    const auto run_start = std::chrono::steady_clock::now();

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

    if (iterations < 0) {
        iterations = 10;
        spdlog::info("  iterations               = {} (default)", iterations);
    }
    else {
        spdlog::info("  iterations               = {}", iterations);
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

#ifdef PRECISION_MESH_HAS_VIEWER
    Viewer* viewer_ptr = nullptr;
#endif

    auto do_pipeline = [&]() -> int {

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
    [[maybe_unused]] std::vector<float> free_edge_lines;    // input BREP free boundary edges, for the viewer (yellow)
    [[maybe_unused]] std::vector<float> healed_edge_lines; // pre-repair open edges, for the viewer (green)

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

        // DEBUG: report non-conformal topology in the INPUT shape -- edges adjacent to !=2
        // faces are free (open shell / un-stitched import) and tessellate to cracks that no
        // amount of meshing can stitch.  Periodic SEAM edges (closed on their face) belong to
        // one face but are NOT boundaries -- exclude them so the count reflects real free edges.
        {
            TopTools_IndexedDataMapOfShapeListOfShape edge_faces;
            TopExp::MapShapesAndAncestors(selected_component->shape, TopAbs_EDGE, TopAbs_FACE, edge_faces);
            int free_edges = 0, free_degenerate = 0, seam_edges = 0, nonmanifold = 0, shown = 0;
            for (int i = 1; i <= edge_faces.Extent(); i++) {
                const TopTools_ListOfShape& faces = edge_faces.FindFromIndex(i);
                int nf = faces.Extent();
                if (nf == 2) continue;
                const TopoDS_Edge& e = TopoDS::Edge(edge_faces.FindKey(i));
                if (nf >= 1 && BRep_Tool::IsClosed(e, TopoDS::Face(faces.First()))) { seam_edges++; continue; }
                if (nf >= 2) { nonmanifold++; continue; }
                Standard_Real f, l; Handle(Geom_Curve) c = BRep_Tool::Curve(e, f, l);
                if (c.IsNull() || l - f < 1e-15) { free_degenerate++; continue; }
                free_edges++;
                if (shown < 12) {
                    gp_Pnt a = c->Value(f), b = c->Value(l);
                    spdlog::debug("    [freeedge] nfaces={} ({:.4f},{:.4f},{:.4f})->({:.4f},{:.4f},{:.4f})",
                                  nf, a.X(),a.Y(),a.Z(), b.X(),b.Y(),b.Z());
                    shown++;
                }
            }
            spdlog::info("    OCCT precision: Confusion={:.3e} PConfusion={:.3e}",
                         Precision::Confusion(), Precision::PConfusion());
            if (free_degenerate)
                spdlog::info("    input topology: {} edges, {} free (real boundary), {} degenerate free (no 3D curve), {} seam, {} non-manifold",
                             edge_faces.Extent(), free_edges, free_degenerate, seam_edges, nonmanifold);
            else
                spdlog::info("    input topology: {} edges, {} free (real boundary), {} seam, {} non-manifold",
                             edge_faces.Extent(), free_edges, seam_edges, nonmanifold);

            // A degenerate-free edge (null curve, 1 face) is only valid as a surface pole,
            // where it loops back to itself (FirstVertex == LastVertex, same object).  If the
            // two endpoint vertices are distinct objects the edge is not a pole -- it connects
            // two different BREP vertices along a single face with no matching opposite face.
            // That face's tessellation will have an open boundary edge there with nothing to
            // stitch against: the BREP is the cause of the tessellation gap.
            for (int i = 1; i <= edge_faces.Extent(); i++) {
                const TopTools_ListOfShape& faces = edge_faces.FindFromIndex(i);
                if (faces.Extent() != 1) continue;
                const TopoDS_Edge& e = TopoDS::Edge(edge_faces.FindKey(i));
                Standard_Real f, l;
                if (!BRep_Tool::Curve(e, f, l).IsNull()) continue;  // not degenerate-free
                TopoDS_Vertex v1 = TopExp::FirstVertex(e);
                TopoDS_Vertex v2 = TopExp::LastVertex(e);
                if (v1.IsNull() || v2.IsNull() || v1.IsSame(v2)) continue;  // proper pole
                gp_Pnt p1 = BRep_Tool::Pnt(v1), p2 = BRep_Tool::Pnt(v2);
                spdlog::warn("    [degenerate-free edge, distinct vertices] "
                             "({:.6f},{:.6f},{:.6f}) vs ({:.6f},{:.6f},{:.6f}) dist={:.3e} -- "
                             "not a pole; will produce open tessellation boundary",
                             p1.X(),p1.Y(),p1.Z(), p2.X(),p2.Y(),p2.Z(), p1.Distance(p2));
            }
        }

        // A degenerate BREP edge (BRep_Tool::Degenerated) collapses to a point -- its two
        // endpoint vertices must be the same topological object (v1.IsSame(v2)).  If they are
        // not, the BREP has a provable inconsistency with no geometric tolerance involved: the
        // topology declares the edge a point but assigns two distinct vertex identities to its
        // endpoints.  In a topologically closed BREP this is the only way a watertight topology
        // can produce a non-watertight tessellation, because BRepMesh evaluates each face's
        // pcurves independently and the two vertex positions will not agree exactly.
        {
            TopTools_IndexedMapOfShape edge_map;
            TopExp::MapShapes(selected_component->shape, TopAbs_EDGE, edge_map);
            int inconsistent = 0;
            for (int i = 1; i <= edge_map.Extent(); i++) {
                const TopoDS_Edge& e = TopoDS::Edge(edge_map.FindKey(i));
                if (!BRep_Tool::Degenerated(e)) continue;
                TopoDS_Vertex v1 = TopExp::FirstVertex(e);
                TopoDS_Vertex v2 = TopExp::LastVertex(e);
                if (v1.IsNull() || v2.IsNull() || v1.IsSame(v2)) continue;
                gp_Pnt p1 = BRep_Tool::Pnt(v1), p2 = BRep_Tool::Pnt(v2);
                spdlog::warn("    [degenerate edge, distinct vertices] ({:.6f},{:.6f},{:.6f}) vs "
                             "({:.6f},{:.6f},{:.6f}) dist={:.3e} -- will cause tessellation gap",
                             p1.X(),p1.Y(),p1.Z(), p2.X(),p2.Y(),p2.Z(),
                             p1.Distance(p2));
                inconsistent++;
            }
            if (inconsistent > 0)
                spdlog::warn("    {} degenerate edge(s) with distinct endpoint vertices",
                             inconsistent);
        }

        // For each edge with a non-null 3D curve, check that the chord from curve(f) to
        // curve(l) bridges the vertex-to-vertex gap within the combined vertex tolerance.
        // Shortfalls beyond tolerance indicate the stored parameter range doesn't correctly
        // connect the edge's bounding vertices.  Diagnostic only -- no shape modification.
        {
            TopTools_IndexedMapOfShape emap;
            TopExp::MapShapes(selected_component->shape, TopAbs_EDGE, emap);
            int inconsistent = 0;
            for (int i = 1; i <= emap.Extent(); i++) {
                const TopoDS_Edge& e = TopoDS::Edge(emap.FindKey(i));
                Standard_Real f, l;
                Handle(Geom_Curve) c = BRep_Tool::Curve(e, f, l);
                if (c.IsNull()) continue;
                TopoDS_Vertex v1 = TopExp::FirstVertex(e);
                TopoDS_Vertex v2 = TopExp::LastVertex(e);
                if (v1.IsNull() || v2.IsNull()) continue;
                gp_Pnt pf = c->Value(f), pl = c->Value(l);
                gp_Pnt q1 = BRep_Tool::Pnt(v1), q2 = BRep_Tool::Pnt(v2);
                double chord    = pf.Distance(pl);
                double vgap     = q1.Distance(q2);
                double tol1     = BRep_Tool::Tolerance(v1), tol2 = BRep_Tool::Tolerance(v2);
                double shortfall = vgap - chord;
                if (shortfall <= tol1 + tol2) continue;
                inconsistent++;
                spdlog::warn("    [edge-vertex mismatch] chord={:.3e} vertex-gap={:.3e} "
                             "shortfall={:.3e} (tol sum={:.3e}) -- "
                             "V1=({:.5f},{:.5f},{:.5f}) V2=({:.5f},{:.5f},{:.5f})",
                             chord, vgap, shortfall, tol1 + tol2,
                             q1.X(),q1.Y(),q1.Z(), q2.X(),q2.Y(),q2.Z());
            }
            if (inconsistent > 0) {
                spdlog::warn("    {} edge(s) with chord-vertex mismatch beyond tolerance",
                             inconsistent);
            }
        }

#ifdef PRECISION_MESH_HAS_VIEWER
        // Sample the input's free boundary edges once (pre-subdivision) for the viewer overlay.
        free_edge_lines = sample_free_edges(selected_component->shape, max_surface_error);
        spdlog::info("    free-edge overlay: {} GL_LINES segments (0 = no renderable free edges)",
                     free_edge_lines.size() / 6);
#endif

        FaceMap face_map;
        if (!no_subdivision && !raw_step_mesh) {
            spdlog::info("  subdividing faces ...");


            std::tie(selected_component->shape, face_map) = 
                subdivide_step_shape(selected_component->shape, min_edge_length, max_edge_length, 
                                     max_surface_error);
        }

        spdlog::info("  tessellating ...");

        auto tessellation = boundary_mesh_mode
            ? boundary_meshes(selected_component->shape)
            : tessellate_shape(selected_component->shape, max_surface_error, !no_tess_repair);

        if (is_step && !boundary_mesh_mode) {
#ifdef PRECISION_MESH_HAS_VIEWER
            // Capture open boundary edges before repair so they can be shown in green.
            if (viewer_ptr) {
                std::map<std::tuple<double,double,double>, int> pid;
                std::vector<std::array<double,3>> pos;
                std::map<std::pair<int,int>, int> ecnt;
                auto id_of = [&](const Mesh::Point& p) {
                    double x=CGAL::to_double(p.x()), y=CGAL::to_double(p.y()), z=CGAL::to_double(p.z());
                    auto key = std::make_tuple(x,y,z);
                    auto it = pid.find(key);
                    if (it != pid.end()) return it->second;
                    int id = (int)pos.size(); pid[key]=id; pos.push_back({x,y,z}); return id;
                };
                for (const auto& [mesh, face] : tessellation) {
                    for (auto f : mesh.faces()) {
                        int ids[3]; int k=0;
                        for (auto v : mesh.vertices_around_face(mesh.halfedge(f)))
                            if (k<3) ids[k++] = id_of(mesh.point(v));
                        if (k!=3) continue;
                        for (int e=0;e<3;e++) {
                            int a=ids[e],b=ids[(e+1)%3]; if(a>b)std::swap(a,b);
                            ecnt[{a,b}]++;
                        }
                    }
                }
                for (const auto& [e, cnt] : ecnt) {
                    if (cnt != 1) continue;
                    const auto& A = pos[e.first]; const auto& B = pos[e.second];
                    healed_edge_lines.push_back((float)A[0]); healed_edge_lines.push_back((float)A[1]); healed_edge_lines.push_back((float)A[2]);
                    healed_edge_lines.push_back((float)B[0]); healed_edge_lines.push_back((float)B[1]); healed_edge_lines.push_back((float)B[2]);
                }
            }
#endif
            repair_open_boundary_loops(tessellation, min_edge_length);
        }

        size_t total_faces = 0;
        for (const auto& [mesh, face]: tessellation) {
            meshes.push_back(mesh);
            segments.push_back(face);
            total_faces += mesh.number_of_faces();
        }

        spdlog::info("  tessellated component into {} faces over {} segments.", total_faces, meshes.size());
#ifdef PRECISION_MESH_HAS_VIEWER
        push_to_viewer(viewer_ptr, meshes, segments, original_faces, component_map, is_step,
                       "After Tessellation", max_surface_error, free_edge_lines, healed_edge_lines);
#endif

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

    // CDT engagement counters, reported by --validate-report.
    size_t cdt_faces_attempted = 0;
    size_t cdt_faces_succeeded = 0;

    if (iterations > 0) {
        split_border_edges(meshes, max_edge_length);

        if (!is_step) {
            split_crease_edges(meshes, crease_angle, max_edge_length);
        }

        if (is_step) {
            snap_border_midpoints_to_brep(meshes, segments, max_edge_length);
        }

#ifdef PRECISION_MESH_HAS_VIEWER
        push_to_viewer(viewer_ptr, meshes, segments, original_faces, component_map, is_step,
                       "After Border Splitting", max_surface_error, free_edge_lines, healed_edge_lines);
#endif

        // UV-grid CDT tessellation for regularly-curved STEP faces (cylinder, and later
        // cone/torus/revolution).  Replaces the BRepMesh interior with a UV-grid CDT
        // seeded by Jacobian-corrected Steiner points, giving uniform normal distribution.
        // Runs once (not iterative — the CDT interior is already at the target density).
        // These faces bypass the standard isotropic remeshing step.
        struct UvFace { size_t idx; int u_steps, v_steps; };
        std::vector<UvFace> uv_faces;
        std::vector<bool> use_uv_tess(meshes.size(), false);
        if (is_step && !no_uv_tess) {
            for (size_t i = 0; i < segments.size(); i++) {
                auto steps = get_face_tessellation_steps(segments[i], min_edge_length,
                                                         max_edge_length, max_surface_error);
                // u=v=1 means the BRepMesh interior is already at target density.
                if (!cdt_eligible(steps) || (steps.u_steps <= 1 && steps.v_steps <= 1))
                    continue;
                uv_faces.push_back({i, steps.u_steps, steps.v_steps});
                use_uv_tess[i] = true;
            }
        }
        if (!uv_faces.empty()) {
            cdt_faces_attempted = uv_faces.size();
            spdlog::info("  UV-grid CDT tessellation for {} regularly-curved faces ...",
                         uv_faces.size());
            for (auto& f : uv_faces) {
                auto res = uv_grid_retessellate(meshes[f.idx], segments[f.idx],
                                                f.u_steps, f.v_steps, min_edge_length, f.idx,
                                                dump_cdt_dir);
                if (res != UvTessResult::Ok) {
                    use_uv_tess[f.idx] = false;
                    if (res == UvTessResult::Failed)
                        spdlog::warn("  seg {}: uv_grid_retessellate failed, falling back to "
                                     "isotropic remeshing", f.idx);
                    else
                        spdlog::debug("  seg {}: UV-grid CDT not applicable, using BRepMesh "
                                      "interior + isotropic remeshing", f.idx);
                }
            }
            cdt_faces_succeeded = std::count(use_uv_tess.begin(), use_uv_tess.end(), true);
            spdlog::info("  {} segments tessellated via UV-grid CDT", cdt_faces_succeeded);

#ifdef PRECISION_MESH_HAS_VIEWER
            push_to_viewer(viewer_ptr, meshes, segments, original_faces, component_map, is_step,
                           "After UV-grid CDT", max_surface_error, free_edge_lines,
                           healed_edge_lines);
#endif
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
#ifdef PRECISION_MESH_HAS_VIEWER
        if (viewer_ptr) {
            rparams.on_iteration_done = [&](int iter) {
                push_to_viewer(viewer_ptr, meshes, segments, original_faces, component_map, is_step,
                               "Iteration " + std::to_string(iter) +
                               "/" + std::to_string(iterations),
                               max_surface_error, free_edge_lines, healed_edge_lines);
            };
        }
#endif
        remesh_and_project(meshes, segments, wire_projectors, surface_projectors, border_projectors,
                           rparams, use_uv_tess);
    } else {
        spdlog::info("  skipping remeshing (--iterations 0).");
    }

#ifdef PRECISION_MESH_HAS_VIEWER
    push_to_viewer(viewer_ptr, meshes, segments, original_faces, component_map, is_step,
                   "Final Result", max_surface_error, free_edge_lines, healed_edge_lines);
    if (viewer_ptr) viewer_ptr->notify_done();
#endif

    if (validate && is_step) {
        spdlog::info("validating tessellation (this may be slow) ...");
        double vtol = std::max(min_edge_length * 0.1, 1e-9);
        // Each segment validates against the REAL edges of the ORIGINAL face it came from
        // (so seams and subdivision cuts aren't treated as boundaries).
        std::vector<TopoDS_Face> edge_faces(segments.size());
        for (size_t i = 0; i < segments.size(); i++) {
            int orig = component_map.count(i) ? component_map.at(i) : (int)i;
            edge_faces[i] = (orig >= 0 && orig < (int)original_faces.size())
                            ? original_faces[orig] : segments[i];
        }
        auto vr = validate_tessellation(meshes, segments, edge_faces, vtol,
                                        validate_surface_error ? 4 : 0);
        spdlog::info("  validation ({} segments, {} triangles, tol={:.4g} {}):",
                     vr.segments, vr.total_tris, vtol * conversion_scale, output_unit);
        spdlog::info("    vertices: {} on-edge, {} interior", vr.border_verts, vr.interior_verts);
        spdlog::info("    on-edge vertex -> STEP edge:    max={:.4g} {} ({} beyond tol)",
                     vr.max_border_edge_dist * conversion_scale, output_unit, vr.misclassified_border);
        spdlog::info("    other vertex   -> STEP surface: max={:.4g} {} ({} beyond tol)",
                     vr.max_interior_face_dist * conversion_scale, output_unit, vr.misclassified_interior);
        if (validate_surface_error)
            spdlog::info("    surface error (tri samples):  max={:.4g} mean={:.4g} {} ({} samples)",
                         vr.max_surface_error * conversion_scale, vr.mean_surface_error * conversion_scale,
                         output_unit, vr.surface_samples);
        else
            spdlog::info("    surface error (tri samples):  not sampled (--validate-surface-error)");
        if (vr.open_boundary_edges > 0)
            spdlog::warn("    watertight: {} open boundary edges (triangle soup, by position)", vr.open_boundary_edges);
        else
            spdlog::info("    watertight: 0 open boundary edges");

        if (!validate_report.empty()) {
            std::ofstream report(validate_report);
            if (!report) {
                spdlog::error("Failed to open validation report file: {}", validate_report);
            }
            else {
                auto now = std::chrono::system_clock::now();
                std::time_t tt = std::chrono::system_clock::to_time_t(now);
                char timestamp[32] = {0};
                std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%dT%H:%M:%S",
                              std::localtime(&tt));
                double runtime_s = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - run_start).count();
                // All distances in output units (same conversion as the console block);
                // min/max edge length etc. were converted to model units earlier, so
                // scale them back.
                const double cs = conversion_scale;
                report << "# precision_mesh validation report\n";
                report << fmt::format("generated: \"{}\"\n", timestamp);
                report << fmt::format("input: \"{}\"\n", input);
                report << fmt::format("component: \"{}\"\n", selected_component->qualified_name);
                report << fmt::format("units: \"{}\"\n", output_unit);
                report << fmt::format("runtime_seconds: {:.2f}\n", runtime_s);
                report << "parameters:\n";
                report << fmt::format("  iterations: {}\n", iterations);
                report << fmt::format("  crease_angle_deg: {:.6g}\n", crease_angle);
                report << fmt::format("  min_edge_length: {:.6g}\n", min_edge_length * cs);
                report << fmt::format("  max_edge_length: {:.6g}\n", max_edge_length * cs);
                report << fmt::format("  max_surface_error: {:.6g}\n", max_surface_error * cs);
                report << fmt::format("  max_boundary_surface_error: {:.6g}\n",
                                      max_boundary_surface_error * cs);
                report << fmt::format("  no_subdivision: {}\n", no_subdivision);
                report << fmt::format("  no_projection: {}\n", no_projection);
                report << fmt::format("  no_uv_tess: {}\n", no_uv_tess);
                report << fmt::format("  no_tess_repair: {}\n", no_tess_repair);
                report << fmt::format("  raw_step_mesh: {}\n", raw_step_mesh);
                report << fmt::format("  unfreeze_step_boundaries: {}\n", unfreeze_step_boundaries);
                report << fmt::format("  validate_surface_error: {}\n", validate_surface_error);
                report << "tessellation:\n";
                report << fmt::format("  segments: {}\n", vr.segments);
                report << fmt::format("  triangles: {}\n", vr.total_tris);
                report << fmt::format("  cdt_faces_attempted: {}\n", cdt_faces_attempted);
                report << fmt::format("  cdt_faces_succeeded: {}\n", cdt_faces_succeeded);
                report << fmt::format("  cdt_faces_fallback: {}\n",
                                      cdt_faces_attempted - cdt_faces_succeeded);
                report << "validation:\n";
                report << fmt::format("  tolerance: {:.6g}\n", vtol * cs);
                report << fmt::format("  border_vertices: {}\n", vr.border_verts);
                report << fmt::format("  interior_vertices: {}\n", vr.interior_verts);
                report << fmt::format("  max_border_edge_dist: {:.6g}\n",
                                      vr.max_border_edge_dist * cs);
                report << fmt::format("  border_beyond_tol: {}\n", vr.misclassified_border);
                report << fmt::format("  max_interior_face_dist: {:.6g}\n",
                                      vr.max_interior_face_dist * cs);
                report << fmt::format("  interior_beyond_tol: {}\n", vr.misclassified_interior);
                report << fmt::format("  open_boundary_edges: {}\n", vr.open_boundary_edges);
                report << fmt::format("  non_manifold_edges: {}\n", vr.non_manifold_edges);
                report << fmt::format("  watertight: {}\n", vr.open_boundary_edges == 0);
                if (validate_surface_error) {
                    report << "  surface_error:\n";
                    report << fmt::format("    max: {:.6g}\n", vr.max_surface_error * cs);
                    report << fmt::format("    mean: {:.6g}\n", vr.mean_surface_error * cs);
                    report << fmt::format("    samples: {}\n", vr.surface_samples);
                }
                else {
                    report << "  surface_error: null\n";
                }
                spdlog::info("  validation report written to: {}", validate_report);
            }
        }
    }
    else if (!validate_report.empty()) {
        spdlog::warn("--validate-report requires a STEP input; no report written.");
    }

    return saveOutput(outputs, meshes, original_faces, component_map, conversion_scale) ? 0 : 1;

    }; // end do_pipeline lambda

#ifdef PRECISION_MESH_HAS_VIEWER
    if (enable_display) {
        Viewer viewer_obj(1280, 800, "PrecisionMesh");
        viewer_ptr = &viewer_obj;
        int result = 0;
        std::thread worker([&] {
            result = do_pipeline();
            // notify_done is called inside do_pipeline when viewer_ptr is set
        });
        viewer_obj.run_event_loop();
        worker.join();
        return result;
    }
#endif
    return do_pipeline();
}
