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
#include <precision_mesh/step_face_util.h>
#include <precision_mesh/step_projection.h>
#include <precision_mesh/step_subdivision.h>
#include <precision_mesh/step_tessellation.h>
#include <precision_mesh/stl.h>

namespace PMP = CGAL::Polygon_mesh_processing;

// property_map() returns std::optional in released CGAL 6.x but std::pair in 5.x and 6.0-dev.
// Detect via the return type rather than version number.
namespace {
    template <typename T, typename = void>
    struct is_optional_result : std::false_type {};
    template <typename T>
    struct is_optional_result<T, std::void_t<decltype(std::declval<T>().value())>> : std::true_type {};
}

template <typename Index, typename Value, typename SurfaceMesh>
auto lookup_property_map(SurfaceMesh& mesh, const std::string& name) {
    auto result = mesh.template property_map<Index, Value>(name);
    if constexpr (is_optional_result<decltype(result)>::value) {
        return result.value();
    } else {
        return result.first;
    }
}

std::unordered_set<std::string> mesh_formats = {".obj", ".off", ".ply", ".stl", ".ts", ".vtp"};
std::unordered_set<std::string> output_units = {"mm", "cm", "m", "in", "ft", "yd"};
std::unordered_map<std::string, float> to_meters = {
    {"mm", 0.001},
    {"cm", 0.01},
    {"m", 1.0},
    {"in", 0.0254},
    {"ft", 0.3048},
    {"yd", 0.9144}
};

struct Component {
    Component(TDF_Label label, TDF_Label reference, const std::string& name,
              const std::string& qualified_name, size_t index, const std::string& reference_name) :
        label(label), reference(reference), name(name), qualified_name(qualified_name),
        index(index), reference_name(reference_name) {}

    TDF_Label label;
    TDF_Label reference;
    std::string name;
    std::string qualified_name;
    int32_t index;
    std::string reference_name;
    TopoDS_Shape shape;
    std::string id;
    double surface_area=0.0;
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

gp_Trsf GetTransform(Handle(XCAFDoc_ShapeTool)& assembly, const TDF_Label& label, int depth = 0) {
    auto transformation = assembly->GetLocation(label);
    TDF_Label parent = label.Father();
    if (!parent.IsNull()) {
        if (depth >= 100) {
            spdlog::warn("Assembly nesting depth limit reached while computing transform.");
            return transformation;
        }
        auto parentTransformation = GetTransform(assembly, parent, depth + 1);
        transformation = transformation.Multiplied(parentTransformation);
    }
    return transformation;
}

void collectComponents(Handle(XCAFDoc_ShapeTool)& assembly,
                       const TDF_LabelSequence& labels,
                       const std::string& parent_qualified_name,
                       std::vector<std::shared_ptr<Component>>& components,
                       int depth = 0) {
    if (depth >= 100) {
        spdlog::warn("Assembly nesting depth limit reached at '{}'.", parent_qualified_name);
        return;
    }
    std::map<std::string, size_t> counts;
    for (Standard_Integer i = 1; i <= labels.Length(); i++) {
        TDF_Label label = labels.Value(i);

        std::string name = getName(label);
        std::string reference_name;

        size_t index = ++counts[name];
        std::string qualified_name = parent_qualified_name.empty()
            ? name
            : parent_qualified_name + "/" + name;
        if (index > 1) {
            qualified_name += "_" + std::to_string(index);
        }

        TDF_Label reference = label;
        if (assembly->GetReferredShape(label, reference)) {
            reference_name = getName(reference);
        }

        auto component = std::make_shared<Component>(label, reference, name, qualified_name,
                                                     index, reference_name);

        if (XCAFDoc_ShapeTool::IsSimpleShape(reference)) {
            gp_Trsf transform = GetTransform(assembly, label);

            component->shape = XCAFDoc_ShapeTool::GetShape(reference);
            if (!component->shape.IsNull()) {
                component->shape = component->shape.Located(TopLoc_Location(transform));
                components.push_back(component);
            }
        }

        TDF_LabelSequence children;
        assembly->GetComponents(label, children);
        collectComponents(assembly, children, qualified_name, components, depth + 1);
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

float getUnitConversionScale(const std::string& from, const std::string& to) {
    auto from_it = to_meters.find(from);
    auto to_it = to_meters.find(to);

    if (from_it == to_meters.end() || to_it == to_meters.end()) {
        return 1.0;
    }

    return from_it->second / to_it->second;
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
            if (output_units.count(normalized_output_unit) == 0) {
                spdlog::error("Output unit {} is not supported.", output_unit);
                return 1;
            }
            output_unit = normalized_output_unit;
            spdlog::info("  output unit: {} ", output_unit);
        }

        if (unit != output_unit) {
            if (output_units.count(unit) == 0) {
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
            double max_surface_error_auto = find_surface_error_param<Mesh>(selected_component->shape,
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

        auto tessellation = tessellate_shape<Mesh>(selected_component->shape, max_surface_error);
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

    spdlog::info("  splitting long border edges ...");
    // TODO(malban): replace split_long_edges on border edges with STEP-curve-based splitting.
    // The current approach uses PMP::split_long_edges which places new vertices by linear
    // interpolation (midpoint of the 3D chord).  The preferred approach is:
    //   1. For each STEP edge shared between adjacent sub-face meshes, compute arc-length-uniform
    //      split points directly on the parametric curve via BRep_Tool::Curve().
    //   2. Insert those vertices into each adjacent mesh using CGAL::Euler::split_edge followed
    //      by explicit repositioning (mesh.point(new_v) = curve_point).
    //   3. Walk each mesh's border halfedges to identify which correspond to a given STEP edge
    //      by proximity (reusing logic similar to get_border_vertex_projector_map).
    // Benefits over the current approach:
    //   - Split vertices land exactly on the STEP curve (no off-surface interpolation error).
    //   - Both adjacent meshes receive identical split points, guaranteeing topological
    //     consistency along shared borders without relying on OCCT's tessellation guarantee.
    size_t border_num_before = 0;
    size_t border_num_after = 0;

    for (size_t i = 0; i < meshes.size(); i++) {
        auto& mesh = meshes[i];
        std::vector<EdgeDescriptor> border_edges;
        PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(
            HalfEdge2Edge(mesh, border_edges)));
        border_num_before += border_edges.size();
        PMP::split_long_edges(border_edges, max_edge_length, mesh);
        // TODO(malban): project newly split border vertices onto STEP boundary curves.
        // split_long_edges places midpoints by linear interpolation, leaving them
        // slightly off the STEP curve.  Snapping them to the wire here (before the
        // remeshing loop) would give isotropic_remeshing accurate border constraints
        // from the first iteration, improving the adaptive sizing field and preventing
        // the interpolation error from being smoothed into adjacent interior faces.
        // Requires initializing wire_projectors before this loop.
        border_edges.clear();
        PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(
            HalfEdge2Edge(mesh, border_edges)));

        border_num_after += border_edges.size();
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

    WireProjectorCachePtr<Mesh> wire_projectors;
    std::vector<std::unique_ptr<StepProjector<Mesh>>> surface_projectors;
    std::vector<std::unique_ptr<StepBorderProjector<Mesh>>> border_projectors;
    if (is_step) {
        spdlog::info("    creating edge projectors ...");
        wire_projectors = get_edge_vertex_wire_projectors<Mesh>(selected_component->shape);

        if (!no_projection) {
            spdlog::info("    initializing surface projectors ...");
            for (const auto& segment : segments) {
                surface_projectors.push_back(std::make_unique<StepProjector<Mesh>>(segment));
                border_projectors.push_back(std::make_unique<StepBorderProjector<Mesh>>(segment));
            }
        }
    }

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
                    try {
                        if (is_step) {
                            PMP::isotropic_remeshing(faces(mesh), sizing_field, mesh,
                                CGAL::parameters::number_of_iterations(1)
                                                 .number_of_relaxation_steps(3)
                                                 .protect_constraints(true));
                        }
                        else {
                            auto crease_map = lookup_property_map<Mesh::Edge_index, bool>(mesh, "crease");
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

        if (is_step && !no_projection) {
            spdlog::info("      projecting ...");
            // Weight increases each iteration (1/N, 1/(N-1), ..., 1/1), reaching
            // full projection (weight=1) on the final iteration so vertices end up
            // exactly on the parametric surface.  The gradual schedule is intentional:
            // aggressively snapping vertices to the surface in early iterations can
            // produce degenerate or inverted faces before the mesh topology has had
            // time to adapt, so we blend toward the surface incrementally.
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, meshes.size()), [&](const tbb::blocked_range<size_t>& r) {
                    for (size_t m=r.begin(); m!=r.end(); ++m) {
                        double weight = 1.0 / (iterations - i);
                        project_to_step<Mesh>(segments[m], meshes[m], wire_projectors,
                                              *surface_projectors[m], *border_projectors[m], weight);
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
    return saveOutput(outputs, meshes, original_faces, component_map, conversion_scale) ? 0 : 1;


    return 0;
}
