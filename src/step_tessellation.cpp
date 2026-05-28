#include <precision_mesh/step_tessellation.h>

#include <spdlog/spdlog.h>

#include <BRep_Tool.hxx>
#include <Geom_Plane.hxx>
#include <STEPControl_Writer.hxx>

size_t get_short_edge_count(const std::vector<std::pair<Mesh, TopoDS_Face>>& tessellation,
                            double min_edge_length)
{
    size_t short_edges = 0;
    double min_sq = min_edge_length * min_edge_length;

    for (const auto& [mesh, face]: tessellation) {
        if (!Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(face)).IsNull()) {
            continue;
        }

        for (auto e: mesh.edges()) {
            if (!mesh.is_border(e)) {
                continue;
            }
            auto h = mesh.halfedge(e);
            auto p_start = mesh.point(mesh.source(h));
            auto p_end   = mesh.point(mesh.target(h));
            if (CGAL::squared_distance(p_start, p_end) < min_sq) {
                short_edges++;
            }
        }
    }

    return short_edges;
}

bool save_shape_as_step(const std::string& path, const TopoDS_Shape& shape) {
    STEPControl_Writer writer;
    auto status = writer.Transfer(shape, STEPControl_AsIs);
    if (status != IFSelect_RetDone) {
        spdlog::error("Error transferring shape to STEP writer.");
        return false;
    }

    status = writer.Write(path.c_str());
    if (status != IFSelect_RetDone) {
        spdlog::error("Error writing STEP file.");
    }

    return true;
}
