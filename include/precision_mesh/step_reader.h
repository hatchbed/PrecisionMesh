#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <gp_Trsf.hxx>
#include <TDF_Label.hxx>
#include <TDF_LabelSequence.hxx>
#include <TopoDS_Shape.hxx>
#include <XCAFDoc_ShapeTool.hxx>

// Represents a single component extracted from a STEP assembly, including its
// resolved shape (with transform applied) and hierarchical qualified name.
struct Component {
    Component(TDF_Label label, TDF_Label reference, const std::string& name,
              const std::string& qualified_name, size_t index,
              const std::string& reference_name)
        : label(label), reference(reference), name(name),
          qualified_name(qualified_name), index(index),
          reference_name(reference_name) {}

    TDF_Label label;
    TDF_Label reference;
    std::string name;
    std::string qualified_name;
    int32_t index;
    std::string reference_name;
    TopoDS_Shape shape;
    std::string id;
    double surface_area = 0.0;
};

// Returns the TDataStd_Name attribute of label, or "<unknown>" if absent.
std::string getName(const TDF_Label& label);

// Returns the cumulative transform for label by walking up the assembly hierarchy.
gp_Trsf GetTransform(Handle(XCAFDoc_ShapeTool)& assembly, const TDF_Label& label,
                     int depth = 0);

// Recursively collects all simple-shape components reachable from labels,
// applying accumulated transforms and building hierarchical qualified names.
void collectComponents(Handle(XCAFDoc_ShapeTool)& assembly,
                       const TDF_LabelSequence& labels,
                       const std::string& parent_qualified_name,
                       std::vector<std::shared_ptr<Component>>& components,
                       int depth = 0);
