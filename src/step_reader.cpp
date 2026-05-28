#include <precision_mesh/step_reader.h>

#include <map>
#include <string>

#include <spdlog/spdlog.h>

#include <TCollection_AsciiString.hxx>
#include <TDataStd_Name.hxx>
#include <TopLoc_Location.hxx>

std::string getName(const TDF_Label& label) {
    std::string name = "<unknown>";
    Handle(TDataStd_Name) name_attribute;
    if (label.FindAttribute(TDataStd_Name::GetID(), name_attribute)) {
        TCollection_AsciiString utf8String(name_attribute->Get().ToExtString(), Standard_False);
        name = utf8String.ToCString();
    }
    return name;
}

gp_Trsf GetTransform(Handle(XCAFDoc_ShapeTool)& assembly, const TDF_Label& label, int depth) {
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
                       int depth)
{
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
