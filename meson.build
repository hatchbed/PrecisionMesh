project('precision_mesh', 'cpp',
  version         : '0.1',
  license         : 'GPLv3',
  default_options : ['cpp_std=c++17', 'warning_level=3', 'buildtype=debugoptimized'])

cmake = import('cmake')

_buildtype_to_cmake = {
  'plain':          'Release',
  'debug':          'Debug',
  'debugoptimized': 'RelWithDebInfo',
  'release':        'Release',
  'minsize':        'MinSizeRel',
}
_cmake_build_type = _buildtype_to_cmake.get(get_option('buildtype'), 'RelWithDebInfo')

boost_dep = dependency('boost', version: '>=1.66')
fontconfig_dep = dependency('fontconfig')
freetype_dep = dependency('freetype2')
gmp_dep = dependency('gmp', version : '>=6.0.0', required : true)
mpfr_dep = dependency('mpfr', version : '>=4.0.0', required : true)
tbb_dep = dependency('tbb')

cgal_modules = ['AABB_tree', 'Algebraic_foundations', 'Arithmetic_kernel', 'BGL',
                'Cartesian_kernel', 'Circulator', 'CGAL_Core', 'Distance_2', 'Distance_3',
                'Filtered_kernel', 'Generator', 'HalfedgeDS', 'Hash_map', 'Homogeneous_kernel',
                'Installation', 'Intersections_2', 'Intersections_3', 'Interval_support',
                'Kernel_23', 'Kernel_d', 'Modifier', 'Modular_arithmetic', 'Number_types',
                'PMP_Mesh_repair', 'PMP_Remeshing', 'Polygon_mesh_processing',
                'Polyhedron', 'Profiling_tools', 'Property_map',
                'Random_numbers', 'Spatial_searching', 'Spatial_sorting', 'Surface_mesh',
                'STL_Extension', 'Stream_support', 'TDS_2', 'Triangulation_2', 'Union_find']

cpp = meson.get_compiler('cpp')
cgal_opts = cmake.subproject_options()
cgal_opts.add_cmake_defines({
  'CMAKE_BUILD_TYPE': _cmake_build_type,
  'WITH_examples': false,
  'WITH_tests': false,
  'WITH_demos': false,
  'WITH_doc': false,
  'WITH_CGAL_Qt5': false,
  'WITH_CGAL_ImageIO': false})
cmake.subproject('cgal', options: cgal_opts)
_cgal_include_dirs = []
foreach module : cgal_modules
  _cgal_include_dirs += include_directories('subprojects/cgal/' + module + '/include')
endforeach
cgal_dep = declare_dependency(include_directories: _cgal_include_dirs)

occt_targets = ['TKBO', 'TKBRep', 'TKCAF', 'TKCDF', 'TKG2d', 'TKG3d', 'TKGeomAlgo', 'TKGeomBase',
                'TKernel', 'TKHLR', 'TKLCAF', 'TKMath', 'TKMesh', 'TKService', 'TKShHealing',
                'TKSTEP', 'TKSTEP209', 'TKSTEPAttr', 'TKSTEPBase', 'TKTopAlgo', 'TKV3d', 'TKVCAF',
                'TKXCAF', 'TKXDESTEP', 'TKXSBase']

# OCCT 7.8+ consolidated TKSTEP/TKSTEP209/TKSTEPAttr/TKSTEPBase/TKXDESTEP into TKDESTEP.
# The subproject (V7_7_2) still uses old names, so keep occt_targets unchanged for that path.
occt_system_targets = ['TKBO', 'TKBRep', 'TKCAF', 'TKCDF', 'TKG2d', 'TKG3d', 'TKGeomAlgo', 'TKGeomBase',
                       'TKernel', 'TKHLR', 'TKLCAF', 'TKMath', 'TKMesh', 'TKService', 'TKShHealing',
                       'TKDESTEP', 'TKTopAlgo', 'TKV3d', 'TKVCAF',
                       'TKXCAF', 'TKXSBase']

occt_dep = dependency('OpenCASCADE', method: 'cmake', version: '>=7.7.2', required: false,
                      modules: occt_system_targets)
if not occt_dep.found()
  occt_opts = cmake.subproject_options()
  occt_opts.add_cmake_defines({
    'CMAKE_BUILD_TYPE': _cmake_build_type,
    'BUILD_LIBRARY_TYPE': 'Shared',
    'BUILD_MODULE_ApplicationFramework': false,
    'BUILD_MODULE_DETools': false,
    'BUILD_MODULE_Draw': false,
    'BUILD_MODULE_Visualization': false})
  occt_proj = cmake.subproject('occt', options: occt_opts)
  _occt_deps = []
  foreach target : occt_targets
    _occt_deps += occt_proj.dependency(target)
  endforeach
  occt_dep = declare_dependency(dependencies: _occt_deps)
endif

eigen_dep = dependency('eigen3', fallback: ['eigen', 'eigen_dep'])

cli11_dep = dependency(['cli11', 'CLI11'], version: '>=2.4.1', fallback: ['cli11', 'CLI11_dep'])

spdlog_dep = dependency('spdlog', fallback: ['spdlog', 'spdlog_dep'])

tinyply_proj = subproject('tinyply')
tinyply_dep = tinyply_proj.get_variable('tinyply_dep')

# ── Optional viewer (GLFW + GLEW + ImGui + OpenGL) ────────────────────────────
viewer_option   = get_option('viewer')
glfw_dep        = dependency('glfw3', version: '>=3.3', required: viewer_option,
                             fallback: ['glfw', 'glfw_dep'])
glew_dep        = dependency('glew',  version: '>=2.1', required: viewer_option,
                             fallback: ['glew', 'glew_dep'])
gl_dep          = dependency('gl', required: viewer_option)
viewer_sources  = []
viewer_deps     = []
viewer_cpp_args = []

if glfw_dep.found() and glew_dep.found() and gl_dep.found()
  imgui_proj = subproject('imgui', default_options: ['glfw=enabled', 'opengl=enabled',
                                                      'vulkan=disabled', 'dx9=disabled',
                                                      'dx10=disabled', 'dx11=disabled',
                                                      'dx12=disabled', 'metal=disabled'])
  imgui_dep       = imgui_proj.get_variable('imgui_dep')
  viewer_sources  = ['src/viewer.cpp']
  viewer_deps     = [glfw_dep, glew_dep, gl_dep, imgui_dep]
  viewer_cpp_args = ['-DPRECISION_MESH_HAS_VIEWER']
endif
# ──────────────────────────────────────────────────────────────────────────────

include_dir = include_directories('include')

precision_mesh = executable('precision_mesh',
  'src/precision_mesh.cpp',
  'src/remeshing.cpp',
  'src/step_face_util.cpp',
  'src/step_projection.cpp',
  'src/step_reader.cpp',
  'src/step_subdivision.cpp',
  'src/step_tessellation.cpp',
  'src/unit_conversion.cpp',
  viewer_sources,
  cpp_args : ['-DCGAL_DISABLE_ROUNDING_MATH_CHECK'] + viewer_cpp_args,
  dependencies : [boost_dep, cgal_dep, cli11_dep, eigen_dep, gmp_dep, fontconfig_dep, freetype_dep,
                  mpfr_dep, occt_dep, spdlog_dep, tbb_dep, tinyply_dep] + viewer_deps,
  include_directories: [include_dir],
  install : true)

test('basic', precision_mesh)

# Install the .desktop file for application menus
install_data('data/precision_mesh.desktop',
  install_dir: get_option('datadir') / 'applications')

# Install the application icon
install_data('data/icon.png',
  rename: 'precision_mesh.png',
  install_dir: get_option('datadir') / 'icons/hicolor/32x32/apps')