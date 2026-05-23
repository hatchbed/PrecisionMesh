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
                'Polygon_mesh_processing', 'Polyhedron', 'Profiling_tools', 'Property_map',
                'Random_numbers', 'Spatial_searching', 'Surface_mesh', 'STL_Extension',
                'Stream_support', 'Union_find']

cpp = meson.get_compiler('cpp')
_cgal_ver = run_command('sh', '-c',
  'find /usr /usr/local -name version.h -path "*/CGAL/version.h" 2>/dev/null' +
  ' | head -1 | xargs grep -oP "(?<=#define CGAL_VERSION )\\S+" 2>/dev/null || echo ""',
  check: false).stdout().strip()
cgal_dep = dependency('', required: false)
if cpp.compiles('''
    #include <CGAL/version.h>
    #if CGAL_VERSION_MAJOR < 5 || (CGAL_VERSION_MAJOR == 5 && CGAL_VERSION_MINOR < 5)
    #error "CGAL version too old (need >= 5.5)"
    #endif
    int main() { return 0; }
  ''', name: 'CGAL ' + (_cgal_ver != '' ? _cgal_ver : '>= 5.5'))
  cgal_dep = declare_dependency()
endif
if not cgal_dep.found()
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
endif

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
    'BUILD_LIBRARY_TYPE': 'Static',
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

include_dir = include_directories('include')

precision_mesh = executable('precision_mesh', 'src/precision_mesh.cpp',
  cpp_args : ['-DCGAL_DISABLE_ROUNDING_MATH_CHECK'],
  dependencies : [boost_dep, cgal_dep, cli11_dep, eigen_dep, gmp_dep, fontconfig_dep, freetype_dep,
                  mpfr_dep, occt_dep, spdlog_dep, tbb_dep, tinyply_dep],
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