project('precision_mesh', 'cpp',
  version         : '0.1',
  license         : 'GPLv3',
  default_options : ['cpp_std=c++17', 'warning_level=3'])
cc = meson.get_compiler('cpp')

cmake = import('cmake')

boost_dep = dependency('boost')
fontconfig_dep = dependency('fontconfig')
freetype_dep = dependency('freetype2')
gmp_dep = dependency('gmp', version : '>=6.0.0', required : true)
mpfr_dep = dependency('mpfr', version : '>=4.0.0', required : true)
tbb_dep = dependency('tbb')

cgal_opts = cmake.subproject_options()
cgal_opts.add_cmake_defines({
  'CMAKE_BUILD_TYPE': 'RelWithDebInfo',
  'WITH_examples': false,
  'WITH_tests': false,
  'WITH_demos': false,
  'WITH_doc': false,
  'WITH_CGAL_Qt5': false,
  'WITH_CGAL_ImageIO': false})
cgal_proj = cmake.subproject('cgal', options: cgal_opts)
cgal_modules = ['AABB_tree', 'Algebraic_foundations', 'Arithmetic_kernel', 'BGL',
                'Cartesian_kernel', 'Circulator', 'CGAL_Core', 'Distance_2', 'Distance_3',
                'Filtered_kernel', 'Generator', 'HalfedgeDS', 'Hash_map', 'Homogeneous_kernel',
                'Installation', 'Intersections_2', 'Intersections_3', 'Interval_support',
                'Kernel_23', 'Kernel_d', 'Modifier', 'Modular_arithmetic', 'Number_types',
                'Polygon_mesh_processing', 'Polyhedron', 'Profiling_tools', 'Property_map', 
                'Random_numbers', 'Spatial_searching', 'Surface_mesh', 'STL_Extension',
                'Stream_support', 'Union_find']
cgal_includes = []
foreach module : cgal_modules
    cgal_includes += 'subprojects/cgal/' + module + '/include'
endforeach

occt_opts = cmake.subproject_options()
occt_opts.add_cmake_defines({
  'CMAKE_BUILD_TYPE': 'RelWithDebInfo',
  'BUILD_LIBRARY_TYPE': 'Static', 
  'BUILD_MODULE_ApplicationFramework': false,
  'BUILD_MODULE_DETools': false,
  'BUILD_MODULE_Draw': false, 
  'BUILD_MODULE_Visualization': false})
occt_proj = cmake.subproject('occt', options: occt_opts)
occt_targets = ['TKBO', 'TKBRep', 'TKCAF', 'TKCDF', 'TKG2d', 'TKG3d', 'TKGeomAlgo', 'TKGeomBase',
                'TKernel', 'TKHLR', 'TKLCAF', 'TKMath', 'TKMesh', 'TKService', 'TKShHealing',
                'TKSTEP', 'TKSTEP209', 'TKSTEPAttr', 'TKSTEPBase', 'TKTopAlgo', 'TKV3d', 'TKVCAF',
                'TKXCAF', 'TKXDESTEP', 'TKXSBase']
occt_deps = []
foreach target : occt_targets
    occt_deps += occt_proj.dependency(target)
endforeach

eigen_proj = subproject('eigen')
eigen_dep = eigen_proj.get_variable('eigen_dep')

cli11_proj = subproject('cli11')
cli11_dep = cli11_proj.get_variable('CLI11_dep')

spdlog_proj = subproject('spdlog')
spdlog_dep = spdlog_proj.get_variable('spdlog_dep')

tinyply_proj = subproject('tinyply')
tinyply_dep = tinyply_proj.get_variable('tinyply_dep')

include_dir = include_directories('include')

precision_mesh = executable('precision_mesh', 'src/precision_mesh.cpp',
  cpp_args : ['-DCGAL_DISABLE_ROUNDING_MATH_CHECK'],
  dependencies : [boost_dep, cli11_dep, eigen_dep, gmp_dep, fontconfig_dep, freetype_dep, 
                  mpfr_dep, spdlog_dep, tbb_dep, tinyply_dep] + occt_deps,
  include_directories: [include_dir, cgal_includes],
  install : true)

test('basic', precision_mesh)
