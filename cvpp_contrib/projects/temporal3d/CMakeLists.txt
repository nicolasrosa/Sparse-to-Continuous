cmake_minimum_required( VERSION 2.8 )
project( regression )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake )
include( CVPP )

#add_executable( cnn
#        src/cnn/main_cnn.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/objects/object.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
#)

#add_executable( cnn_single
#        src/cnn/main_cnn_single.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/objects/object.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
#)

# add_executable( cnn_olorin_carina
#         src/cnn_hilbert/main_cnn_olorin_carina.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/objects/object.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
# )

# add_executable( cnn_olorin_kitti
#         src/cnn_hilbert/main_cnn_olorin_kitti.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/objects/object.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
# )

# add_executable( cnn_olorin_kitti_nick
#         src/cnn_hilbert/main_cnn_olorin_kitti_nick.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/objects/object.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
# )

add_executable( cnn_olorin_kitti_nick_fix_missing_images
        src/cnn_hilbert/main_cnn_olorin_kitti_nick_fix_missing_images.cpp
        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
        ${CPP_INCLUDE_DIR}/cvpp/objects/object.cpp
        ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
        ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
)

# add_executable( cnn_hilbert
#         src/cnn_hilbert/main_cnn_olorin.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/objects/object.cpp
#         ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
# )

#add_executable( cnn_comp
#        src/cnn/main_cnn_comp.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/objects/object.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
#)

#add_executable( cnn_comp2
#        src/cnn/main_cnn_comp2.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/objects/object.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
#)

#add_executable( temporal3d
#        src/update/main_temporal3d.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#)

#add_executable( comparison
#        src/update/main_comparison.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#)

#add_executable( vector
#        src/main_vector.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#)

#add_executable( process
#        src/main_process.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#)

#add_executable( incremental
#        src/main_incremental.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#)

#add_executable( motion
#        src/main_motion.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#)

#add_executable( marching
#        src/main_marching.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#)

#add_executable( marching2
#        src/main_marching2.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/interfaces/cpplot.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/marching_cubes/marching3D.cpp
#        ${CPP_INCLUDE_DIR}/cvpp/algorithms/hilbert_maps/real/hm_real.cpp
#)
