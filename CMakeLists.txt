cmake_minimum_required(VERSION 3.5)

project (mongo-hash-table-profile)

include_directories (
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/../abseil-cpp"
    "${PROJECT_SOURCE_DIR}/../flat_hash_map"
    "${PROJECT_SOURCE_DIR}/../../lib/boost_1_68_0"
)

add_executable(mongo-hash-table-profile 
    main.cpp 
    "murmurhash3/MurmurHash3.cpp"
    "${PROJECT_SOURCE_DIR}/../abseil-cpp/absl/hash/internal/city.cc"
    "${PROJECT_SOURCE_DIR}/../abseil-cpp/absl/hash/internal/hash.cc"
    "${PROJECT_SOURCE_DIR}/../abseil-cpp/absl/container/internal/raw_hash_set.cc"
     )