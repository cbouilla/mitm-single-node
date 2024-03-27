#ifndef MITM_UTIL_FILE_SYSTEM
#define MITM_UTIL_FILE_SYSTEM
#include <iostream>
#include <fstream>
#include <filesystem>

inline /* since it's defined in a header */
void create_folder_if_not_exist(const std::string& folder_path) {
  if (!std::filesystem::exists(folder_path)) {
    try {
      std::filesystem::create_directories(folder_path);
      std::cout << "Folder created: " << folder_path << std::endl;
    } catch (const std::filesystem::filesystem_error& e) {
      std::cerr << "Error creating folder: " << e.what() << std::endl;
    }
  }
}

inline /* since it's defined in a header */
void create_file_if_not_exist(const std::string& file_path) {
  if (!std::filesystem::exists(file_path)) {
    std::ofstream outfile(file_path);
    if (outfile.is_open()) {
      std::cout << "File created: " << file_path << std::endl;
      outfile.close();
    } else {
      std::cerr << "Error creating file: Unable to open file." << std::endl;
    }
  }
}
#endif 
